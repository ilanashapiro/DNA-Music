import requests
import sys
import ensembl_rest

def read_transcript_file(list):
    with open(list) as f:
        for line in f:
            yield line.strip()


def lookup_transcript(transcript):
    # server = "http://rest.ensembl.org"
    # ext = "/lookup/id/"+transcript+"?expand=1"

    # r = requests.get(server+ext, headers={"Content-Type": "application/json"})

    # if not r.ok:
    #     r.raise_for_status()
    #     sys.exit()

    # lookup = r.json()['Transcript'][0]

    data = ensembl_rest.symbol_lookup(
        species="Homo sapiens",
        symbol="TP53",
        params={'expand':True}
    )

    for entry in data['Transcript']:
        if entry['is_canonical']:
            return entry['Exon']
    return []


def intron_coords(exons, i):
    if exons[i]['strand'] == -1:
        intron_start = str(exons[i+1]['end']+1)
        intron_end = str(exons[i]['start']-1)
    else:
        intron_start = str(exons[i]['end']+1)
        intron_end = str(exons[i+1]['start']-1)

    coords = str(exons[i]['seq_region_name'])+":"+intron_start+".."+intron_end+":"+str(exons[i]['strand'])
    print(coords)
    return coords

def exon_coords(exon):
    exon_start = str(exon['start'])
    exon_end = str(exon['end'])

    coords = str(exon['seq_region_name'])+":"+exon_start+".."+exon_end+":"+str(exon['strand'])
    print(coords)
    return coords


if __name__ == "__main__":
    transcripts = [
        "ENSG00000141510"]

    for t in transcripts:
        exons = lookup_transcript(t)
        for exon in exons:
            print(exon)
        for i in range(0, len(exons)-1):
            icoords = intron_coords(exons, i)
            ecoords = exon_coords(exons[i])
            server = "http://rest.ensembl.org"
            iext = "/sequence/region/human/"+icoords
            eext = "/sequence/region/human/"+ecoords

            ir = requests.get(server+iext, headers={"Content-Type": "text/x-fasta"})

            if not ir.ok:
                print("IR NOT OK")
                ir.raise_for_status()
                sys.exit()

            er = requests.get(server+eext, headers={"Content-Type": "text/x-fasta"})

            if not er.ok:
                print("ER NOT OK")
                er.raise_for_status()
                sys.exit()
                
            # create a header like >ENST00000412061.3.Intron_1 chromosome:GRCh38:17:43094861:43095845:-1
            # use headline = ">"+t+".Intron_"+str(i+1) if you don't want the chromosomal position
            iheadline = ">"+t+".Intron_"+str(i+1)+" "+ir.text.split("\n")[0][1:]
            isequence = ir.text.split("\n", 1)[1]

            eheadline = ">"+t+".Exon_"+str(i+1)+" "+er.text.split("\n")[0][1:]
            esequence = er.text.split("\n", 1)[1]

            print(eheadline, esequence, sep="\n")
            print(iheadline, isequence, sep="\n")