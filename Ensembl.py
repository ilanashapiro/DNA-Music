import ensembl_rest
import requests
from requests.exceptions import HTTPError

def lookup_transcript(species_name, gene_name):
    data = ensembl_rest.symbol_lookup(
        species=species_name,
        symbol=gene_name,
        params={'expand':True}
    )

    for entry in data['Transcript']:
        if entry['is_canonical']:
            return (entry['Exon'], entry['id'])
    return []

def lookup_multiple_regions(regions):
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    try:
        response = requests.post(server+ext, headers=headers, data='{ "regions" : ' + regions + ' }')
        response.raise_for_status()
        decoded = response.json()

        regions = []
        for entry in decoded:
            regions.append(entry['seq'])

        return regions
    except HTTPError as exc:
        raise Exception(exc.response.status_code)
    
def get_intron_coords(exons, i):
    if exons[i]['strand'] == -1:
        intron_start = str(exons[i+1]['end']+1)
        intron_end = str(exons[i]['start']-1)
    else:
        intron_start = str(exons[i]['end']+1)
        intron_end = str(exons[i+1]['start']-1)

    coords = str(exons[i]['seq_region_name'])+":"+intron_start+".."+intron_end+":"+str(exons[i]['strand'])
    return coords

def get_exon_coords(exon):
    exon_start = str(exon['start'])
    exon_end = str(exon['end'])

    coords = str(exon['seq_region_name'])+":"+exon_start+".."+exon_end+":"+str(exon['strand'])
    return coords

# DEPRECATED -- now using optimized multi-region query
def lookup_single_region_from_coords(coords):
    server = "http://rest.ensembl.org"
    generic_ext = "/sequence/region/human/"

    try:
        response = requests.get(server+generic_ext+coords, headers={"Content-Type": "text/plain"})
        response.raise_for_status()
        seq = response.text.split("\n", 1)[0]
        return seq
    except HTTPError as exc:
        raise Exception(exc.response.status_code)

def lookup_5prime_UTR_from_transcript_id(ID):
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/" + ID + "?mask_feature=1;type=cdna"

    def get_lowercase_prefix(input_str):
        c = input_str[0]
        idx = 1
        res = ""
        while idx < len(input_str) and c.islower(): 
            res += c
            c = input_str[idx]
            idx += 1
        return res

    try:
        response = requests.get(server+ext, headers={"Content-Type": "text/plain"})
        response.raise_for_status()
        seq = response.text.split("\n", 1)[0]
        return get_lowercase_prefix(seq).upper()
    except HTTPError as exc:
        raise Exception(exc.response.status_code)

def get_sequence(species_name, gene_name):
    (exons_info_list, canonical_transcript_id) = lookup_transcript(species_name, gene_name)
    CDS_list = []
    UTR_5prime_exons_list = []
    UTR_5prime = lookup_5prime_UTR_from_transcript_id(canonical_transcript_id)

    exon_regions_coords = []
    intron_regions_coords = []
    for i in range(0, len(exons_info_list)):
        exon_regions_coords.append(get_exon_coords(exons_info_list[i]))
        if i < len(exons_info_list) - 1:
            intron_regions_coords.append(get_intron_coords(exons_info_list, i))
    
    exon_regions = lookup_multiple_regions(str(exon_regions_coords).replace("\'", "\"" ))
    intron_regions = lookup_multiple_regions(str(intron_regions_coords).replace("\'", "\"" ))

    for i in range(0, len(exon_regions) - 1):
        exon_seq = exon_regions[i]
        intron_seq = intron_regions[i]

        if UTR_5prime.startswith(exon_seq): # the entire exon is part of the 5' UTR
            UTR_5prime_split_on_curr_exon = UTR_5prime.split(exon_seq)
            UTR_5prime_exons_list.append(exon_seq)
            UTR_5prime = ''.join(UTR_5prime_split_on_curr_exon[1:])
            exon_seq = ""
        elif len(UTR_5prime) > 0 and exon_seq.startswith(UTR_5prime): # part of the exon constitutes the last part of the 5' UTR
            curr_exon_split_on_UTR_5prime = exon_seq.split(UTR_5prime)
            UTR_5prime_exons_list.append(UTR_5prime)
            exon_seq = ''.join(curr_exon_split_on_UTR_5prime[1:])
            UTR_5prime = ""

        # exon_seq could be emptry str if entire exon is UTR, which we just stripped, but we need it as placeholder so exon indexing in the list holds
        CDS_list.append(exon_seq)  
        CDS_list.append(intron_seq)
            
    # create a header like >ENST00000412061.3.Intron_1 chromosome:GRCh38:17:43094861:43095845:-1
    # use headline = ">"+t+".Intron_"+str(i+1) if you don't want the chromosomal position
    # print(headline, sequence, sep="\n")

    final_exon_coords = get_exon_coords(exons_info_list[-1])
    final_exon_seq = lookup_single_region_from_coords(final_exon_coords)
    CDS_list.append(final_exon_seq)

    return (UTR_5prime_exons_list, CDS_list)

# def write_sequences_file_to_folder(class_name, class_folder_path, species_scientific_name, gene_name):
#     '''write the DNA sequence for the desired gene and species into a folder that has the name of the species' biological class'''
#     sequence = query_ensembl_sequence(species_scientific_name, gene_name).text

#     if(not os.path.isdir(os.path.join(class_folder_path, class_name))):
#         os.mkdir(os.path.join(class_folder_path, class_name))
#     full_file_name = os.path.join(class_folder_path, class_name, species_scientific_name)  + ".txt"
    
#     with open(full_file_name, 'w') as f:
#         f.write(sequence)

# get_sequence("Homo sapiens", "TP53")
print("HERE")
for entry in get_sequence("Homo sapiens", "TP53"):
    print(entry, "\n")