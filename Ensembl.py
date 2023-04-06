import ensembl_rest
import os, os.path
import requests
from requests.exceptions import HTTPError

'''a module to obtain the DNA sequence (the coding strand) of a gene from Ensembl. Obtains FASTA-formatted text file exons (including 5' and 3' UTR) that are capitalized and introns in lower case'''
def query_ensembl_sequence(species_scientific_name, gene_name):
    '''queries the nucleotide sequence from Ensembl. Exons are in caps, introns lowercase. The sequence represents the coding strand of the DNA (i.e. matches RNA, 1-1 T-U valid substitution)'''
    # id = ensembl_rest.symbol_lookup(
    #     species=species_scientific_name,
    #     symbol=gene_name
    # ).get('id')

    data = ensembl_rest.symbol_lookup(
        species=species_scientific_name,
        symbol=gene_name,
        params={'expand':True}
    )
    
    for entry in data['Transcript']:
        if entry['is_canonical']:
            for exon in entry['Exon']:
                print(exon)
            break
    # server = "https://rest.ensembl.org"
    # ext = "/sequence/id/" + id + "?mask_feature=1" # mask_feature makes it so the exons are in caps and introns in lowercase

    # ext1 = "/lookup/id/ENST00000269305.9?object_type=transcript;expand=1"

    # print(server+ext)
    # try:
    #     response = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    #     response.raise_for_status()
    #     return response
    # except HTTPError as exc:
    #     raise Exception(exc.response.status_code)
    
def write_sequences_file_to_folder(class_name, class_folder_path, species_scientific_name, gene_name):
    '''write the DNA sequence for the desired gene and species into a folder that has the name of the species' biological class'''
    sequence = query_ensembl_sequence(species_scientific_name, gene_name).text

    if(not os.path.isdir(os.path.join(class_folder_path, class_name))):
        os.mkdir(os.path.join(class_folder_path, class_name))
    full_file_name = os.path.join(class_folder_path, class_name, species_scientific_name)  + ".txt"
    
    with open(full_file_name, 'w') as f:
        f.write(sequence)

query_ensembl_sequence("Homo sapiens", "TP53")