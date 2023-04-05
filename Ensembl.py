import ensembl_rest
import os, os.path
import requests
from http import HTTPError

'''a module to obtain the cDNA sequence of a gene from Ensembl. Obtains FASTA-formatted text file exons (including 5' and 3' UTR) that are capitalized and introns in lower case'''
def query_ensembl_sequence(species_scientific_name, gene_name):
    '''formats the nucleotide sequence from the NCBI by removing the information header so we are just left with the nucleotides'''
    id = ensembl_rest.symbol_lookup(
        species=species_scientific_name,
        symbol=gene_name
    ).get('id')

    server = "https://rest.ensembl.org"
    ext = "/sequence/id/" + id + "?mask_feature=1" # mask_feature makes it so the exons are in caps and introns in lowercase

    try:
        response = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
        response.raise_for_status()
    except HTTPError as exc:
        raise Exception(exc.response.status_code)
    
def write_sequences_file_to_folder(class_name, class_folder_path, species_scientific_name, gene_name):
    '''using the ID list from the generate_ID_list function above, search for each ID using Entrez.efetch and show the user the informational header of each result.
    If the user types False, move on to the next ID in the list. If the user types True, write the formatted sequence to a .txt file with the organism's scientific name.
    If the user types False and the end of the ID list is reached, or if the ID list is originally empty (e.g. a bad request), and exception is raised.'''
    sequence = query_ensembl_sequence(species_scientific_name, gene_name).text

    if(not os.path.isdir(os.path.join(class_folder_path, class_name))):
        os.mkdir(os.path.join(class_folder_path, class_name))
    fullfilename = os.path.join(class_folder_path, class_name, species_scientific_name)  + ".txt"
    
    with open(fullfilename, 'w') as f:
        f.write(sequence.text)
