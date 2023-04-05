import os, os.path
import math

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from urllib.error import HTTPError
import time

import Ensembl, DNA_MIDI, Sequence_Classifier, Harmonic_Sequence

'''These are the 5 runner functions for this program that can be called as desired. If starting with an empty directory tree, 
you must first populate class directories with nucleotide files and their associated midi and musical keys .txt file by calling write_nucleotides_and_generate_midi_txt_files 
(with the desired parameters, see function) for at least 1 class and 1 organism/species.
Then, if you want to get the musicality score of an organism (get_species_musicality_score) or classify it into a class (classify_species), populate the root 
directory with at least 1 nucleotide .txt file for this organism and its associated midi and musical keys .txt file by calling write_nucleotides_and_generate_midi_txt_files 
again, but this time having the root directory "." as your folder. 
    Alternatively, if you want, you can create an instance of the NCBI class on your own through terminal and call its write_sequences_file_to_folder function if you just want to 
generate nucleotide .txt files.
    Then, if you want to do it this way, you will need to populate the directory in the same way described above (including the files in the root directory to be used for classification 
and musicality score analysis) with nucleotide .txt files. Then, you can either call the generate_single_midi_and_txt for each nucleotide .txt file in order to generate the MIDI and
associated musical keys .txt file, or you can call generate_midi_and_associated_txt_files, which generates the MIDI and musical keys .txt files for the nucleotide .txt files of the entire
directory tree.
    Then, in the same way as described above, classify organisms in the root directory into classes or get musicality scores by calling classify_species or get_species_musicality_score.
'''

def generate_single_midi_and_txt(folderpath, filename, num_tracks, track_num):
    '''(used as helper function for write_nucleotides_and_generate_midi_txt_files. However, if you on your own create nucleotide .txt files through an instance of the NCBI class, 
    you can generate the MIDI and musical keys .txt files on your own through this method. ASSUMES YOU ALREADAY HAVE NUCLEOTIDE .TXT FILES GENERATED. IF NOT USE write_nucleotides_and_generate_midi_txt_files. 
    Description:
    Generate the MIDI and .txt file of associated musical key changes for a given filename. This filename should represent a .txt file of 
    nucleotides for a desired organism/species through the DNA_to_MIDI class.'''
    fullfilename = os.path.join(folderpath, filename)
    with open(fullfilename, 'r') as file:
        nucleotides = file.read().replace('\n', '')
        DNA_MIDI.create_track(folderpath, filename[:-4], track_num, nucleotides)
        DNA_MIDI.add_notes(track_num)

    DNA_MIDI.write_to_disk(fullfilename)

def write_nucleotides_and_generate_midi_txt_files(class_name, species_scientific_name, gene_name, class_folder_path = "."):
    '''runner function that takes in a class name, the scientific name for the desired species, the desired gene name, and optionally a class folder path 
    (if you don't want it to be created in the root directory. Creates an instance of the NCBI class, and then uses Biopython to access the NCBI (National Center for Biotechnology Information
    Then the user is asked to approve the DNA selection from the NCBI through terminal prompt, and once the data is approved, it is written to a .txt file with the species 
    scientific name to a folder with the class name. Then, the MIDI file and .txt file of musical key changes are generated for that DNA .txt file'''
    if(not os.path.isdir(os.path.join(class_folder_path, class_name))):
        os.mkdir(os.path.join(class_folder_path, class_name))
        
    Ensembl.write_sequences_file_to_folder(class_name, class_folder_path, species_scientific_name, gene_name)
    generate_single_midi_and_txt(os.path.join(class_folder_path, class_name), species_scientific_name + ".txt", 1, 0)

def classify_species(organism_scientific_name):
    '''runner function that classifies a given species into a class (e.g. "Mammals"). Creates an instance of the Harmonic_Sequences class to get the harmonic sequences list (percent of DNA corresponding to 
    each of the 4 harmonic sequences). Then creates an instance of the Sequence_Classifier class, and uses the classify_species method to assign a class to the given organism repsenting a specific species.'''
    harmseq_to_classify = Harmonic_Sequence.get_harmonic_sequences(organism_scientific_name + "_translation_keys.txt") # e.g. "delphinapterus leucas.txt"
    assigned_class = Sequence_Classifier.classify_species(harmseq_to_classify)
    return assigned_class

def get_species_musicality_score(organism_scientific_name):
    '''runner function that gives a species a musicality score on a scale of 0-100 (higher score means more musical) through the Harmonic_Sequences class.'''
    musicality_score = Harmonic_Sequence.get_musicality_rating(organism_scientific_name + "_translation_keys.txt")
    return musicality_score

def generate_midi_and_associated_txt_files():
    '''if you, on your own, create multiple nucleotide .txt files through an instance of the NCBI class, possibly in differnt folders,  
    you can generate the MIDI and musical keys .txt files for all of the nucleotide .txt files on your own through this method.
    ASSUMES YOU ALREADAY HAVE NUCLEOTIDE .TXT FILES GENERATED. IF NOT USE write_nucleotides_and_generate_midi_txt_files. Description:
    Generate the MIDI and .txt file of associated musical key changes for all the .txt files of nucleotides in the entire directory tree.'''
    AllFiles = list(os.walk("."))
    for item in AllFiles:
        (foldername, LoDirs, LoFiles) = item  
        isRoot = foldername == "."
        if not isRoot:
            for filename in LoFiles:
                isHarmonicKeys = filename[-len("translation_keys.txt"):] == "translation_keys.txt"
                if filename[-3:] == "txt" and not isHarmonicKeys:
                    generate_single_midi_and_txt(foldername, filename, 1, 0)

class_species_dict_TP53 = {
    "Amphibia":["Rhinatrema bivittatum", "Nanorana parkeri", "Rana temporaria", "Bufo gargarizans", "Geotrypetes seraphini", "Xenopus tropicalis", "Microcaecilia unicolor"],
    "Mammalia":["Bos taurus", "Canis lupus familiaris", "Delphinapterus leucas", "Homo sapiens", "Mus musculus", "Oryctolagus cuniculus", "Ovis aries"],
    "Actinopterygii":["Chanos chanos", "Danio rerio", "Gobiocypris rarus", "Oncorhynchus mykiss", "Salmo salar", "Oncorhynchus keta", "Chelmon rostratus"],
    "Lepidosauria":["Pogona vitticeps", "Varanus komodoensis", "Sceloporus undulatus", "Zootoca vivipara", "Podarcis muralis", "Gekko japonicus", "Crotalus tigris"],
    "Aves":["Gallus gallus", "Anas platyrhynchos", "Haliaeetus leucocephalus", "Lonchura striata domestica", "Taeniopygia guttata", "Pseudopodoces humilis", "Sturnus vulgaris"],
    "Testudines":["Dermochelys coriacea", "Caretta caretta", "Chelonia mydas", "Mauremys reevesii", "Mauremys mutica", "Terrapene carolina triunguis", "Chelonoidis abingdonii"]
}
write_nucleotides_and_generate_midi_txt_files("Mammalia", "Homo sapiens", "TP53")
num_classified_correctly = 0
total_species = 0
highest_musicality_score = 0
highest_musicality_species = ""
# for (class_name, species_list) in class_species_dict_TP53.items():
#     for species in species_list:
#         # write_nucleotides_and_generate_midi_txt_files(class_name, species, "TP53")
#         print("SPECIES: ", species)
#         classification = classify_species(class_name + "/" + species)
#         if classification == class_name:
#             num_classified_correctly += 1
#         print("CLASSIFICATION: ", classification)
#         print("CORRECT CLASSIFICATION: ", class_name)
#         musicality_score = get_species_musicality_score(class_name + "/" + species)
#         if musicality_score > highest_musicality_score:
#             highest_musicality_score = musicality_score
#             highest_musicality_species = species
#         print("MUSICALITY SCORE: ", get_species_musicality_score(class_name + "/" + species))
#         print()
#         total_species += 1

print("PERCENT CLASSIFIED CORRECTLY: ", num_classified_correctly / total_species * 100)
print("HIGHEST MUSICALITY SPECIES AND SCORE: ", highest_musicality_species, highest_musicality_score)



