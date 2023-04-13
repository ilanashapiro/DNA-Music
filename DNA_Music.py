import os, os.path
import Ensembl, DNA_MIDI, Sequence_Classifier, Harmonic_Sequence

'''These are the 5 runner functions for this program that can be called as desired. If starting with an empty directory tree, 
you must first populate class directories with nucleotide files and their associated midi and musical keys .txt file by calling write_nucleotides_and_generate_midi_txt_files 
(with the desired parameters, see function) for at least 1 class and 1 organism/species.
Then, if you want to get the musicality score of an organism (get_species_musicality_score) or classify it into a class (classify_species), populate the root 
directory with at least 1 nucleotide .txt file for this organism and its associated midi and musical keys .txt file by calling write_nucleotides_and_generate_midi_txt_files 
again, but this time having the root directory "." as your folder. 
'''

def write_nucleotides_and_generate_midi_txt_files(class_name, species_scientific_name, gene_name, class_folder_path = "."):
    '''runner function that takes in a class name, the scientific name for the desired species, the desired gene name, and optionally a class folder path 
    (if you don't want it to be created in the root directory. Creates an instance of the NCBI class, and then uses Biopython to access the NCBI (National Center for Biotechnology Information
    Then the user is asked to approve the DNA selection from the NCBI through terminal prompt, and once the data is approved, it is written to a .txt file with the species 
    scientific name to a folder with the class name. Then, the MIDI file and .txt file of musical key changes are generated for that DNA .txt file'''
    folder_path = os.path.join(class_folder_path, class_name)
    full_file_name = os.path.join(folder_path, species_scientific_name + ".txt").replace(" ", "_" )
    if(not os.path.isdir(folder_path)):
        os.mkdir(folder_path)
        
    (UTR_5prime_exons_list, CDS_list) = Ensembl.get_sequence(species_scientific_name, gene_name)
    DNA_MIDI.add_notes(folder_path, species_scientific_name, full_file_name, UTR_5prime_exons_list, CDS_list)

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

class_species_dict_TP53 = {
    "Amphibia":["Rhinatrema bivittatum", "Nanorana parkeri", "Rana temporaria", "Bufo gargarizans", "Geotrypetes seraphini", "Xenopus tropicalis", "Microcaecilia unicolor"],
    "Mammalia":["Bos taurus", "Canis lupus familiaris", "Delphinapterus leucas", "Homo sapiens", "Mus musculus", "Oryctolagus cuniculus", "Ovis aries"],
    "Actinopterygii":["Chanos chanos", "Danio rerio", "Gobiocypris rarus", "Oncorhynchus mykiss", "Salmo salar", "Oncorhynchus keta", "Chelmon rostratus"],
    "Lepidosauria":["Pogona vitticeps", "Varanus komodoensis", "Sceloporus undulatus", "Zootoca vivipara", "Podarcis muralis", "Gekko japonicus", "Crotalus tigris"],
    "Aves":["Gallus gallus", "Anas platyrhynchos", "Haliaeetus leucocephalus", "Lonchura striata domestica", "Taeniopygia guttata", "Pseudopodoces humilis", "Sturnus vulgaris"],
    "Testudines":["Dermochelys coriacea", "Caretta caretta", "Chelonia mydas", "Mauremys reevesii", "Mauremys mutica", "Terrapene carolina triunguis", "Chelonoidis abingdonii"]
}
write_nucleotides_and_generate_midi_txt_files("Mammalia", "Mouse", "TP53")
num_classified_correctly = 0
total_species = 0
highest_musicality_score = 0
highest_musicality_species = ""
musicality_score = get_species_musicality_score("Mammalia" + "/" + "Homo sapiens")
print(musicality_score)
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

# print("PERCENT CLASSIFIED CORRECTLY: ", num_classified_correctly / total_species * 100)
# print("HIGHEST MUSICALITY SPECIES AND SCORE: ", highest_musicality_species, highest_musicality_score)



