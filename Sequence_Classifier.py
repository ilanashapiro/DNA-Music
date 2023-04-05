import math
import os, os.path

import Harmonic_Sequence
    
def num_files_specify_filename_ending(path, file_criterion):
    """ walks a whole directory structure
        and returns how many files are in it of a specified filename ending (which must be a STRING)"""
    AllFiles = list(os.walk(path))

    totalFiles = 0

    for item in AllFiles:
        (_, _, LoFiles) = item   
        for filename in LoFiles:
            if filename[-len(file_criterion):] == file_criterion:
                totalFiles += 1
    return totalFiles 

def get_sequence_percentages(root_directory):
    '''returns a dictionary matching a class name (e.g. "Mammals") to a list of length 4 consisting of the average percentages of the 4 harmonic 
    sequences contained in the DNA of the species/organisms in that class'''
    classes_dictionary = {}
    AllFiles = list(os.walk(root_directory))

    for item in AllFiles:
        (foldername, _, LoFiles) = item  
        isRoot = foldername == root_directory
        invalidDataDir = isRoot or any(invalidStr in foldername for invalidStr in ['.git', '-']) # things like ICCC-author-kit and the .git folder
        total_sequences = [0, 0, 0, 0]

        if num_files_specify_filename_ending(foldername, "translation_keys.txt") == 0 and not invalidDataDir:
            raise Exception("A class folder must have some harmonic sequence text files for organisms representing species within that class. Populate the class folder " + foldername + " and try again")

        if not invalidDataDir:
            for filename in LoFiles:
                if filename[-len("translation_keys.txt"):] == "translation_keys.txt":
                    fullfilename = os.path.join(foldername, filename)
                    total_sequences = [current_total + new_value for current_total, new_value in zip(total_sequences, Harmonic_Sequence.get_harmonic_sequences())]
            average_class_sequences =  [x/num_files_specify_filename_ending(foldername, "translation_keys.txt") for x in total_sequences]
            classes_dictionary[foldername] = average_class_sequences
    
    return classes_dictionary

def four_dim_distance(four_list1, four_list2):
    '''Find the "4-dimensional" distance between two lists, where each list is treated like a point in 4-dimensional space.
    Applies the distance formula in 4 dimensions: distance between (x1, y1, z1, w1) and (x2, y2, z2, w2) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 + (w2-w1)^2)'''
    if len(four_list1) != 4 or len(four_list2) != 4:
        raise Exception("Cannot compute 4-dimensional distance if length of each input vector is not 4!")
    
    sum_of_squared_distances = 0
    for i in range(4): #length of tuples should be 4
        sum_of_squared_distances += (four_list2[i] - four_list1[i]) ** 2
    return math.sqrt(sum_of_squared_distances)

def classify_species(species_harmseq_to_classify):
    '''classifies a given species into a class (e.g. "Mammals") based on the harmonic sequences derived from its DNA.
    The percentages of this species/organism’s DNA that correspond to the four harmonic sequences are determined and packaged into a list. 
    Then, using the “four dimensional distance formula,” the “distance” between the new organism’s harmonic sequence and the averages for each 
    class is computed. The organism is then classified into the class whose average harmonic sequence percentages are “closest” to that of the organism’s.'''
    classes_sequences_dictionary = get_sequence_percentages(".")
    if len(classes_sequences_dictionary) == 0:
        raise Exception("Cannot classify a species to an empty class folder. Populate the class folder with harmonic sequence data and try again.")

    if not isinstance(species_harmseq_to_classify, list):
        raise Exception("Species harmonic sequence should be a list of length 4 consisting of the percentages of the 4 harmonic sequences contained in the DNA")
    
    index = 0
    min_distance = 0
    min_class = ""
    for class_name, class_sequences in classes_sequences_dictionary.items():
        if index == 0: #initialize minimum distance to any element in the dictionary
            min_distance = four_dim_distance(class_sequences, species_harmseq_to_classify)
            min_class = class_name
            index += 1
        else:
            current_distance = four_dim_distance(class_sequences, species_harmseq_to_classify)
            if current_distance < min_distance:
                min_distance = current_distance
                min_class = class_name
    return min_class[2:] #remove the root directory text ".\"