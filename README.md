# DNA to Music: Conversion and Analysis
Converts DNA to MIDI and analyzes harmonic sequences.
DNA_Music.py has the code for the project
The midi files are musical representations of the DNA of the organisms they are named for.
The filenames ending with "translation_keys.txt" are a record of the musical key changes that occur during translation regions in the DNA (when converted to RNA) for later use in harmonic analysis.

# Instructions on Running the Program
At the bottom of DNA_Music.py are the 5 runner functions for this program that can be called as desired. If starting with an empty directory tree, you must first populate class directories with nucleotide files and their associated midi and musical keys .txt file by calling write_nucleotides_and_generate_midi_txt_files (with the desired parameters, see function) for at least 1 class and 1 organism/species.

Then, if you want to get the musicality score of an organism (get_species_musicality_score) or classify it into a class (classify_species), populate the root directory with at least 1 nucleotide .txt file for this organism and its associated midi and musical keys .txt file by calling write_nucleotides_and_generate_midi_txt_files again, but this time having the root directory "." as your folder. 

Alternatively, if you want, you can create an instance of the NCBI class on your own through terminal and call its write_sequences_file_to_folder function if you just want to generate nucleotide .txt files.

Then, if you want to do it this way, you will need to populate the directory in the same way described above (including the files in the root directory to be used for classification and musicality score analysis) with nucleotide .txt files. Then, you can either call the generate_single_midi_and_txt for each nucleotide .txt file in order to generate the MIDI and associated musical keys .txt file, or you can call generate_midi_and_associated_txt_files, which generates the MIDI and musical keys .txt files for the nucleotide .txt files of the entire directory tree.

Then, in the same way as described above, classify organisms in the root directory into classes or get musicality scores by calling classify_species or get_species_musicality_score.
