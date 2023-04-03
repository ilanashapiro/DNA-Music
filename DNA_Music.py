import os, os.path
import math
import functools
from midiutil.MidiFile import MIDIFile
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

class NCBI:
    '''a class to encapsulate the abilities of the Biopython library to access and interpret information from the NCBI (National Institute for Biotechnology Information)'''
    def __init__(self, email = 'liiro.ana@gmail.com'):
        Entrez.email = email

    def generate_ID_list(self, organism_scientific_name, gene_name):
        '''generates an ID list from the NCBI for a given organism and gene search. Intended for later use with Biopython's Entrez.efetch to attain the genetic information for the given organism and gene '''
        handle = Entrez.esearch(db='nucleotide', term = [organism_scientific_name + "[Orgn] AND " + gene_name + "[Gene]"])
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]

    def get_formatted_sequence(self, sequence):
        '''formats the nucleotide sequence from the NCBI by removing the information header so we are just left with the nucleotides'''
        sequence = "\n".join(sequence.split("\n")[1:]) #remove the header
        sequence = sequence.replace("N", "").replace("\n", "") #remove nucleotides labeled "N", which means they are unconfirmed/unidentified by the NCBI website. Then remove all newline characters
        return sequence

    def get_sequence_header(self, sequence):
        '''strip the informational header from the nucleotide sequence from the NCBI and return it for informative purposes'''
        header = sequence.split("\n")[0] #get the header
        return header
        
    def write_sequences_file_to_folder(self, class_name, class_folder_path, species_scientific_name, gene_name):
        '''using the ID list from the generate_ID_list function above, search for each ID using Entrez.efetch and show the user the informational header of each result.
        If the user types False, move on to the next ID in the list. If the user types True, write the formatted sequence to a .txt file with the organism's scientific name.
        If the user types False and the end of the ID list is reached, or if the ID list is originally empty (e.g. a bad request), and exception is raised.'''
        ID_list = self.generate_ID_list(species_scientific_name, gene_name)

        if len(ID_list) == 0: #in case of bad request, in which case the ID list will be empty
            raise Exception("Bad request! Make sure the FULL scientific organism name and FULL gene name are correct")
    
        if(not os.path.isdir(os.path.join(class_folder_path, class_name))):
            os.mkdir(os.path.join(class_folder_path, class_name))
        fullfilename = os.path.join(class_folder_path, class_name, species_scientific_name)  + ".txt"
        file_written = False
        for seq_id in ID_list:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=seq_id)
            record = handle.read()
            write_file = (input("Is this the correct sequence? \n Type \"True\" or \"False.\" (case sensitive) \n" + self.get_sequence_header(record) + "\n") == "True")
            if write_file:
                with open(fullfilename, 'w') as f:
                    f.write(self.get_formatted_sequence(record))
                handle.close()
                file_written = True
                break

        if not file_written:
            raise Exception("Search for another gene that has the documentation you want in the NCBI.")
            
class DNA_to_MIDI:
    '''a class to convert DNA (in the form of a text file of nucleotides) into a MIDI file'''
    CHANNEL = 0 #set instrument to a piano

    #pitches correspondence according to MIDIUtil
    PITCH_DICTIONARY = {"C": 60,
                        "C#": 61,
                        "D": 62,
                        "Eb": 63,
                        "E": 64,
                        "F": 65,
                        "F#": 66,
                        "G": 67,
                        "Ab": 68,
                        "A": 69,
                        "Bb": 70,
                        "B": 71}
    KEYS = ["C", "C#", "D", "Eb", "E", "F", "F#", "G", "Ab", "A", "Bb", "B"] 

    DNA_TO_CHROMATIC = {"C": "",  
                        "G": "",
                        "A": "",
                        "U": ""}

    START_CODONS = ["AUG"]
    STOP_CODONS = ["UAA", "UGA", "UAG"]
    AMINO_ACIDS= {#essential
                  "AUG":"C", #methionine/start codon
                  "AUU":"C#",
                  "AUC":"C#",
                  "AUA":"C#", #isoleucine
                  "AAA":"D",
                  "AAG":"D", #lysine
                  "ACU":"Eb",
                  "ACC":"Eb",
                  "ACA":"Eb",
                  "ACG":"Eb", #threonine
                  "UUU":"E",
                  "UUC":"E", #phenylalanine
                  "UGG":"F", #tryptophan
                  "UUA":"F#",
                  "UUG":"F#",
                  "CUU":"F#",
                  "CUC":"F#",
                  "CUA":"F#",
                  "CUG":"F#", #leucine
                  "CAU":"G",
                  "CAC":"G", #histidine
                  "GUU":"Ab",
                  "GUC":"Ab",
                  "GUA":"Ab",
                  "GUG":"Ab", #valine

                  #nonessential
                  "AAU":"A",
                  "AAC":"A", #asparagine
                  "GAU":"Bb",
                  "GAC":"Bb", #aspartate 
                  "GCU":"B",
                  "GCC":"B",
                  "GCA":"B",
                  "GCG":"B", #alanine
                  
                  #conditionally essential
                  "UAU":"C",
                  "UAC":"C", #tyrosine
                  "UGU":"D",
                  "UGC":"D", #cysteine
                  "UCC":"E",
                  "UCU":"E",
                  "UCA":"E",
                  "UCG":"E", 
                  "AGU":"E",
                  "AGC":"E", #serine
                  "AGA":"F",
                  "AGG":"F", 
                  "CGU":"F",
                  "CGC":"F",
                  "CGA":"F",
                  "CGG":"F", #arginine
                  "CCU":"G",
                  "CCC":"G",
                  "CCA":"G",
                  "CCG":"G", #proline
                  "CAA":"A",
                  "CAG":"A",
                  "GAA":"A",
                  "GAG":"A", #glutamine/glutamic acid
                  "GGU":"B",
                  "GGC":"B",
                  "GGA":"B",
                  "GGG":"B", #glycine
                  
                  #stop codons
                  "UAA":"C",
                  "UAG":"E",
                  "UGA":"G"
                }

    def __init__(self, numTracks):
        '''initialize the instance of DNA_to_MIDI by setting the number of tracks in the MIDI file'''
        self.numTracks = numTracks
        self.mf = MIDIFile(self.numTracks) 
        
    def create_track(self, folder_name, track_name, track_num, nucleotides=None, tonic = 0, time = 0, duration = 1, repeat = False, 
                    codon = False, start_volume = 50):
        '''create_track creates the specified track of the midi file with the above parameters specified. Nucleotides refers to the string
        of DNA nucleotides for the desired organism and gene to be turned into music'''
        self.duration = duration
        self.repeat = repeat
        self.codon = codon
        self.volume = start_volume
        self.folder_name = folder_name

        if nucleotides == None:
            raise Exception("Be sure to pass in nucleotides!")
        self.nucleotidesDNA = nucleotides
        self.nucleotidesRNA = self.DNA_to_RNA(nucleotides) #convert to RNA for analysis
        self.base_pairsRNA = self.DNA_to_RNA(self.get_sequence_complement(nucleotides))
        self.compliment = self.get_sequence_complement(nucleotides)  

        self.tonic = tonic
        self.mediant = self.tonic + 4
        self.dominant = self.tonic + 7

        self.time = time
        self.volume = start_volume

        self.duration = duration
        self.repeat = repeat
        self.codon = codon

        self.mf.addTrackName(track_num, 0, track_name)
        self.mf.addTempo(track_num, self.time, 200)
        
        self.tonic_note = "C"
        self.mediant_note = "Eb"
        self.dominant_note = "G"

        DNA_to_MIDI.DNA_TO_CHROMATIC["C"] = self.tonic_note
        DNA_to_MIDI.DNA_TO_CHROMATIC["G"] = self.mediant_note
        DNA_to_MIDI.DNA_TO_CHROMATIC["A"] = self.mediant_note
        DNA_to_MIDI.DNA_TO_CHROMATIC["U"] = self.dominant_note	

        self.translation_keys_file = open(os.path.join(self.folder_name, track_name + '_translation_keys.txt'), 'w') #create file to write the keys to during translation

    def change_key(self, new_key, isMajor):
        '''changes the key given a new key and a mode (major or minor). The intervallic distance between the root of the current key and the root of the new key (e.g. C and A) is determined, 
        and, along with the mode, the tonic, mediant, and dominant notes of the new key (which are the 3 notes of the major/minor triad of that key) are determined'''
        semitones_up_scale = abs(DNA_to_MIDI.PITCH_DICTIONARY[self.tonic_note] - DNA_to_MIDI.PITCH_DICTIONARY[new_key])
        self.tonic += semitones_up_scale 
        if self.tonic >= len(DNA_to_MIDI.KEYS):
            self.tonic = self.tonic - len(DNA_to_MIDI.KEYS)

        if isMajor: 
            self.mediant = self.tonic + 4

        else: #minor 
            self.mediant = self.tonic + 3

        if self.mediant >= len(DNA_to_MIDI.KEYS):
            self.mediant = self.mediant - len(DNA_to_MIDI.KEYS) 

        self.dominant = self.tonic + 7
        if self.dominant >= len(DNA_to_MIDI.KEYS):
            self.dominant = self.dominant - len(DNA_to_MIDI.KEYS)

        self.tonic_note = DNA_to_MIDI.KEYS[self.tonic]
        self.mediant_note = DNA_to_MIDI.KEYS[self.mediant]
        self.dominant_note = DNA_to_MIDI.KEYS[self.dominant]

        DNA_to_MIDI.DNA_TO_CHROMATIC["C"] = self.tonic_note
        DNA_to_MIDI.DNA_TO_CHROMATIC["G"] = self.mediant_note
        DNA_to_MIDI.DNA_TO_CHROMATIC["A"] = self.mediant_note
        DNA_to_MIDI.DNA_TO_CHROMATIC["U"] = self.dominant_note

    def add_notes(self, track_num):
        '''adds notes to the MIDI file. Initial key is always C minor. Before the start codon AUG is encountered, individual nucleotides outline
        the 3 notes of minor triad (base pairs form dyads, or 2-note chords). Then, during translation, the key is major and the volume doubles. The
        key is determined by the codon (see the AMINO_ACIDS dictionary above). Individual nucleotides no longer determine the notes; this is now dictated
        by the codons. Each time the key changes based on the codon, the major triad of that key is played, until a stop codon is reached. Then we remain in
        the same minor key as the previous stop codon, but the volume is halved again and individual nucleotides once again outline the notes of the new
        minor triad. Then this process can repeat if another start codon is then encountered'''
        tripleStop = -1
        triple = False
        translation = False

        for nucleotideNum in range(len(self.nucleotidesRNA)):
            nucleotide = self.nucleotidesRNA[nucleotideNum]
            base_pair = self.base_pairsRNA[nucleotideNum]

            if tripleStop == nucleotideNum:
                triple = False

            if not triple and nucleotideNum < len(self.nucleotidesRNA) - 2:
                triple = True
                tripleStop = nucleotideNum + 3
                triplet_codon = nucleotide + self.nucleotidesRNA[nucleotideNum + 1] + self.nucleotidesRNA[nucleotideNum + 2]

                if not translation and triplet_codon in DNA_to_MIDI.START_CODONS:
                    self.volume = int(self.volume * 2)
                    self.change_key(DNA_to_MIDI.AMINO_ACIDS[triplet_codon], True)
                    self.duration = 2
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.tonic_note], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.mediant_note], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.dominant_note], self.time, self.duration, self.volume)
                    self.time += self.duration
                    
                    translation = True
                    
                    self.translation_keys_file.write(DNA_to_MIDI.AMINO_ACIDS[triplet_codon] + " ")
                
                elif not translation:
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[DNA_to_MIDI.DNA_TO_CHROMATIC[nucleotide]], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[DNA_to_MIDI.DNA_TO_CHROMATIC[base_pair]], self.time, self.duration, self.volume)
                    self.time += self.duration

                elif translation and triplet_codon in DNA_to_MIDI.STOP_CODONS:
                    self.volume = int(self.volume / 2)
                    self.change_key(DNA_to_MIDI.AMINO_ACIDS[triplet_codon], False) 

                    #longer (minor) chord to signify end of translation.
                    self.duration = 2
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.tonic_note], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.mediant_note], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.dominant_note], self.time, self.duration, self.volume)
                    self.time += self.duration

                    translation = False

                    self.translation_keys_file.write(DNA_to_MIDI.AMINO_ACIDS[triplet_codon] + " ")

                elif translation:
                    self.change_key(DNA_to_MIDI.AMINO_ACIDS[triplet_codon], True)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.tonic_note], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.mediant_note], self.time, self.duration, self.volume)
                    self.mf.addNote(track_num, DNA_to_MIDI.CHANNEL, DNA_to_MIDI.PITCH_DICTIONARY[self.dominant_note], self.time, self.duration, self.volume)
                    self.time += self.duration

                    self.translation_keys_file.write(DNA_to_MIDI.AMINO_ACIDS[triplet_codon] + " ")
            self.duration = 1
        self.translation_keys_file.close()

    def get_sequence_complement(self, sequence):
        '''get the genetic sequence complement (base pairs) with Biopython'''
        sequence = Seq(sequence)
        return str(sequence.complement())

    def DNA_to_RNA(self, dna_sequence):
        '''convert DNA to RNA with Biopython'''
        dna_sequence = Seq(dna_sequence)
        return str(dna_sequence.transcribe())

    def write_to_disk(self, fullfileName):
        '''write the MIDI file to the disk'''
        with open(fullfileName[:-4] + ".mid", 'wb') as outf:
            self.mf.writeFile(outf)

class Harmonic_Sequences:
    '''a class used for determining harmonic sequences of a given set of musical key, ostentiably derived from an organism
    in the creation of its MIDI sound file from its DNA. This class can also determine the intervallic differences between 
    a list of keys, as well as give a musicality rating to a particular organism based on its harmonic sequences'''
    semitones_from_C = {"C": 0,
                        "C#": 1,
                        "D": 2,
                        "Eb": 3,
                        "E": 4,
                        "F": 5,
                        "F#": 6,
                        "G": 7,
                        "Ab": 8,
                        "A": 9,
                        "Bb": 10,
                        "B": 11}

    def __init__(self, fullfilename):
        '''initialize the instance variable fullfilename to the user's choice. This is the name of the translation keys file for a desired organism
        that was generated during MIDI file generation for that organism's DNA'''
        self.fullfilename = fullfilename

    def get_harmonic_intervals(self):
        '''read in keys from the translation keys file. Each key is separated by a space, so we split the string of keys at the " " in order to get 
        a list of keys. Then, we go through the list and find the intervallic difference between each adjacent pairs of keys, based on the number of semitones each
        key is from C. This list of harmonic intervals is returned.'''
        with open(self.fullfilename, 'r') as f:
            DNA_keys = f.read().split(" ")[:-1] #remove the trailing empty string

        nucleotide_intervals = []
        key_num = 0
        for key_num in range(len(DNA_keys) - 1):
            start_key = Harmonic_Sequences.semitones_from_C[DNA_keys[key_num]]
            end_key = Harmonic_Sequences.semitones_from_C[DNA_keys[key_num + 1]]
            interval = end_key - start_key
            if interval < 0:
                interval = (11-start_key) + (end_key + 1)
            nucleotide_intervals.append(interval)

        return nucleotide_intervals

    def get_harmonic_sequences(self):
        '''Returns a list of length 4. Each element of this list corresponds to the percentage of the entire list of harmonic intervals corresponding to an organism's DNA that 
        belongs to each of the 4 main harmonic sequences: ascending 5ths (continually moving up a perfect fifth interval, or 7 semitones up), descending 5ths (continually moving 
        down a perfect fifth interval, or 7 semitones down, which, for our purposes, we are considering to be a perfect fourth up/5 semitones up since the range of our key changes 
        spans 1 octave so moving up a fifth and down a fourth become enharmonically equivalent/harmonically the same in our model), ascending 5-6 (down 5 semitones, up 7 semitones, 
        and repeat), and descending 5-6 (up 7 semitones, up 2 semitones, and repeat). For any of these sequences to be valid, they must be at least of length 3 (so have a minimum of 
        two intervals that fit the above criteria). In addition, sequences cannot overlap.'''

        desc_5_6_semitones_CMAJ_MODEL = ["C", "G", "A", "E", "F", "C", "D", "A", "B", "F", "G", "D", "E", "B"]
        asc_5_6_semitones_CMAJ_MODEL = ["C", "A", "D", "B", "E", "C", "F", "D", "G", "E", "A", "F", "B", "G"]
        #asc_5ths_semitones_CMAJ_MODEL = ["C", "G", "D", "A", "E", "B", "F#", "C#", "Ab", "Eb", "Bb", "F"]
        #desc_5ths_semitones_CMAJ_MODEL = ["C", "F", "Bb", "Eb", "Ab", "C#", "F#", "B", "E", "A", "D", "G"]

        desc_5_6_intervals = []
        for key_num in range(len(desc_5_6_semitones_CMAJ_MODEL) - 1):
            start_key = Harmonic_Sequences.semitones_from_C[desc_5_6_semitones_CMAJ_MODEL[key_num]]
            end_key = Harmonic_Sequences.semitones_from_C[desc_5_6_semitones_CMAJ_MODEL[key_num + 1]]
            interval = end_key - start_key
            if interval < 0:
                interval = (11-start_key) + (end_key + 1)
            desc_5_6_intervals.append(interval)

        asc_5_6_intervals = []
        for key_num in range(len(asc_5_6_semitones_CMAJ_MODEL) - 1):
            start_key = Harmonic_Sequences.semitones_from_C[asc_5_6_semitones_CMAJ_MODEL[key_num]]
            end_key = Harmonic_Sequences.semitones_from_C[asc_5_6_semitones_CMAJ_MODEL[key_num + 1]]
            interval = end_key - start_key
            if interval < 0:
                interval = (11-start_key) + (end_key + 1)
            asc_5_6_intervals.append(interval)

        nucleotide_intervals = self.get_harmonic_intervals()

        desc_5ths_interval = 5  #order is relative here (e.g. in our model we care about the key itself and not the octave/frequency of the note), so we do not care if the interval is up a fourth or down a fifth, these are equivalent for our purposes. we will interpret these as ascending 4ths. also perfect 5ths do not have modes
        asc_5ths_interval = 7 

        desc_5ths = False
        asc_5ths = False
        desc_5_6 = False
        asc_5_6 = False

        desc_5ths_counter = 0
        asc_5ths_counter = 0
        desc_5_6_counter = 0
        asc_5_6_counter = 0
        sequence_index = 0
        current_sequence_length = 0
        
        #next, analyze the text file for key changes
        for intervalNum in range(len(nucleotide_intervals)):
            interval = nucleotide_intervals[intervalNum]
            if sequence_index > 14: #(length - 1) of both asc 5 6 and desc 5 6 arrays
                sequence_index = 0
            
            #not in any sequence
            if not (desc_5ths or asc_5ths or desc_5_6 or asc_5_6):
                if intervalNum == len(nucleotide_intervals) - 1:
                    break #we don't consider sequences of length 1 (as it would be here, since this is the case when we are not currently in any sequence)
                if interval == desc_5ths_interval:
                    desc_5ths = True
                    current_sequence_length += 1
                elif interval == asc_5_6_intervals[sequence_index]:
                    asc_5_6 = True
                    current_sequence_length += 1
                    sequence_index += 1
                elif interval == asc_5ths_interval: # or interval == desc_5_6_intervals[sequence_index] (they start the same way)
                    asc_5ths = True
                    desc_5_6 = True
                    sequence_index += 1

            #the first interval is the same for desc 5-6 and asc 5ths		
            elif desc_5_6 and asc_5ths:
                #the next interval determines the sequence is desc 5-6
                if interval == desc_5_6_intervals[sequence_index]:
                    asc_5ths = False
                    current_sequence_length += 2
                    sequence_index += 1
                    if intervalNum == len(nucleotide_intervals) - 1: #we know, here, that the current sequence length is 2, which is a valid sequence
                        desc_5_6_counter += current_sequence_length
                #the next interval determines the sequence is asc 5ths
                elif interval == asc_5ths_interval: #asc 5ths
                    desc_5_6 = False
                    sequence_index = 0
                    current_sequence_length += 2
                    if intervalNum == len(nucleotide_intervals) - 1: #we know, here, that the current sequence length is 2, which is a valid sequence
                        asc_5ths_counter += current_sequence_length
                #the next interval does not correspond to either sequence, so this is not actually a sequence as sequences cannot have length 1
                else:
                    desc_5_6 = False
                    asc_5ths = False
                    sequence_index = 0
                    current_sequence_length = 0

            elif desc_5_6:
                if interval == desc_5_6_intervals[sequence_index]:
                    sequence_index += 1
                    current_sequence_length += 1
                    if intervalNum == len(nucleotide_intervals) - 1 and current_sequence_length > 1:
                        desc_5_6_counter += current_sequence_length 
                else:
                    if current_sequence_length > 1:
                        desc_5_6_counter += current_sequence_length
                    current_sequence_length = 0
                    desc_5_6 = False
                    sequence_index = 0
            elif asc_5ths:
                if interval == asc_5ths_interval:
                    current_sequence_length += 1
                    if intervalNum == len(nucleotide_intervals) - 1 and current_sequence_length > 1:
                        asc_5ths_counter += current_sequence_length 
                else:
                    if current_sequence_length > 1:
                        asc_5ths_counter += current_sequence_length
                    current_sequence_length = 0
                    asc_5ths = False
            elif desc_5ths:
                if interval == desc_5ths_interval:
                    current_sequence_length += 1
                    if intervalNum == len(nucleotide_intervals) - 1 and current_sequence_length > 1:
                        desc_5ths_counter += current_sequence_length 
                else:
                    if current_sequence_length > 1:
                        desc_5ths_counter += current_sequence_length
                    current_sequence_length = 0
                    desc_5ths = False
            else: #asc 5 6
                if interval ==  asc_5_6_intervals[sequence_index]:
                    sequence_index += 1
                    current_sequence_length += 1
                    if intervalNum == len(nucleotide_intervals) - 1 and current_sequence_length > 1:
                        asc_5_6_counter += current_sequence_length 
                else:
                    if current_sequence_length > 1:
                        asc_5_6_counter += current_sequence_length
                    current_sequence_length = 0
                    asc_5_6 = False
                    sequence_index = 0

        #determine the percentage of the total list of harmonic intervals for the given organism that corresponds to each of the 4 sequences.
        percent_desc_5_6 = desc_5_6_counter / len(nucleotide_intervals)
        percent_asc_5_6 = asc_5_6_counter / len(nucleotide_intervals)
        percent_desc_5ths = desc_5ths_counter / len(nucleotide_intervals)
        percent_asc_5ths = asc_5ths_counter / len(nucleotide_intervals)
        
        return [percent_desc_5_6, percent_asc_5_6, percent_desc_5ths, percent_asc_5ths]
    
    def get_musicality_rating(self):
        '''returns a "musicality rating" on a scale of 0-100, where a higher score means more musical. This is equal to the sum of the percentages
        of the list of harmonic intervals for the given organism that correspond to the 4 harmoinic sequences, now on a scale of 0-100. Thus, the greater
        amount of an organism's DNA that corresponds to a harmonic sequence in its musical representation, the more musical it is.'''
        total_percent_sequences = functools.reduce(lambda a,b : a+b, self.get_harmonic_sequences())
        return total_percent_sequences * 100

class Sequence_Classifier:
    def __init__(self, species_sequences):
        '''set the instance variable species_sequences to the user input. The input should be a list of length 4 consisting of the percentages of the 4 harmonic sequences contained in the DNA'''
        if not isinstance(species_sequences, list):
            raise Exception("Species harmonic sequence should be a list of length 4 consisting of the percentages of the 4 harmonic sequences contained in the DNA")
        self.species_sequences = species_sequences 
        
    def num_files_specify_filename_ending(self, path, file_criterion):
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

    def get_sequence_percentages(self, root_directory):
        '''returns a dictionary matching a class name (e.g. "Mammals") to a list of length 4 consisting of the average percentages of the 4 harmonic 
        sequences contained in the DNA of the species/organisms in that class'''
        classes_dictionary = {}
        AllFiles = list(os.walk(root_directory))

        for item in AllFiles:
            (foldername, _, LoFiles) = item  
            isRoot = foldername == root_directory
            invalidDataDir = isRoot or any(invalidStr in foldername for invalidStr in ['.git', '-']) # things like ICCC-author-kit and the .git folder
            total_sequences = [0, 0, 0, 0]

            if self.num_files_specify_filename_ending(foldername, "translation_keys.txt") == 0 and not invalidDataDir:
                raise Exception("A class folder must have some harmonic sequence text files for organisms representing species within that class. Populate the class folder " + foldername + " and try again")

            if not invalidDataDir:
                for filename in LoFiles:
                    if filename[-len("translation_keys.txt"):] == "translation_keys.txt":
                        fullfilename = os.path.join(foldername, filename)
                        organism_sequences = Harmonic_Sequences(fullfilename)
                        total_sequences = [current_total + new_value for current_total, new_value in zip(total_sequences, organism_sequences.get_harmonic_sequences())]
                average_class_sequences =  [x/self.num_files_specify_filename_ending(foldername, "translation_keys.txt") for x in total_sequences]
                classes_dictionary[foldername] = average_class_sequences
        
        return classes_dictionary

    def four_dim_distance(self, four_list1, four_list2):
        '''Find the "4-dimensional" distance between two lists, where each list is treated like a point in 4-dimensional space.
        Applies the distance formula in 4 dimensions: distance between (x1, y1, z1, w1) and (x2, y2, z2, w2) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 + (w2-w1)^2)'''
        if len(four_list1) != 4 or len(four_list2) != 4:
            raise Exception("Cannot compute 4-dimensional distance if length of each input vector is not 4!")
        
        sum_of_squared_distances = 0
        for i in range(4): #length of tuples should be 4
            sum_of_squared_distances += (four_list2[i] - four_list1[i]) ** 2
        return math.sqrt(sum_of_squared_distances)

    def classify_species(self):
        '''classifies a given species into a class (e.g. "Mammals") based on the harmonic sequences derived from its DNA.
        The percentages of this species/organism’s DNA that correspond to the four harmonic sequences are determined and packaged into a list. 
        Then, using the “four dimensional distance formula,” the “distance” between the new organism’s harmonic sequence and the averages for each 
        class is computed. The organism is then classified into the class whose average harmonic sequence percentages are “closest” to that of the organism’s.'''
        classes_sequences_dictionary = self.get_sequence_percentages(".")
        if len(classes_sequences_dictionary) == 0:
            raise Exception("Cannot classify a species to an empty class folder. Populate the class folder with harmonic sequence data and try again.")
        
        index = 0
        min_distance = 0
        min_class = ""
        for class_name, class_sequences in classes_sequences_dictionary.items():
            if index == 0: #initialize minimum distance to any element in the dictionary
                min_distance = self.four_dim_distance(class_sequences, self.species_sequences)
                min_class = class_name
                index += 1
            else:
                current_distance = self.four_dim_distance(class_sequences, self.species_sequences)
                if current_distance < min_distance:
                    min_distance = current_distance
                    min_class = class_name
        return min_class[2:] #remove the root directory text ".\"



'''At the bottom of the program are the 5 runner functions for this program that can be called as desired. If starting with an empty directory tree, 
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
    DNA_MIDI_converter = DNA_to_MIDI(num_tracks)
    with open(fullfilename, 'r') as file:
        nucleotides = file.read().replace('\n', '')
        DNA_MIDI_converter.create_track(folderpath, filename[:-4], track_num, nucleotides)
        DNA_MIDI_converter.add_notes(track_num)

    DNA_MIDI_converter.write_to_disk(fullfilename)

def write_nucleotides_and_generate_midi_txt_files(class_name, species_scientific_name, gene_name, class_folder_path = "."):
    '''runner function that takes in a class name, the scientific name for the desired species, the desired gene name, and optionally a class folder path 
    (if you don't want it to be created in the root directory. Creates an instance of the NCBI class, and then uses Biopython to access the NCBI (National Center for Biotechnology Information
    Then the user is asked to approve the DNA selection from the NCBI through terminal prompt, and once the data is approved, it is written to a .txt file with the species 
    scientific name to a folder with the class name. Then, the MIDI file and .txt file of musical key changes are generated for that DNA .txt file'''
    if(not os.path.isdir(os.path.join(class_folder_path, class_name))):
        os.mkdir(os.path.join(class_folder_path, class_name))
        
    ncbi = NCBI()
    ncbi.write_sequences_file_to_folder(class_name, class_folder_path, species_scientific_name, gene_name)
    generate_single_midi_and_txt(os.path.join(class_folder_path, class_name), species_scientific_name + ".txt", 1, 0)

def classify_species(organism_scientific_name):
    '''runner function that classifies a given species into a class (e.g. "Mammals"). Creates an instance of the Harmonic_Sequences class to get the harmonic sequences list (percent of DNA corresponding to 
    each of the 4 harmonic sequences). Then creates an instance of the Sequence_Classifier class, and uses the classify_species method to assign a class to the given organism repsenting a specific species.'''
    harmseq = Harmonic_Sequences(organism_scientific_name + "_translation_keys.txt") #"delphinapterus leucas.txt"
    harmseq_to_classify = harmseq.get_harmonic_sequences()
    seqclassifier = Sequence_Classifier(harmseq_to_classify)
    assigned_class = seqclassifier.classify_species()
    return assigned_class

def get_species_musicality_score(organism_scientific_name):
    '''runner function that gives a species a musicality score on a scale of 0-100 (higher score means more musical) through the Harmonic_Sequences class.'''
    harmseq = Harmonic_Sequences(organism_scientific_name + "_translation_keys.txt")
    musicality_score = harmseq.get_musicality_rating()
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

num_classified_correctly = 0
total_species = 0
for (class_name, species_list) in class_species_dict_TP53.items():
    for species in species_list:
        # write_nucleotides_and_generate_midi_txt_files(class_name, species, "TP53")
        print("SPECIES: ", species)
        classification = classify_species(class_name + "/" + species)
        if classification == class_name:
            num_classified_correctly += 1
        print("CLASSIFICATION: ", classification)
        print("CORRECT CLASSIFICATION: ", class_name)
        print("MUSICALITY SCORE: ", get_species_musicality_score(class_name + "/" + species))
        print()
        total_species += 1

print("PERCENT CLASSIFIED CORRECTLY: ", num_classified_correctly / total_species * 100)



