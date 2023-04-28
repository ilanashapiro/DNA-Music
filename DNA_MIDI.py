from Bio.Seq import Seq
from midiutil.MidiFile import MIDIFile
import os, os.path
import enum

'''a module to convert DNA (in the form of a text file of nucleotides) into a MIDI file'''
CHANNEL = 0 #set instrument to a piano

class ChordQuality(enum.Enum):
   MAJOR = 1
   MINOR = 2
   DIMINISHED = 3

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

STOP_CODONS = ["TAA", "TGA", "TAG"]
AMINO_ACIDS = { #essential
                "ATG":"C", #methionine/start codon
                "ATT":"C#",
                "ATC":"C#",
                "ATA":"C#", #isoleucine
                "TTA":"C#",
                "TTG":"C#",
                "CTT":"C#",
                "CTC":"C#",
                "CTA":"C#",
                "CTG":"C#", #leucine -- group with isoleucine since they're most similar by Grantham's distance
                "AAA":"D",
                "AAG":"D", #lysine
                "ACT":"Eb",
                "ACC":"Eb",
                "ACA":"Eb",
                "ACG":"Eb", #threonine
                "TTT":"E",
                "TTC":"E", #phenylalanine
                "TGG":"F", #tryptophan
                "CAT":"F#",
                "CAC":"F#", #histidine
                "GTT":"G",
                "GTC":"G",
                "GTA":"G",
                "GTG":"G", #valine

                #nonessential
                "AAT":"Ab",
                "AAC":"Ab", #asparagine
                "GAT":"A",
                "GAC":"A", #aspartate 
                "GCT":"Bb",
                "GCC":"Bb",
                "GCA":"Bb",
                "GCG":"Bb", #alanine
                "GAA":"B",
                "GAG":"B", #glutamate
                
                #conditionally essential
                "TAT":"C",
                "TAC":"C", #tyrosine
                "TGT":"D",
                "TGC":"D", #cysteine
                "TCC":"E",
                "TCT":"E",
                "TCA":"E",
                "TCG":"E", 
                "AGT":"E",
                "AGC":"E", #serine
                "AGA":"F",
                "AGG":"F", 
                "CGT":"F",
                "CGC":"F",
                "CGA":"F",
                "CGG":"F", #arginine
                "CCT":"G",
                "CCC":"G",
                "CCA":"G",
                "CCG":"G", #proline
                "CAA":"A",
                "CAG":"A", #glutamine
                "GGT":"B",
                "GGC":"B",
                "GGA":"B",
                "GGG":"B", #glycine
                
                #stop codons
                "TAA":"C",
                "TAG":"E",
                "TGA":"G" }

COMPLEMENT = {
    "A":"T",
    "T":"A",
    "C":"G",
    "G":"C",
}

CONDITIONALLY_ESSENTIAL = ["TAT",
                "TAC", #tyrosine
                "TGT",
                "TGC", #cysteine
                "TCC",
                "TCT",
                "TCA",
                "TCG", 
                "AGT",
                "AGC", #serine
                "AGA",
                "AGG", 
                "CGT",
                "CGC",
                "CGA",
                "CGG", #arginine
                "CCT",
                "CCC",
                "CCA",
                "CCG", #proline
                "CAA",
                "CAG", #glutamine
                "GGT",
                "GGC",
                "GGA",
                "GGG"] #glycine

def change_key(curr_tonic, curr_tonic_note_name, new_key, chord_quality, dna_to_chromatic_dict):
    '''changes the key given a new key and a mode (major or minor). The intervallic distance between the root of the current key and the root of the new key (e.g. C and A) is determined, 
    and, along with the mode, the tonic, mediant, and dominant notes of the new key (which are the 3 notes of the major/minor triad of that key) are determined'''
    
    semitones_diff = PITCH_DICTIONARY[curr_tonic_note_name] - PITCH_DICTIONARY[new_key]
    semitones_up_scale = abs(semitones_diff) if semitones_diff < 0 else len(KEYS) - semitones_diff

    new_tonic = (curr_tonic + semitones_up_scale) % len(KEYS)
    mediant = (new_tonic + (4 if chord_quality == ChordQuality.MAJOR else 3)) % len(KEYS) 
    dominant = (new_tonic + (6 if chord_quality == ChordQuality.DIMINISHED else 7)) % len(KEYS) 
    
    tonic_note_name = KEYS[new_tonic]
    mediant_note_name = KEYS[mediant]
    dominant_note_name = KEYS[dominant]

    dna_to_chromatic_dict["A"] = tonic_note_name
    dna_to_chromatic_dict["C"] = mediant_note_name
    dna_to_chromatic_dict["T"] = mediant_note_name
    dna_to_chromatic_dict["G"] = dominant_note_name

    return (new_tonic, tonic_note_name, mediant_note_name, dominant_note_name)

def get_submediant_note_name(curr_tonic):
    submediant = curr_tonic + 9
    if submediant >= len(KEYS):
        submediant = submediant - len(KEYS)
    return KEYS[submediant]

def process_UTR_frag(UTR_frag, midi_file, track_num, dna_to_chromatic_dict, time, duration, volume):
    for nucleotide in UTR_frag:
        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dna_to_chromatic_dict[nucleotide.upper()]], time, duration, volume)
        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dna_to_chromatic_dict[COMPLEMENT[nucleotide].upper()]], time, duration, volume)
        time += duration
    return time
        
def add_notes(folder_name, track_name, full_file_name, UTR_5prime_exons_list, CDS_list):
    '''adds notes to the MIDI file. Initial key is always C minor. Before the start codon AUG is encountered, individual nucleotides outline
    the 3 notes of minor triad (base pairs form dyads, or 2-note chords). Then, during translation, the key is major and the volume doubles. The
    key is determined by the codon (see the AMINO_ACIDS dictionary above). Individual nucleotides no longer determine the notes; this is now dictated
    by the codons. Each time the key changes based on the codon, the major triad of that key is played, until a stop codon is reached. Then we remain in
    the same minor key as the previous stop codon, but the volume is halved again and individual nucleotides once again outline the notes of the new
    minor triad. Then this process can repeat if another start codon is then encountered'''

    '''create_track creates the specified track of the midi file with the above parameters specified. Nucleotides refers to the string
    of DNA nucleotides for the desired organism and gene to be turned into music'''

    tonic_note_name = "C"
    mediant_note_name = "Eb"
    dominant_note_name = "G"

    dna_to_chromatic_dict = {}
    dna_to_chromatic_dict["A"] = tonic_note_name
    dna_to_chromatic_dict["C"] = mediant_note_name
    dna_to_chromatic_dict["T"] = mediant_note_name
    dna_to_chromatic_dict["G"] = dominant_note_name	

    volume = 50
    duration = 1
    time = 0
    tonic = 0
    track_num = 0

    midi_file = MIDIFile(1) 
    midi_file.addTrackName(track_num, 0, track_name)
    midi_file.addTempo(track_num, 0, 200)

    translation_keys_file = open(os.path.join(folder_name, track_name + '_translation_keys.txt'), 'w') #create file to write the keys to during translation

    UTR_frag_index = 0
    translation_initiated = False
    translation_completed = False
    start_codon_encountered = False
    remaining_length_in_codon = 0 # for codons across splice regions
    is_spliced = False
    regionIdx = 0
    
    while regionIdx < len(CDS_list):
        region = CDS_list[regionIdx]

        if UTR_frag_index < len(UTR_5prime_exons_list):
            change_key(tonic, tonic_note_name, tonic_note_name, ChordQuality.MINOR, dna_to_chromatic_dict)
            time = process_UTR_frag(UTR_5prime_exons_list[UTR_frag_index], midi_file, track_num, dna_to_chromatic_dict, time, 1, volume)
            UTR_frag_index += 1
            if UTR_frag_index == len(UTR_5prime_exons_list): # translation begins when we finish the UTR
                translation_initiated = True

        # then the entire region is part of the 5' UTR based on how the list is preprocessed
        # so the current exon is empty since we just processed the entire exon (UTR) 
        # and so we want to skip it and immediately process the next intron now which is at the next index
        if len(region) == 0: 
            regionIdx += 1
            region = CDS_list[regionIdx]

        nucleotideIdx = 0
        is_exon = regionIdx % 2 == 0
        
        while nucleotideIdx < len(region):
            nucleotide = region[nucleotideIdx]
            base_pair = COMPLEMENT[nucleotide]

            if is_exon:
                # we are in the middle of translation
                if translation_initiated and not translation_completed:
                    if nucleotideIdx == 0 and remaining_length_in_codon > 0: # if we're in a new exon and the codon began in the previous region (end of splice)
                        prev_exon = CDS_list[regionIdx - 2] # bc the next region in the list is an intron
                        initial_codon_frag_length = 3 - remaining_length_in_codon
                        initial_codon_frag = prev_exon[-initial_codon_frag_length:]
                        remaining_codon_frag = region[:remaining_length_in_codon]
                        triplet_codon = initial_codon_frag + remaining_codon_frag
                        nucleotideIdx += remaining_length_in_codon
                        remaining_length_in_codon = 0
                        # is_spliced remains true from the else clause
                    elif nucleotideIdx < len(region) - 2: # normal codon
                        triplet_codon = region[nucleotideIdx:nucleotideIdx+3]
                        nucleotideIdx += 3
                        is_spliced = False
                    else: # the codon begins across a splice region
                        initial_codon_frag = region[nucleotideIdx:len(region)]
                        remaining_length_in_codon = 3 - len(initial_codon_frag)
                        next_exon = CDS_list[regionIdx + 2] # bc the next region in the list is an intron
                        remaining_codon_frag = next_exon[:remaining_length_in_codon]
                        triplet_codon = initial_codon_frag + remaining_codon_frag
                        nucleotideIdx = len(region)
                        is_spliced = True

                    if not start_codon_encountered: # i.e. we're just beginning translation, the next codon should be ATG or the file is corrupt
                        if triplet_codon != "ATG":
                            raise Exception("Invalid sequence, 5' UTR is not followed by start codon, it is followed by", triplet_codon)
                        start_codon_encountered = True
                        volume = 100
                        duration = 0.5 if is_spliced else 2 
                        (tonic, tonic_note_name, mediant_note_name, dominant_note_name) = change_key(tonic, tonic_note_name, AMINO_ACIDS[triplet_codon], ChordQuality.MAJOR, dna_to_chromatic_dict)
                        dominant_pitch = PITCH_DICTIONARY[dominant_note_name] + (1 if triplet_codon in CONDITIONALLY_ESSENTIAL else 0) # we want augmented chord in this case
                        
                        #longer (major) chord to signify start of translation.
                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note_name], time, duration, volume)
                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration, volume)
                        midi_file.addNote(track_num, CHANNEL, dominant_pitch, time, duration, volume)
                        
                        translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")
                        time += duration

                    elif triplet_codon in STOP_CODONS:
                        volume = 50
                        duration = 0.5 if is_spliced else 2 
                        
                        (tonic, tonic_note_name, mediant_note_name, dominant_note_name) = change_key(tonic, tonic_note_name, AMINO_ACIDS[triplet_codon], ChordQuality.MINOR, dna_to_chromatic_dict) 
                        
                        #longer (minor) chord to signify end of translation.
                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note_name], time, duration, volume)
                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration, volume)
                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dominant_note_name], time, duration, volume)

                        translation_completed = True
                        translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")
                        time += duration
                        
                        # process the remainder of the exon as 3' UTR
                        duration = 1
                        time = process_UTR_frag(region[nucleotideIdx:], midi_file, track_num, dna_to_chromatic_dict, time, duration, volume)
                        nucleotideIdx = len(region)
                    else: # we are actively translating
                        volume = 100
                        duration = 0.5 if is_spliced else 1

                        (tonic, tonic_note_name, mediant_note_name, dominant_note_name) = change_key(tonic, tonic_note_name, AMINO_ACIDS[triplet_codon], ChordQuality.MAJOR, dna_to_chromatic_dict)
                        dominant_pitch = PITCH_DICTIONARY[dominant_note_name] + (1 if triplet_codon in CONDITIONALLY_ESSENTIAL else 0) # we want augmented chord in this case

                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note_name], time, duration, volume)
                        midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration, volume)
                        midi_file.addNote(track_num, CHANNEL, dominant_pitch, time, duration, volume)
                        
                        translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")
                        time += duration

                # we are dealing with entire exons in the 3' UTR following translation
                else: 
                    time = process_UTR_frag(region, midi_file, track_num, dna_to_chromatic_dict, time, duration, volume)
                    nucleotideIdx += 1 # since individual nucleotides, rather than codons, dictate the harmony here

            else: # we are in an intron
                duration = 1
                volume = 50
                (tonic, _, _, _) = change_key(tonic, tonic_note_name, tonic_note_name, ChordQuality.DIMINISHED, dna_to_chromatic_dict)
                midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dna_to_chromatic_dict[nucleotide.upper()]], time, duration, volume)
                midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dna_to_chromatic_dict[base_pair.upper()]], time, duration, volume)
                midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[get_submediant_note_name(tonic)], time, duration, volume)
                nucleotideIdx += 1
                time += duration

            duration = 1
        regionIdx += 1
    write_to_disk(full_file_name, midi_file)

def get_sequence_complement(sequence):
    '''get the genetic sequence complement (base pairs) with Biopython'''
    sequence = Seq(sequence)
    return str(sequence.complement())

def get_uppercase_prefix(input_str):
    c = input_str[0]
    idx = 0
    res = ""
    while idx < len(input_str) and c.isupper(): 
        res += c
        c = input_str[idx]
        idx += 1
    return res
        

def write_to_disk(full_file_name, midi_file):
    '''write the MIDI file to the disk'''
    with open(full_file_name[:-4] + ".mid", 'wb') as outf:
        midi_file.writeFile(outf)

