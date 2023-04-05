from Bio.Seq import Seq
from midiutil.MidiFile import MIDIFile
import os, os.path

'''a module to convert DNA (in the form of a text file of nucleotides) into a MIDI file'''
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

STOP_CODONS = ["TAA", "TGA", "TAG"]
AMINO_ACIDS = { #essential
                "ATG":"C", #methionine/start codon
                "ATT":"C#",
                "ATC":"C#",
                "ATA":"C#", #isoleucine
                "AAA":"D",
                "AAG":"D", #lysine
                "ACT":"Eb",
                "ACC":"Eb",
                "ACA":"Eb",
                "ACG":"Eb", #threonine
                "TTT":"E",
                "TTC":"E", #phenylalanine
                "TGG":"F", #tryptophan
                "TTA":"F#",
                "TTG":"F#",
                "CTT":"F#",
                "CTC":"F#",
                "CTA":"F#",
                "CTG":"F#", #leucine
                "CAT":"G",
                "CAC":"G", #histidine
                "GTT":"Ab",
                "GTC":"Ab",
                "GTA":"Ab",
                "GTG":"Ab", #valine

                #nonessential
                "AAT":"A",
                "AAC":"A", #asparagine
                "GAT":"Bb",
                "GAC":"Bb", #aspartate 
                "GCT":"B",
                "GCC":"B",
                "GCA":"B",
                "GCG":"B", #alanine
                
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
                "CAG":"A",
                "GAA":"A",
                "GAG":"A", #glutamine/glutamic acid
                "GGT":"B",
                "GGC":"B",
                "GGA":"B",
                "GGG":"B", #glycine
                
                #stop codons
                "TAA":"C",
                "TAG":"E",
                "TGA":"G" }

def change_key(curr_tonic, new_key, isMajor, dna_to_chromatic_dict):
    '''changes the key given a new key and a mode (major or minor). The intervallic distance between the root of the current key and the root of the new key (e.g. C and A) is determined, 
    and, along with the mode, the tonic, mediant, and dominant notes of the new key (which are the 3 notes of the major/minor triad of that key) are determined'''
    semitones_up_scale = abs(PITCH_DICTIONARY[curr_tonic] - PITCH_DICTIONARY[new_key])
    new_tonic = curr_tonic + semitones_up_scale 
    if new_tonic >= len(KEYS):
        new_tonic = new_tonic - len(KEYS)

    if isMajor: 
        mediant = new_tonic + 4

    else: #minor 
        mediant = new_tonic + 3

    if mediant >= len(KEYS):
        mediant = mediant - len(KEYS) 

    dominant = new_tonic + 7
    if dominant >= len(KEYS):
        dominant = dominant - len(KEYS)

    dna_to_chromatic_dict["A"] = KEYS[new_tonic]
    dna_to_chromatic_dict["C"] = KEYS[mediant]
    dna_to_chromatic_dict["T"] = KEYS[mediant]
    dna_to_chromatic_dict["G"] = KEYS[dominant]

    return new_tonic

def add_notes(folder_name, track_name, full_file_name, nucleotides=None):
    '''adds notes to the MIDI file. Initial key is always C minor. Before the start codon AUG is encountered, individual nucleotides outline
    the 3 notes of minor triad (base pairs form dyads, or 2-note chords). Then, during translation, the key is major and the volume doubles. The
    key is determined by the codon (see the AMINO_ACIDS dictionary above). Individual nucleotides no longer determine the notes; this is now dictated
    by the codons. Each time the key changes based on the codon, the major triad of that key is played, until a stop codon is reached. Then we remain in
    the same minor key as the previous stop codon, but the volume is halved again and individual nucleotides once again outline the notes of the new
    minor triad. Then this process can repeat if another start codon is then encountered'''

    '''create_track creates the specified track of the midi file with the above parameters specified. Nucleotides refers to the string
    of DNA nucleotides for the desired organism and gene to be turned into music'''

    if nucleotides == None:
        raise Exception("Be sure to pass in nucleotides!")

    tonic_note_name = "C"
    mediant_note_name = "Eb"
    dominant_note_name = "G"

    dna_to_chromatic_dict = {}
    dna_to_chromatic_dict["A"] = tonic_note_name
    dna_to_chromatic_dict["C"] = mediant_note_name
    dna_to_chromatic_dict["T"] = mediant_note_name
    dna_to_chromatic_dict["G"] = dominant_note_name	

    nucleotides_complement = get_sequence_complement(nucleotides)
    volume = 50
    duration = 1
    time = 0
    tonic = 0
    track_num = 0

    midi_file = MIDIFile(1) 
    midi_file.addTrackName(track_num, 0, track_name)
    midi_file.addTempo(track_num, 0, 200)

    translation_keys_file = open(os.path.join(folder_name, track_name + '_translation_keys.txt'), 'w') #create file to write the keys to during translation

    nucleotideIdx = 0
    translation_initiated = False
    translation_completed = False
    codon_splice_fragment = ""

    while nucleotideIdx < len(nucleotides):
        nucleotide = nucleotides[nucleotideIdx]
        base_pair = nucleotides_complement[nucleotideIdx]

        # if nucleotideIdx < len(nucleotides) - 2:
        if nucleotide.isupper():
            # we are in the middle of translation
            if translation_initiated and not translation_completed:
                triplet_codon = nucleotides[nucleotideIdx:nucleotideIdx+2] # should never go out of bounds since the data should always contain a stop codon
                
                if nucleotide.isupper() and triplet_codon in STOP_CODONS:
                    volume = int(volume / 2)
                    tonic = change_key(tonic, AMINO_ACIDS[triplet_codon], False, dna_to_chromatic_dict) 

                    #longer (minor) chord to signify end of translation.
                    midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note_name], time, duration*2, volume)
                    midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration*2, volume)
                    midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration*2, volume)

                    translation_completed = True
                    translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")

                else:
                    tonic = change_key(tonic, AMINO_ACIDS[triplet_codon], True, dna_to_chromatic_dict)
                    midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note_name], time, duration, volume)
                    midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration, volume)
                    midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration, volume)

                    translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")

                nucleotideIdx += 3

            # we are in an exon, we have not begun translating, and we encounter a start codon for the first time. translation begins
            elif not translation_initiated and nucleotideIdx < len(nucleotides) - 2 and nucleotides[nucleotideIdx:nucleotideIdx+2] == "ATG": 
                triplet_codon = nucleotides[nucleotideIdx:nucleotideIdx+2]
                volume = int(volume * 2)
                tonic = change_key(tonic, AMINO_ACIDS[triplet_codon], True, dna_to_chromatic_dict)

                #longer (major) chord to signify start of translation.
                midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note_name], time, duration*2, volume)
                midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration*2, volume)
                midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note_name], time, duration*2, volume)
                
                translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")
                translation_initiated = True
                nucleotideIdx += 3
                
            # we are in either the 5' or 3' UTR of an exon
            else:
                nucleotideIdx += 1

        else: # lowercase, we are in an intron
            midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dna_to_chromatic_dict[nucleotide]], time, duration, volume)
            midi_file.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dna_to_chromatic_dict[base_pair]], time, duration, volume)
            nucleotideIdx += 1
                
        time += duration

    translation_keys_file.close()

    write_to_disk(full_file_name, midi_file)

def get_sequence_complement(sequence):
    '''get the genetic sequence complement (base pairs) with Biopython'''
    sequence = Seq(sequence)
    return str(sequence.complement())

def write_to_disk(full_file_name, midi_file):
    '''write the MIDI file to the disk'''
    with open(full_file_name[:-4] + ".mid", 'wb') as outf:
        midi_file.writeFile(outf)

