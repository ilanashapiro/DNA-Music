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

DNA_TO_CHROMATIC = {"A": "",  
                    "C": "",
                    "G": "",
                    "T": ""}

START_CODONS = ["ATG"]
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

TONIC_NOTE = "C"
MEDIANT_NOTE = "Eb"
DOMINANT_NOTE = "G"

def __init__(numTracks):
    '''initialize the instance of DNA_to_MIDI by setting the number of tracks in the MIDI file'''
    numTracks = numTracks
    mf = MIDIFile(numTracks) 
    
def create_track(folder_name, track_name, track_num, nucleotides=None, tonic = 0, time = 0, duration = 1, repeat = False, 
                codon = False, start_volume = 50):
    '''create_track creates the specified track of the midi file with the above parameters specified. Nucleotides refers to the string
    of DNA nucleotides for the desired organism and gene to be turned into music'''
    duration = duration
    repeat = repeat
    codon = codon
    volume = start_volume
    folder_name = folder_name

    if nucleotides == None:
        raise Exception("Be sure to pass in nucleotides!")
    nucleotidesDNA = nucleotides
    nucleotidesRNA = DNA_to_RNA(nucleotides) #convert to RNA for analysis
    base_pairsRNA = DNA_to_RNA(get_sequence_complement(nucleotides))

    tonic = tonic
    mediant = tonic + 4
    dominant = tonic + 7

    time = time
    volume = start_volume

    duration = duration
    repeat = repeat
    codon = codon

    mf.addTrackName(track_num, 0, track_name)
    mf.addTempo(track_num, time, 200)

    DNA_TO_CHROMATIC["C"] = TONIC_NOTE
    DNA_TO_CHROMATIC["G"] = MEDIANT_NOTE
    DNA_TO_CHROMATIC["A"] = MEDIANT_NOTE
    DNA_TO_CHROMATIC["U"] = DOMINANT_NOTE	

    translation_keys_file = open(os.path.join(folder_name, track_name + '_translation_keys.txt'), 'w') #create file to write the keys to during translation

def change_key(new_key, isMajor):
    '''changes the key given a new key and a mode (major or minor). The intervallic distance between the root of the current key and the root of the new key (e.g. C and A) is determined, 
    and, along with the mode, the tonic, mediant, and dominant notes of the new key (which are the 3 notes of the major/minor triad of that key) are determined'''
    semitones_up_scale = abs(PITCH_DICTIONARY[tonic_note] - PITCH_DICTIONARY[new_key])
    tonic += semitones_up_scale 
    if tonic >= len(KEYS):
        tonic = tonic - len(KEYS)

    if isMajor: 
        mediant = tonic + 4

    else: #minor 
        mediant = tonic + 3

    if mediant >= len(KEYS):
        mediant = mediant - len(KEYS) 

    dominant = tonic + 7
    if dominant >= len(KEYS):
        dominant = dominant - len(KEYS)

    tonic_note = KEYS[tonic]
    mediant_note = KEYS[mediant]
    dominant_note = KEYS[dominant]

    DNA_TO_CHROMATIC["C"] = tonic_note
    DNA_TO_CHROMATIC["G"] = mediant_note
    DNA_TO_CHROMATIC["A"] = mediant_note
    DNA_TO_CHROMATIC["U"] = dominant_note

def add_notes(track_num):
    '''adds notes to the MIDI file. Initial key is always C minor. Before the start codon AUG is encountered, individual nucleotides outline
    the 3 notes of minor triad (base pairs form dyads, or 2-note chords). Then, during translation, the key is major and the volume doubles. The
    key is determined by the codon (see the AMINO_ACIDS dictionary above). Individual nucleotides no longer determine the notes; this is now dictated
    by the codons. Each time the key changes based on the codon, the major triad of that key is played, until a stop codon is reached. Then we remain in
    the same minor key as the previous stop codon, but the volume is halved again and individual nucleotides once again outline the notes of the new
    minor triad. Then this process can repeat if another start codon is then encountered'''
    tripleStop = -1
    triple = False
    translation = False

    for nucleotideNum in range(len(nucleotidesRNA)):
        nucleotide = nucleotidesRNA[nucleotideNum]
        base_pair = base_pairsRNA[nucleotideNum]

        if tripleStop == nucleotideNum:
            triple = False

        if not triple and nucleotideNum < len(nucleotidesRNA) - 2:
            triple = True
            tripleStop = nucleotideNum + 3
            triplet_codon = nucleotide + nucleotidesRNA[nucleotideNum + 1] + nucleotidesRNA[nucleotideNum + 2]

            if not translation and triplet_codon in START_CODONS:
                volume = int(volume * 2)
                change_key(AMINO_ACIDS[triplet_codon], True)
                duration = 2
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dominant_note], time, duration, volume)
                time += duration
                
                translation = True
                
                translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")
            
            elif not translation:
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[DNA_TO_CHROMATIC[nucleotide]], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[DNA_TO_CHROMATIC[base_pair]], time, duration, volume)
                time += duration

            elif translation and triplet_codon in STOP_CODONS:
                volume = int(volume / 2)
                change_key(AMINO_ACIDS[triplet_codon], False) 

                #longer (minor) chord to signify end of translation.
                duration = 2
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dominant_note], time, duration, volume)
                time += duration

                translation = False

                translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")

            elif translation:
                change_key(AMINO_ACIDS[triplet_codon], True)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[tonic_note], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[mediant_note], time, duration, volume)
                mf.addNote(track_num, CHANNEL, PITCH_DICTIONARY[dominant_note], time, duration, volume)
                time += duration

                translation_keys_file.write(AMINO_ACIDS[triplet_codon] + " ")
        duration = 1
    translation_keys_file.close()

def get_sequence_complement(sequence):
    '''get the genetic sequence complement (base pairs) with Biopython'''
    sequence = Seq(sequence)
    return str(sequence.complement())

def DNA_to_RNA(dna_sequence):
    '''convert DNA to RNA with Biopython'''
    dna_sequence = Seq(dna_sequence)
    return str(dna_sequence.transcribe())

def write_to_disk(fullfileName):
    '''write the MIDI file to the disk'''
    with open(fullfileName[:-4] + ".mid", 'wb') as outf:
        mf.writeFile(outf)

