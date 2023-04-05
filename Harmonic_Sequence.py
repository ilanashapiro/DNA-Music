import functools

'''a module used for determining harmonic sequences of a given set of musical key, ostentiably derived from an organism
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

def get_harmonic_intervals(fullfilename):
    '''read in keys from the translation keys file. Each key is separated by a space, so we split the string of keys at the " " in order to get 
    a list of keys. Then, we go through the list and find the intervallic difference between each adjacent pairs of keys, based on the number of semitones each
    key is from C. This list of harmonic intervals is returned.'''
    with open(fullfilename, 'r') as f:
        DNA_keys = f.read().split(" ")[:-1] #remove the trailing empty string

    nucleotide_intervals = []
    key_num = 0
    for key_num in range(len(DNA_keys) - 1):
        start_key = semitones_from_C[DNA_keys[key_num]]
        end_key = semitones_from_C[DNA_keys[key_num + 1]]
        interval = end_key - start_key
        if interval < 0:
            interval = (11-start_key) + (end_key + 1)
        nucleotide_intervals.append(interval)

    return nucleotide_intervals

def get_harmonic_sequences(fullfilename):
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
        start_key = semitones_from_C[desc_5_6_semitones_CMAJ_MODEL[key_num]]
        end_key = semitones_from_C[desc_5_6_semitones_CMAJ_MODEL[key_num + 1]]
        interval = end_key - start_key
        if interval < 0:
            interval = (11-start_key) + (end_key + 1)
        desc_5_6_intervals.append(interval)

    asc_5_6_intervals = []
    for key_num in range(len(asc_5_6_semitones_CMAJ_MODEL) - 1):
        start_key = semitones_from_C[asc_5_6_semitones_CMAJ_MODEL[key_num]]
        end_key = semitones_from_C[asc_5_6_semitones_CMAJ_MODEL[key_num + 1]]
        interval = end_key - start_key
        if interval < 0:
            interval = (11-start_key) + (end_key + 1)
        asc_5_6_intervals.append(interval)

    nucleotide_intervals = get_harmonic_intervals(fullfilename)

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

def get_musicality_rating(fullfilename):
    '''returns a "musicality rating" on a scale of 0-100, where a higher score means more musical. This is equal to the sum of the percentages
    of the list of harmonic intervals for the given organism that correspond to the 4 harmoinic sequences, now on a scale of 0-100. Thus, the greater
    amount of an organism's DNA that corresponds to a harmonic sequence in its musical representation, the more musical it is.'''
    total_percent_sequences = functools.reduce(lambda a,b : a+b, get_harmonic_sequences(fullfilename))
    return total_percent_sequences * 100