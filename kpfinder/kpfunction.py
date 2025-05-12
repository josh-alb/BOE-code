import re
import collections

def kp(file_path, tonic):
    """
    Analyze a music file and return the degrees array.
    
    Args:
        file_path: Path to the music file
        tonic: Tonic note (a-g with optional accidentals)
    
    Returns:
        degrees: Array of 12 values representing scale degree counts
    """
    # token pattern
    pattern = re.compile(r'([a-g])\1*([#\-n]*)')
    
    # count each note
    notes = collections.Counter()
    with open(file_path, 'r') as file:
        for line in file:
            # remove comments and metadata
            if line and line[0] in '!*':
                continue

            # handle barlines
            if line and line[0] == '=':
                if len(line) > 1 and line[1] == '9':
                    break
                else:
                    continue

            # update note counts
            notes.update(re.findall(pattern, line.lower()))

    # note degrees in C major
    note_degree = {'c': 0, 'd': 2, 'e': 4, 'f': 5, 'g': 7, 'a': 9, 'b': 11}
    offsets = {'#': 1, '-': -1, '': 0, 'n': 0}

    # transpose to tonic
    tonic_note, tonic_accidental = re.match(r'([a-g])([#-]*)', tonic).groups()
    note_degree = {
        note: (degree - note_degree[tonic_note] - sum(offsets[a] for a in tonic_accidental)) % 12
        for note, degree in note_degree.items()
    }

    # map note to tonic chromatic scale degrees
    degrees = [0 for _ in range(12)]
    for (name, accidental), count in notes.items():
        degree = (note_degree[name] + sum(offsets[a] for a in accidental)) % 12
        degrees[degree] += count

    return degrees
