# Codon dictionary

codons = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

stop_codons = {'TAA', 'TAG', 'TGA'}

# Store sequence names and sequences in a dictionary from FASTA file

def parse_FASTA(FASTA_file):
    sequences = {} # Empty dictionary to store sequences (names are keys, sequences are values)
    name, seq = '', '' # Initialize empty name and sequence
    for line in FASTA_file:
        if line.startswith('>'):
            if name:
                sequences[name] = seq # Save the previous sequence
            name = line[1:] # Remove '>'
            seq = ''  # Reset seq to empty string, ready to build the next sequence
        else: # If the line does not start with > → it is part of the sequence
            seq += line.strip()  # Append sequence lines, removing whitespace
    if name:
        sequences[name] = seq # Save the last sequence 
    return sequences
    

# Translate from a start codon

def translate(seq, start):  # start → where to begin translating (the starting position, like 0, 1, or 2).
    amino_acid_seq = '' 
    for i in range(start, len(seq) - 2, 3):  # Move through the sequence three letters at a time.     #  Slicing beyond the end doesn't crash, but you would get an incomplete codon So we stop earlier at len(seq) - 2, to make sure seq[i:i+3] will always have 3 full letters.
        codon = seq[i:i+3]
        if codon in stop_codons:
            break
        if codon in codons:
            amino_acid_seq += codons[codon]
        else:
            amino_acid_seq += '-'  # '-' for unknown codon
    return amino_acid_seq



# Argparse for command line

import argparse  # helps your program accept input from the command line.

parser = argparse.ArgumentParser() # Create a parser object to handle command-line arguments.
parser.add_argument('input') #  We defined
parser.add_argument('output') # We defined
args = parser.parse_args()


# Read input FASTA file

with open(args.input, 'r') as f:
    fasta_seqs = parse_FASTA(f)

# Process and translate

translated = {} # empty dictionary (keys are headers of proteins and values are aa sequences)
for name in fasta_seqs:
    seq = fasta_seqs[name]
    for i in range(len(seq) - 2):
        if seq[i:i+3] == 'ATG':
            amino_acid = translate(seq, i) # starts at position i and translates codons into aa
            protein_name = name + '-' + str(i)
            translated[protein_name] = amino_acid
            

# Print output

with open(args.output, 'w' ) as out:
    for header in translated:
        protein = translated[header]
        print(">" + header, file = out) # >header
        print(protein, file = out) # Write the actual protein sequence
