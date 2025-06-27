## UCL_Decoding_Project

This project reads one or more nucleotide sequences from a FASTA file and translates them into amino acid sequences based on the standard genetic code.

For each sequence:

- All ATG start codons (methionine) are identified.
- From each start codon, translation continues codon-by-codon until a stop codon (TAA, TAG, or TGA) is found.
- If no stop codon is encountered before the end of the sequence, translation proceeds to the end, ignoring any incomplete codon (1–2 bases)
