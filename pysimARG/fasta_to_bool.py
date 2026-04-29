import numpy as np
from Bio import SeqIO


def fasta_to_bool(fasta_file):
    '''Provide a fasta file path as the input.'''
    sequences = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_chars = list(str(record.seq).upper())
        sequences.append(seq_chars)

    char_matrix = np.array(sequences)

    reference_seq = char_matrix[0]
    bool_matrix = (char_matrix != reference_seq)
    
    return bool_matrix
