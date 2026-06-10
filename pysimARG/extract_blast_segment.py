import pandas as pd
import numpy as np
from Bio import SeqIO


def extract_blast_segments(blast_df, fasta_path, query_id, subject_ids):
    """
    Extracts sequence segments from a FASTA file based on BLAST output.
    Return a np array of extracted sequences from subject genomes.
    
    Assumes blast_df has standard column names:
    ['Query_ID', 'Subject_ID', 'Identity', 'Alignment_Length', 'Mismatches', 'Gap_Openings',
     'Query_Start', 'Query_End', 'Subject_Start', 'Subject_End', 'E-value', 'Bit_Score']
    """
    
    df_filtered = blast_df[blast_df['Query_ID'] == query_id].copy()

    # Handle Duplicates: Keep only the best hit for each subject genome
    df_sorted = df_filtered.sort_values(by='Bit_Score', ascending=False)
    df_best = df_sorted.drop_duplicates(subset='Subject_ID', keep='first')
    
    # Load the genomes into memory as a dictionary
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    
    extracted_segments = []
    for subj_id in subject_ids:
        match_row = df_best[df_best['Subject_ID'] == subj_id]
        
        # Handle missing matches
        if match_row.empty:
            print(f"Missing Result: No match found for genome '{subj_id}' with query '{query_id}'.")
            continue

        if subj_id not in genome_dict:
            print(f"Missing Genome: '{subj_id}' was in BLAST results, but not found in the FASTA file.")
            continue
            
        # Get coordinates
        sstart = int(match_row['Subject_Start'].iloc[0])
        send = int(match_row['Subject_End'].iloc[0])
        
        # Get the full genome sequence
        full_seq = genome_dict[subj_id].seq
        
        # Extract the segment
        if sstart < send:
            segment = full_seq[sstart - 1 : send]
        else:
            raw_segment = full_seq[send - 1 : sstart]
            segment = raw_segment.reverse_complement()

        # Append the extracted sequence string to our results
        extracted_segments.append(list(str(segment).upper()))
        
    return np.array(extracted_segments)
