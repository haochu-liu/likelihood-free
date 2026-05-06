import re
import io
from Bio import Phylo

def load_local_trees(local_tree_file):
    genome_blocks = []

    # ^\[(\d+)\]  --> Matches [number] at the start of the line and captures the number
    # (.*)        --> Captures everything else (the Newick string)
    pattern = re.compile(r'^\[(\d+)\](.*)')
    
    print("Parsing local trees...")
    
    with open(local_tree_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
                
            match = pattern.match(line)
            if match:
                span_length = int(match.group(1))

                newick_str = match.group(2)

                tree_obj = Phylo.read(io.StringIO(newick_str), "newick")

                genome_blocks.append((span_length, tree_obj))
                
    return genome_blocks


def get_tree_at_position(genome_blocks, target_site):
    """
    Finds the local tree that covers a specific site on the genome.
    
    Parameters:
    - genome_blocks: List of tuples (span_length, tree_object)
    - target_site: Integer representing the genomic position (1-based indexing)
    
    Returns:
    - tree_object, block_start, block_end
    """
    if target_site < 1:
        print("Error: Target site must be 1 or greater.")
        return None, None, None

    current_start = 1
    
    for span, tree in genome_blocks:

        current_end = current_start + span - 1

        if current_start <= target_site <= current_end:
            return tree, current_start, current_end

        current_start += span

    print(f"Error: Site {target_site} is beyond the simulated genome length (Max: {current_start - 1}).")
    return None, None, None
