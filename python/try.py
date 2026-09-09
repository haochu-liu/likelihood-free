import numpy as np
import pandas as pd
from Bio import Phylo
from tqdm import tqdm
import sys
from pathlib import Path
current_path = Path(__file__).resolve()
project_root = current_path.parent.parent
pysimARG_path = project_root / 'pysimARG'
data_path = project_root / 'data'

sys.path.append(str(pysimARG_path))
from segment_summary_stats import segment_summary_stats
from clonal_genealogy import ClonalTree



clonal_tree = Phylo.read(str(data_path / "klebsiella" / "klebsiella_clonal.nwk"), "newick")
Phylo.draw_ascii(clonal_tree)

clonal_edge = np.loadtxt(str(data_path / "klebsiella" / "clonal_edge.csv"), delimiter=",", dtype=float)
clonal_node_height = np.loadtxt(str(data_path / "klebsiella" / "clonal_node_height.csv"), delimiter=",", dtype=float)
print(clonal_edge.shape, clonal_node_height.shape)

bool_mat = np.loadtxt(str(data_path / "klebsiella" / "genomes_bool.csv"), delimiter=",", dtype=bool)
print(bool_mat.shape, bool_mat.dtype)

np.random.seed(100)
clonal_tree = ClonalTree(n=100)
clonal_tree.edge = clonal_edge
clonal_tree.node_height = clonal_node_height
clonal_tree.height = np.max(clonal_node_height)
clonal_tree.length = np.sum(clonal_edge[:, 2])

rand_seg_df = pd.read_csv(str(data_path / "klebsiella" / "rand_seg_info.csv"))
print(rand_seg_df.head())

rand_seg1000_summary = np.full((2000, 46), np.nan)

for i in tqdm(range(2000), desc="Processing genome segments"):
    start_pos = rand_seg_df.loc[i, 'Start_pos']
    end_pos = rand_seg_df.loc[i, 'End_pos']
    
    seg_matrix = bool_mat[:, start_pos-1:end_pos]  # Adjust for 0-based indexing
    summary_stats = segment_summary_stats(clonal_tree, seg_matrix)
    rand_seg1000_summary[i, :] = summary_stats

    if (i + 1) % 100 == 0:
        np.savetxt(str(data_path / "klebsiella" / "rand_seg1000.csv"), rand_seg1000_summary, delimiter=",")
