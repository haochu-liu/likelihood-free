import numpy as np
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
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


clonal_edge = np.loadtxt(str(data_path / "staph" / "clonal_edge.csv"), delimiter=",", dtype=float)
clonal_node_height = np.loadtxt(str(data_path / "staph" / "clonal_node_height.csv"), delimiter=",", dtype=float)
genomes_bool = np.loadtxt(str(data_path / "staph" / "genomes_bool.csv"), delimiter=",", dtype=bool)

core_gene_info = pd.read_csv(str(data_path / "staph" / "core_gene_info.csv"))
core_lengths_array = core_gene_info['Gene_Length'].to_numpy()
end_pos_array = core_gene_info['End_pos'].to_numpy()

num_genes = len(core_lengths_array)
genomes_length = genomes_bool.shape[1]

seg_length = [200, 500, 1000, 2000, 5000]
back_seg_df = pd.DataFrame({'Gene_Length': np.hstack([np.repeat(seg_len, num_genes) for seg_len in seg_length])})
back_seg_df['End_pos'] = np.hstack([end_pos_array for seg_len in seg_length])
back_seg_df['Start_pos'] = back_seg_df['End_pos'] - back_seg_df['Gene_Length'] + 1
back_seg_df['Alignment'] = True
back_seg_df

np.random.seed(100)
clonal_tree = ClonalTree(n=110)

clonal_tree.edge = clonal_edge
clonal_tree.node_height = clonal_node_height
clonal_tree.height = np.max(clonal_node_height)
clonal_tree.length = np.sum(clonal_edge[:, 2])

def process_segment(i, start_pos, end_pos):
    if start_pos < 1:
        return i, None 

    seg_matrix = genomes_bool[:, start_pos-1:end_pos] 
    summary_stats = segment_summary_stats(clonal_tree, seg_matrix)
    return i, summary_stats


if __name__ == "__main__":
    num_segments = back_seg_df.shape[0]
    back_seg_summary_stats = np.full((num_segments, 46), np.nan)


    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(process_segment, i, row['Start_pos'], row['End_pos']): i 
            for i, row in back_seg_df.iterrows()
        }
        for future in tqdm(as_completed(futures), total=num_segments, desc="Processing genome segments"):
            i, summary_stats = future.result()
            if summary_stats is not None:
                back_seg_summary_stats[i, :] = summary_stats

    back_seg_df.to_csv(str(data_path / "staph" / "back_seg_df.csv"), index=False)
    np.savetxt(str(data_path / "staph" / "back_seg_summary_stats.csv"), back_seg_summary_stats, delimiter=",")


