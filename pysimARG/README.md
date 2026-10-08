This is a small package that translate and modify R package [simARG](https://github.com/haochu-liu/simARG) into python.

`birth_death_sim.py`: simulate birth death process as an example.

`tree.py`: define a parental class for tree structures.

`clonal_genealogy.py`: a subclass `ClonalTree` and simulation method for clonal genealogy tree.

`ClonalOrigin_ARG.py`: a subclass `ARG` and ClonalOrigin simulation with pair or seq model.

`ClonalOrigin_nodes.py`: find the row index for recombination nodes in `ARG` simulator.

`pair_simulator.py`: ClonalOrigin pair model called by `ARG`.

`seq_simulator.py`: ClonalOrigin seq model called by `ARG`.

`add_mutation.py`: simulate mutations for the ARG object.

`add_mutation_truncated.py`: simulate mutations for the ARG object given the sites all polymorphic.

`localtree.py`: pick the local tree from the ARG object.

`localtree_simbac.py`: functions to load and select local trees.

`G3_test.py`: compute three-gamete test.

`G4_test.py`: compute four-gamete test.

`LD_r.py`: compute the square of correlation coefficient for LD.

`homoplasy_index.py`: compute the homoplasy index for a given ARG and leaf node data.

`ClonalOrigin_pair_sim.py`: simulate summary statistics by ClonalOrigin pair models with mutations.

`ClonalOrigin_seq_sim.py`: simulate summary statistics by ClonalOrigin seq models with mutation.

`discrete_uniform.py`: a function to simulate from discrete uniform in torch settings.

`fasta_to_bool.py`: convert `.fasta` sequences to a boolean matrix.

`newick_to_tree.py`: convert a Newick tree to a ClonalTree object.

`Watterson_theta.py`: compute Watterson's theta estimator.

`Tajima_pi.py`: compute Tajima's pi and Wakeley's pi^2.

`Tajima_D.py`: compute Tajima's D as a normalized difference between Tajima's pi and Watterson's theta.

`LD.py`: compute D, D', r^2 for linkage disequilibrium.

`Kelly_Z.py`: compute Kelly's Z_nS for the given sequence (repeat the same loop in `seq_sim`).

`Hudson_Rm.py`: compute Hudson's R_M estimator as the minimal number of recombinations

`Wall_BQ.py`: compute Wall's B and Q statistics.

`exp_regression.py`: fit an exponential regression model to the given data and provide coefficients.

`segment_summary_stats.py`: provide summary statistics given a gene matrix.

`extract_blast_segment.py`: extract sequence segments from a FASTA file based on BLAST output.

`evaluate_posterior_metrics.py`: compute evaluation metrics for posterior samples across multiple runs.

`LeaveLengthOut_NN.py`: provide embedding network for SBI while keeping length dimension out.

`hotspot_DoG.py`: a tool to detect hotspot candidates using DoG.

`hotspot_HMM.py`: a tool to detect hotspot candidates using HMM.

## Simulation functions

`ClonalTree`
-> `clonal_genealogy.py`, `tree.py`



`ClonalOrigin_seq_sim`
-> `ClonalOrigin_ARG.py`, `add_mutation.py`, `segment_summary_stats.py`,
   `seq_simulator.py`, `ClonalOrigin_nodes.py`,
   `G4_test.py`, `LD.py`, `homoplasy_index.py`, `Watterson_theta.py`,
   `Tajima_pi.py`, `Tajima_D.py`, `Wall_BQ.py`, `exp_regression.py`
