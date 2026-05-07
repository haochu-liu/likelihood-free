This is a small package that translate and modify R package [simARG](https://github.com/haochu-liu/simARG) into python.

`birth_death_sim.py`: simulate birth death process as an example.

`tree.py`: define a parental class for tree structures.

`clonal_genealogy`: a subclass `ClonalTree` and simulation method for clonal genealogy tree.

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

`homoplasy_index_simbac.py`: compute the homoplasy index from SimBac simulated data.

`ClonalOrigin_pair_sim.py`: simulate summary statistics by ClonalOrigin pair models with mutations.

`ClonalOrigin_seq_sim.py`: simulate summary statistics by ClonalOrigin seq models with mutation.

`discrete_uniform.py`: a function to simulate from discrete uniform in torch settings.

`fasta_to_bool.py`: convert `.fasta` sequences to a boolean matrix.

`newick_to_tree.py`: convert a Newick tree to a ClonalTree object.
