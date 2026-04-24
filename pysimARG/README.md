This is a small package that translate and modify R package [simARG](https://github.com/haochu-liu/simARG) into python.

`tree.py`: define a parental class for tree structures.

`clonal_genealogy`: a subclass `ClonalTree` and simulation method for clonal genealogy tree.

`ClonalOrigin_ARG.py`: a subclass `ARG` and ClonalOrigin simulation method for pair of sites.

`ClonalOrigin_nodes.py`: find the row index for recombination nodes in `ARG` simulator.

`pair_simulator.py`: ClonalOrigin pair model called by `ARG`.

`seq_simulator.py`: ClonalOrigin seq model called by `ARG`.

`add_mutation.py`: simulate mutations for the ARG object.

`add_mutation_truncated.py`: simulate mutations for the ARG object given the sites all polymorphic.

`localtree.py`: pick the local tree from the ARG object.

`G3_test.py`: compute three-gamete test.

`LD_r.py`: compute the square of correlation coefficient for LD.

`homoplasy_index.py`: compute the homoplasy index for a given ARG and leaf node data.

`ClonalOrigin_pair_sim.py`: simulate summary statistics by ClonalOrigin pair models with mutations.

`ClonalOrigin_seq_sim.py`: simulate summary statistics by ClonalOrigin seq models with mutation.

`discrete_uniform.py`: a function to simulate from discrete uniform in torch settings.
