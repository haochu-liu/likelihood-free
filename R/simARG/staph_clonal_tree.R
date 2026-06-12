library(ape)
library(phangorn)


# Load tree from ClonalFrameML and convert to clonal tree
t=read.tree('data/staph/CFML_tree/cfml_results.labelled_tree.newick')
d=cophenetic.phylo(t)
t2=upgma(d)
e=dist.nodes(t2)[Ntip(t2),1]/2#current root height
e2=2*(1-1/Ntip(t2))#coalescent expectation of root height
t2$edge.length=t2$edge.length/e*e2

# Plot both trees
plot(t, show.tip.label = F)
plot(t2, show.tip.label = F)

# Save t2
write.tree(phy=t2, file="data/staph/saureus_clonal.nwk")
