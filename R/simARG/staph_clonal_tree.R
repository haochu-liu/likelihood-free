library(ape)
library(phangorn)
t=read.tree('data/staph/CFML_tree/cfml_results.labelled_tree.newick')
d=cophenetic.phylo(t)
t2=upgma(d)
e=dist.nodes(t2)[Ntip(t2),1]/2#current root height
e2=2*(1-1/Ntip(t2))#coalescent expectation of root height
t2$edge.length=t2$edge.length/e*e2

plot(t, show.tip.label = F)
plot(t2, show.tip.label = F)
