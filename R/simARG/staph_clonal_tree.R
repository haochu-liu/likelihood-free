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
plot(t2, show.tip.label = T)

n_leaf = Ntip(t2)
edge = t2$edge
leaf_height = rep(0, n_leaf)
names(leaf_height) = c(1:n_leaf)
node_height = c(leaf_height, branching.times(t2))
node_height = sort(node_height)

# Change node index
# Node 1 - 110: unchanged
# Node 111 - 219: ordered by height
node_order <- as.integer(names(node_height))
edge[] <- match(edge, node_order)

edge_full = matrix(NA, nrow=nrow(edge), ncol=3)
edge_full[, 1:2] = edge
edge_full[, 3] = t2$edge.length
edge_full <- edge_full[order(edge_full[, 1]), ]

write.table(edge_full, file="data/staph/clonal_edge_new.csv",
            sep = ",", row.names=FALSE, col.names=FALSE)
write.table(node_height, file="data/staph/clonal_node_height_new.csv",
            sep = ",", row.names=FALSE, col.names=FALSE)
write.table(t2$tip.label, file="data/staph/tip_names.csv",
            sep = ",", row.names=FALSE, col.names=FALSE)

# Save t2
write.tree(phy=t2, file="data/staph/saureus_clonal.nwk")
