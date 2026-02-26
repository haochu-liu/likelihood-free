edge <- matrix(NA, nrow = 9, ncol = 3)
edge[, 1] <- c(6, 6, 7, 7, 8, 8, 9, 10, 10)
edge[, 2] <- c(4, 5, 3, 6, 1, 2, 7, 8, 9)

node_vec <- 1:10
node_site_vec <- c(0, 0, 1, 0, 1, NA, NA, NA, NA, NA)
n_leaf <- 5
# Compute actual changes
s_site <- 0
site_list <- setNames(rep(list(NA), length(node_vec)), node_vec)
site_list[1:n_leaf] <- as.integer(node_site_vec[1:n_leaf])
for (i in (n_leaf+1):length(node_vec)) {
  parent_node <- node_vec[i]
  node_index <- which(edge[, 1] == parent_node)
  if (length(node_index) == 2) {
    # Coalescent structure
    children_node <- edge[node_index, 2]
    child_1 <- site_list[[as.character(children_node[1])]]
    child_2 <- site_list[[as.character(children_node[2])]]
    intersec <- intersect(child_1, child_2)
    if (length(intersec) == 0) {
      # intersection is empty -> mutation
      site_list[[as.character(parent_node)]] <- union(child_1, child_2)
      s_site <- s_site + 1
    } else {
      # intersection is not empty -> no mutation
      site_list[[as.character(parent_node)]] <- intersec
    }
  } else if (length(node_index) == 1) {
    # Recombination structure
    children_node <- edge[node_index, 2]
    site_list[[as.character(parent_node)]] <- site_list[[as.character(children_node)]]
  }
}
s_site
