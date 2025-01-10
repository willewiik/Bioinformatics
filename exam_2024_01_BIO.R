

# Q1

# a) #

# (S1, S1, S2), 1/54
# (S1, S1, S1), 


# b) #
# (S1, S1, S1), 2/81
# (S1, S1, S2), 2/162
# (S1, S2, S2), 1/6
# Total prob = 2/162 + 1/6 + 2/81 = 31/162



# c) #




# Q2
load("//filur01.it.liu.se/students/wilwi856/Downloads/gene_expression_measurements.RData")


loess_fit <- loess(gene_expression_measurements[,1] ~ gene_expression_measurements[,2])
plot(gene_expression_measurements[,1],gene_expression_measurements[,2])
lines(gene_expression_measurements[,1], predict(loess_fit), col = "blue", lwd = 2)

M <- gene_expression_measurements[,1] - gene_expression_measurements[,2]

# Upregulated genes: M > 2 (arbitrary threshold for example)
upregulated <- which(M > 2)

# Downregulated genes: M < -2 (arbitrary threshold for example)
downregulated <- which(M < -2)


# Q3

# Lecture 7

# 1.
adj_matrix <- matrix(0, nrow = 8, ncol = 8)

# Define edges based on the graph structure
edges <- list(
  c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 2),
  c(5, 7), c(5, 8), c(7, 8), c(6, 7)
)

# Populate the adjacency matrix
for (edge in edges) {
  adj_matrix[edge[1], edge[2]] <- 1
  adj_matrix[edge[2], edge[1]] <- 1 # Symmetric for an undirected graph
}

# Print the adjacency matrix
print(adj_matrix)

# 2
# A graph is bipartite if its nodes can be divided into two disjoint sets such 
# that no two nodes within the same set are adjacent. Alternatively, a graph is 
# bipartite if it does not contain any odd-length cycles.
# therefore, not biparte


# 3
# 1 -> 2 -> 3 -> 4 -> 5 -> 1


# 4
# The diameter of a graph is the longest shortest path between any two nodes.
# 1 to 8, lentgh of 4


# 5
library(igraph)
#bla bla, betweeness()





# EXAM


# Q2

seq1 <- "CGTTC"
seq2 <- "CGATC"

alignment1 <- align_sequences(seq1, seq2, d = -1, mismatch = -1, match = 1, method="needleman")
plot(alignment1)
alignment1
