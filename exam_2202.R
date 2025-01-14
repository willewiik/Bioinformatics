



tree <- ape::read.nexus("~/Downloads/phyltree_20220207.nex")
plot(tree)

# Q1
ape::plot.phylo(tree, cex = 0.5, adj = 0.1)

# Q2
ape::degree(tree)


# Q3
ape::dist.nodes(tree)[99,100]

# Q4
ape::dist.nodes(tree)[3,68]

# Q5
ape::dist.nodes(tree)[71,80]

ape::nodepath(tree, 3, 68)
