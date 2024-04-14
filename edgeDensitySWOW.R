swow <- read.delim("SampleDataFilesBNS/strength.SWOW-EN.R1.csv", sep = "\t")

# Edge list of forward association strengths
swow <- tibble(swow)
# Require at least x people to produce associations, this can be changed
swowT <- swow %>% filter(R1 >=2)
# Make graph from two columns of edge list
sg <- graph_from_edgelist(as.matrix(swowT[,1:2]), directed=TRUE)

# remove all nodes with 0 outdegree, so we only consider the ~12000 cues

delnod <- which(degree(sg, mode = "out")<=0)

sgnoout <- delete.vertices(sg, delnod)

sgnoout <- as.undirected(
  sgnoout,
  mode = "collapse")
  
igraph::edge_density(sgnoout) # .0015
