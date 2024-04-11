# Sample free associations for each cue, then produce cue x cue matrix (cosine similarity) and get back networks of cues to compute transitivity, degree, aspl

# Simulation to recreate environment and cognition and entropy/similarity measures

set.seed(1)
##### Rescorla Wagner function ######

rescorlaWagner <- function(vmat = vmat, cue, outcome, alpha=1, beta=.2) {
  # cue can be compound: e.g., cue = c("A", "B")
  # outcome can be compound: e.g., outcome = c("A", "B")
  cuevalues <- matrix(vmat[cue, ], nrow = length(cue))
  valueSums <- apply(cuevalues,2,sum )     # cumulative cue value
  # lambda = 1 when present 0 when absent
  lambda = as.numeric(rownames(vmat) %in% outcome)
  pError <- lambda - valueSums     # prediction error = observed - expected
  valueChange <- alpha * beta * pError # value change 
  vmat[cue, ] <- t(t(vmat[cue, ]) + valueChange)                   # update value
  return(vmat)
}


# World parameters (500/2000/200 for)
wordsInWorld=500
Associates = 2000 # Edges in environment
# Set learning events and age classes
learningEvents <- 200 
ageEpochs = 4
Worlds =50 
# Here we fix beta at .1 and use r^-1 (alpha = 1) for rank-based network
Betas = .1
alphas <- 1
alphai = 1
bis = 1
# Data storage for network metrics by age
# trans
transdev <- matrix(NA, nrow=4, ncol=Worlds)
# degree
degreedev <- matrix(NA, nrow=4, ncol=Worlds)
# aspl 
distancedev <- matrix(NA, nrow=4, ncol=Worlds)

# Each World creates a new environment and a new developmental learning trajectory with associations
for(woi in 1:Worlds){
  # Build the world
  x <- 1:wordsInWorld
  a = alphas[alphai] # set to 1 for association network demonstration
  pairs <- c()
  for(i in 1:Associates){
    pairs <- rbind(pairs, sample(x, size = 2, prob=x^(-a), replace = FALSE))
  }
  # Make graph from pairs
  ii <- graph_from_edgelist(pairs,directed=FALSE) 
  
  # Add isolates if number of vertices is not 500
  difamt = 500-length(V(ii))
  if (difamt > 0){
    ii <-  add_vertices(ii,difamt)
  }
  
  # Rename graph for learning representation 
  iis <- ii
  
  #### Prepare Representation Matrix ####
  
  n = length(V(iis))+1 # number of cues (words) + context cue
  
  # Initialize zero value matrix for learning
  vmat <- matrix(0, nrow = n, ncol = n)
  rownames(vmat) <- 1:n
  colnames(vmat) <- 1:n
  
  # Initialize metrics
  edgeE <- rep(NA, ageEpochs)
  nodeCount <- rep(NA, ageEpochs)
  ageNetworkList <- list(NULL)
  
  ##### Learn Representations across epochs #####
  
  for(lage in 1:ageEpochs){
    # take fraction of environment that is learnable 
    #ageWords <- round(wordsInWorld*y[lage]/y[length(y)])
    ageWords <- wordsInWorld
    #iis <- subgraph(ii, 1:ageWords) 
    # make training data set
    traindata <- c()
    # sample edges from environment
    for(i in 1:learningEvents){
      traindata <- rbind(traindata, cue_outcome <- ends(iis, sample(E(iis), 1, prob=E(iis)$weight))) 
    }
    # keep list of first words learned in year 1
    if(lage == 1){
      firsttraindata <- traindata
    }
    # train on cue-outcome pairs
    for(i in 1:learningEvents){
      cue_outcome <- traindata[i,]
      cue_outcome<- sample(as.vector(cue_outcome))
      vmat <- rescorlaWagner(vmat, cue=c(n, cue_outcome[1]), outcome=cue_outcome[2], beta = Betas[bis]) 
    }
    # Make undirected graph from representation 
    gle <- graph_from_adjacency_matrix(vmat, weighted=TRUE, diag=FALSE, mode = "undirected")
    # Remove context cue
    gle <- igraph::delete_vertices(gle, n)
    # take subgraph of learnable words (only nec. if growing)
    gle <- igraph::subgraph(gle, 1:ageWords)
    nodeCount[lage] <- length(V(gle))
    # Remove negative edges
    negEdges <- which(E(gle)$weight < 0)
    gle <- igraph::delete_edges(gle, negEdges)
    # save learned representation
    ageNetworkList[[lage]] <- gle
  } # End development representation construction
  
  #### Compute associations for each network
  # Number of cues
  numcues = 30 
  # Make cue list from nodes with degree 3 or more across all four epochs
  deg3list.1 <- which(degree(ageNetworkList[[1]]) > 3)
  deg3list.2 <- which(degree(ageNetworkList[[2]]) > 3)
  deg3list.3 <- which(degree(ageNetworkList[[3]]) > 3)
  deg3list.4 <- which(degree(ageNetworkList[[4]]) > 3)
  # Take intersection of all lists
  gclist <- intersect(intersect(intersect(deg3list.1,deg3list.2), deg3list.3),deg3list.4)
  # Make cue list from sample of intersection of all lists
  cuelist <- sample(gclist, numcues)
  # store networks and matrices
  gcs <- list(NULL)
  gmats <- list(NULL)
  # For each age network 
  for(nits in 1:4){
    # Get weighted adjacency matrix for producing associates at each age
    getassmat <- get.adjacency(ageNetworkList[[nits]], attr="weight", sparse = FALSE)
    # Make empty cue x target matrix
    cxtm <- matrix(0, nrow= length(cuelist), ncol = length(V(ageNetworkList[[4]])))
    ## For number of participants 
    numparticipants = 10
    ## Generate up to 3 associates each 
    for(cuei in 1:length(cuelist)){
      # For number of participants
      for(partis in 1:numparticipants){
        # Sample from row of adjacency matrix in proportion to weight
        threeasss <- sample(1:ncol(cxtm), 3, prob=getassmat[cuelist[cuei],])
        # Add to cue x target matrix
        cxtm[cuei,threeasss] <- cxtm[cuei,threeasss] + 1
      } 
    }
    #### Compute cue x cue network
    cuesims <-  lsa::cosine(t(cxtm)) 
    # set diagonal to 0
    diag(cuesims) <- 0
    #### Compute stats on network
    gcs[[nits]] <- cuesims
    #### 
    gmats[[nits]] <- graph_from_adjacency_matrix(cuesims, weighted = TRUE, mode="undirected")
  } # end production of each cue x cue network across development
  
  ## compute median sim across all worlds as network threshold  
  medw <- median(unlist(gcs))
  # for each network threshold then compute values
  for(ip in 1:length(gcs)){
    nettodo <- gmats[[ip]]
    nettodo <- delete_edges(nettodo, which(E(nettodo)$weight < medw))
    E(nettodo)$weight <- 1
    # transitivity
    transdev[ip, woi] <- mean(igraph::transitivity(nettodo, type="localundirected"), na.rm=TRUE)
    # degree
    degreedev[ip, woi] <- mean(igraph::degree(nettodo))
    # aspl
    distancedev[ip, woi] <- igraph::mean_distance(nettodo)
  }
}
  
pdf(file="NetworkAssociations.pdf", width = 7, height = 4) 
par(mfrow=c(1,3)) 
par(mar=c(5,5,2,2))
trm <- rowMeans(transdev)
dgrm <- rowMeans(degreedev)
dirm <- rowMeans(distancedev)
sdtrm <- apply(transdev, 1, sd)
sddgrm <- apply(degreedev, 1, sd)
sddirm <- apply(distancedev, 1, sd)
plot(dgrm, cex.lab=1.5, xaxt="n", ylab = "Degree", xlab = "Epoch", ylim=c(12, 18), xlim = c(.5, 4.5))
arrows(1:4, dgrm+sddgrm, 1:4, dgrm-sddgrm, angle=90, code=3, length = .1)
axis(1, at=1:4, labels=c("1","2","3","4") )
plot(dirm, cex.lab=1.5, xaxt="n", ylab="Average shortest path length", xlab="Epoch", ylim = c(1.3, 1.6), xlim = c(.5, 4.5))
axis(1, at=1:4, labels=c("1","2","3","4") )
arrows(1:4, dirm+sddirm, 1:4, dirm-sddirm, angle=90, code=3, length = .1)
plot(trm, cex.lab=1.5, xaxt="n", ylab = "Clustering Coefficient", xlab = "Epoch", ylim = c(0.5,.8), xlim=c(.5, 4.5))
arrows(1:4, trm+sdtrm, 1:4, trm-sdtrm, angle=90, code=3, length = .1)
axis(1, at=1:4, labels=c("1","2","3","4") )
dev.off()
