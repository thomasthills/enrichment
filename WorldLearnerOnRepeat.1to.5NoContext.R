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
Worlds = 10
Betas = seq(.1,.5, .1)

alphas <- c(0,1)

Elist <- list(NULL)
Slist <- list(NULL)

for(alphai in 1:length(alphas)){
 
  EEB <- matrix(NA, nrow=length(Betas), ncol = 4)
  SSB <- matrix(NA, nrow=length(Betas), ncol = 4)
  
  for(bis in 1:length(Betas)){
    # Entropy keeper
    EE <- matrix(NA, nrow=Worlds, ncol = ageEpochs)
    # Similarity keeper
    SS <- matrix(NA, nrow=Worlds, ncol = ageEpochs)
    
    for(woi in 1:Worlds){
        
      # Build the world
      x <- 1:wordsInWorld
      a = alphas[alphai] # set to 1 for ranking and 0 for ER with fixed number
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
      
      ##### Learn Representation #####
      
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
          #vmat <- rescorlaWagner(vmat, cue=c(n, cue_outcome[1]), outcome=cue_outcome[2], beta = Betas[bis]) 
          # Do above without context
          vmat <- rescorlaWagner(vmat, cue=c(cue_outcome[1]), outcome=cue_outcome[2], beta = Betas[bis]) 
        }
        # make undirected graph from representation 
        gle <- graph_from_adjacency_matrix(vmat, weighted=TRUE, diag=FALSE, mode = "undirected")
        # remove context cue
        gle <- igraph::delete_vertices(gle, n)
        # take subgraph of learnable words (only nec. if growing)
        gle <- igraph::subgraph(gle, 1:ageWords)
        nodeCount[lage] <- length(V(gle))
        # Remove negative edges
        negativeEdges <- which(E(gle)$weight < 0)
        gle <- igraph::delete_edges(gle, negativeEdges)
        # save learned representation
        ageNetworkList[[lage]] <- gle
        # copy network
        g2 <- gle
        # get symmetric weighted network 
        weightMatrix <- as_adjacency_matrix(g2, attr="weight", sparse=FALSE)
        # set negative values to zero
        weightMatrix[weightMatrix < 0] <- 0
        # normalize rows for entropy
        ww<-weightMatrix/rowSums(weightMatrix, na.rm=TRUE)
        # entropy function
        entropy <- function(p) rowSums(-(p * log(p)), na.rm = TRUE)
        # compute entropy for each node
        node_entropy <- entropy(ww) 
        
        # only compute entropy for words learned about at time 1 (removes isolates with 0 entropy)
        ftlist <- unique(c(firsttraindata[,1], firsttraindata[,2]))
        # take median entropy of list
        edgeE[lage] <- mean(node_entropy[ftlist])
        # remove edges with 0 or less weight
        g2 <- igraph::delete.edges(gle, which(E(gle)$weight <=0.0))
        # set remaining weights to 1
        E(g2)$weight <- 1
        # plot learned representation
        # plot(g2, vertex.size = 1, edge.arrow.size = 0, vertex.label=NA, layout=layout_with_fr(g2))
        # label first one in 'learned'
        # if(lage == 1){
         #  text(0, 1.5, "Learned lexicon")
        # }
        # label all with iterations
        # its <- lage*learningEvents
        # text(0, -1.5, paste("t =",its))
      }
      # Entropy values to save 
      EE[woi,] <- edgeE
      
      #### Compute and Plot Similarity
      
      # choose target pairs, from first training data (firsttraindata)
      learningPairs = 20 # this 50 random pairs
      simpair <- firsttraindata[sample(seq(1:nrow(firsttraindata)), learningPairs ),]
      # simpair <- firsttraindata # this takes all pairs
      # set initial activation levels
      simpair <- data.frame(simpair, activation = 100)
      
      # time for spreading activation
      tspr <- 10 
      # assign column names
      names(simpair) <- c("node", "node", "activation")
      ## create datastore for similarity judgments
      simJudge <- c()
      ## for each age network
      for(sage in 1:length(ageNetworkList)){
        ## initiate activation at each node and measure activation at the other 
        for(testrow in 1:nrow(simpair)){
          # initiate activation at the cue
          df1 <- spreadr::spreadr(start_run = simpair[testrow,c(1,3)], decay = 0,
                                  retention = 0, suppress = 0,
                                  network = ageNetworkList[[sage]], time = tspr)
          # measure max activation at the target node
          maxActivation12 <- max(subset(df1, node == simpair[testrow,2])$activation)
          # initiate activation at the other cue
          df2 <- spreadr::spreadr(start_run = simpair[testrow,c(2,3)], decay = 0,
                                  retention = 0, suppress = 0,
                                  network = ageNetworkList[[sage]], time = tspr)
          # measure max activation at the target node
          maxActivation21 <- max(subset(df2, node == simpair[testrow,1])$activation)
          # add the max activations
          simval <- maxActivation12 + maxActivation21
          # add results to the data frame
          simJudge <- rbind(simJudge, c(simpair[testrow,1],simpair[testrow,2], simval, sage))
        }
      }
      # ready data frame with results
      simJudge <- data.frame(simJudge)
      # label columns
      names(simJudge) <- c("node1", "node2", "similarity", "age")
      # get stats by age
      sed <- summarySE(simJudge, measurevar="similarity", groupvars="age")
      SS[woi,] <- sed$similarity
    }
    
    msee <- apply(EE, 2, mean)
    sdee <- apply(EE, 2, sd)
    sdee <- sdee/sqrt(nrow(EE))
    mses <- apply(SS, 2, mean)
    sdes <- apply(SS, 2, sd)
    sdes <- sdes/sqrt(nrow(SS))
   
    EEB[bis,] <- msee 
    SSB[bis,] <- mses 
  }
  
Elist[[alphai]] <- EEB
Slist[[alphai]] <- SSB
}

pdf(file="EntropySim1000.1to.5NoContext.pdf", width=9, height=6)
par(mfrow=c(1,2))
plot(1:4, Elist[[1]][1,], ylim = c(0, 1.2), cex = 0, xlab = "Epoch", ylab = "Entropy", cex.lab=1.5, xaxt="n")
axis(1, at=1:4, labels=c("1","2","3","4") )
for(i in 1:nrow(Elist[[1]])){
  lines(1:4, Elist[[1]][i,], lty = i, lwd = 1.5)
}
for(i in 1:nrow(Elist[[2]])){
  lines(1:4, Elist[[2]][i,], lty = i, col = "red", lwd=1.5)
}


plot(1:4, Slist[[1]][1,], cex = 0, ylim =c(0, 150), xlab = "Epoch", ylab = "Similarity", cex.lab = 1.5, xaxt="n")
axis(1, at=1:4, labels=c("1","2","3","4") )
for(i in 1:nrow(Slist[[1]])){
  lines(1:4, Slist[[1]][i,], lty = i, lwd=1.5)
}
for(i in 1:nrow(Slist[[2]])){
  lines(1:4, Slist[[2]][i,], lty = i, col = "red", lwd=1.5)
}
legend(2.9, 150, legend=c(TeX('0'),TeX('1')), title=TeX('a'),col = c("black", "red"), lty = 1, bty="n", lwd = 1.5, cex = .9)
legend(1, 50, legend=c(TeX('.1'), TeX('.2'), TeX('.3'),TeX('.4'), TeX('.5')), title=TeX('\\beta'), lty = 1:4, bty="n", lwd=1.4, cex = .9)
dev.off()

#save(Elist, Slist, file="WorldLearnerFrom.1to.5.RData")
