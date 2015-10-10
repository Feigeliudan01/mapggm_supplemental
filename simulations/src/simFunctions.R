require(mvtnorm)
require(Matrix)
require(igraph)
require(huge)


generateCov <- function(n.paired=10, n.unpaired=10,  plabel='', viz=TRUE, 
                       corr.in=0.6, corr.out=0.2, corr.out.off=NA,
                       diag.add=0.05, id=NULL, type='er', n.clust=1, 
                       er.prob=.2, nn.s=4, sb.probin=.2, sb.probout=.1,
                       png=NULL, min.diag.add=1){  
  # Garbage collection
#   stop()
  gcinfo(FALSE)

  # Handling off-diagonals (option added 2/6)
  if(is.na(corr.out.off)){ corr.out.off <- corr.out }
  
  # Building a gene/protein joint node set
  if ( is.null(id) ){
    d <- n.paired*2 + n.unpaired # d is total length of vector
    g <- n.paired + n.unpaired # g is total number of nodes
    id <- rep(1:g, times=c(rep(2,n.paired), rep(1, n.unpaired)))
  }
  d <- length(id) # d = total number of variables
  g <- length(unique(id)) # g = number of joint nodes
  g.list <- as.numeric(table(id)[unique(id)]) 
  g.ind <- rep(c(1:g), g.list)
  
  
  # Templates for adjacency matrices
  Theta.in <- Theta.out <- matrix(0, d, d)

  # Establish nodes
  for (i in 1:g) {
    tmp = which(g.ind == i)
    Theta.in[tmp, tmp] = 1
    rm(tmp)
    gc()
  }
  
  # Establish between-node edges
  # First: cluster assignment
#   clust.assign <- sample(rep(1:n.clust, g)[1:g], g)
  clust.assign <- sort(rep(1:n.clust, g)[1:g])
  
  # Second: within each cluster, assignment based on method
  if (type=='er') {
    # Erdos-Renyi for cross-node (but within-cluster) edges
    tmp <- matrix(runif(g^2, 0, .5), g, g)
    tmp <- ( (tmp + t(tmp)) < er.prob )
    Theta.out <- tmp[g.ind, g.ind]
    same <- ((sqrt(tcrossprod(clust.assign)) %% 1) ==0 )
    Theta.out[Theta.in==1] <- 0
    Theta.out  <- same[g.ind, g.ind] * Theta.out
    rm(tmp)
  } else if (type=='sb'){
    # Stochastic blockmodel
    subclust.assign <- rep(0, length(clust.assign))
    for(ind in 1:max(clust.assign)){
      subind <- which(clust.assign==ind)
      subclust.assign[subind[1:ceiling(length(subind)/2)]] <- 1
      subclust.assign[which(clust.assign==ind) & subclust.assign==0] <- 2
    }
    tmp <- matrix(runif(g^2, 0, .5), g, g)
    tmp1 <- ( (tmp + t(tmp)) < sb.probin ) # within-block
    tmp2 <- ( (tmp + t(tmp)) < sb.probout ) # out-block
    Theta1.out <- tmp1[g.ind, g.ind] 
    Theta2.out <- tmp2[g.ind, g.ind]
    same <- ((sqrt(tcrossprod(clust.assign)) %% 1)  ==0 )
    same.subin <- ((sqrt(tcrossprod(subclust.assign)) %% 1) ==0 ) * same
    same.subout <- ((sqrt(tcrossprod(subclust.assign)) %% 1) !=0 ) * same
    Theta.out  <- (same.subin[g.ind, g.ind] * Theta1.out) + 
      (same.subout[g.ind, g.ind] * Theta2.out)                                
    Theta.out[Theta.in==1] <- 0                 
    rm(tmp)               
  }else if (type=='chain'){
    # Chain simulation as described in Kolar et al
    tmp = matrix(0, g, g)
    for(clust in 1:n.clust){
       chain <- sample(which(clust.assign==clust), sum(clust.assign==clust))
       for(ind in 1:(length(chain)-1)){
         tmp[chain[ind], chain[ind+1] ] <- tmp[chain[ind+1], chain[ind] ] <- 1
       }
    }
    Theta.out <- tmp[g.ind, g.ind]
  } else if (type=='nn'){
    # Nearest-neighbor as described in Kolar et al
    tmp = matrix(0, g, g)
    coords <- cbind(runif(g), runif(g))
    for(clust in 1:n.clust){
      inclust <- which(clust.assign==clust)
      dists <- as.matrix(dist(coords[inclust,]))
      tmp[inclust, inclust] <- apply(dists, 2, function(vec){rank(abs(vec))<=nn.s})
    }
    tmp = ceiling((tmp + t(tmp))/2)
    Theta.out <- tmp[g.ind, g.ind]
    Theta.out[Theta.in==1] <- 0
  } else {
    stop('Error: invalid type') 
  }
  gc()
#   stop()
  # Build a precision matrix
  Theta <- Theta.in + Theta.out
  diag(Theta) <- 0
  Theta.out.off <- Theta.out.on <- matrix(0, nrow(Theta), ncol(Theta))
  dup <- duplicated(id)
  Theta.out.on[dup,dup] <- Theta.out[dup,dup]
  Theta.out.on[!dup,!dup] <- Theta.out[!dup,!dup]
  Theta.out.off[dup,!dup ] <- Theta.out[dup,!dup]
  Theta.out.off[!dup,dup] <- Theta.out[!dup,dup]
  Omega <- -1*(Theta.in*corr.in + Theta.out.on*corr.out +
                 Theta.out.off*corr.out.off)
  diag(Omega) <- 1
#   browser()
  # Tweak to ensure positive-definite
  mineval <- min(eigen(Omega)$values)
  add.counter <- 0
  while( mineval < 0.5 | add.counter < min.diag.add){
    diag(Omega) =   diag(Omega) + rexp(nrow(Omega)) * diag.add
#     print(diag(Omega)[1])
    mineval <- min(eigen(Omega)$values)
    add.counter <- add.counter + 1
  }
  Omega <- Omega/mean(diag(Omega))
  Sigma <- solve(Omega)
#   Sigma <- solve(Omega)
#   sc <- mean(diag(Sigma))
#   Sigma <- Sigma/sc
#   Omega <- Omega*sc

  
  if (viz == TRUE) {
    par(cex.main=1.5)
#     if(!is.null(png)){png(sprintf('%s.png',label))}
    fullfig = par(mfrow = c(2,2), pty = "s", 
                  omi = c(0.3, 0.3, 0.3, 0.3), mai = c(0.3, 0.3, 0.3, 0.3))
        g = graph.adjacency(Theta, mode = "undirected", diag = FALSE)
    layout.grid = layout.circle(g)
    fullfig[1] = plot(g, layout = layout.grid, edge.color = "black", 
                      vertex.color = ifelse(duplicated(id), 'SkyBlue','white'),
                      vertex.size = 15, vertex.label = id, 
                      main = plabel)
    fullfig[2] = image(Theta, col = gray.colors(256), main = "Adjacency Matrix")

    fullfig[3] = image(Omega, col = gray.colors(256), main = "Precision Matrix")
    fullfig[4] = image(Sigma, col = gray.colors(256), main = "Covariance Matrix")

#     if(!is.null(png)){dev.off()}
    rm(fullfig, g, layout.grid)
    gc()
   }

  out = list(Sigma=Sigma, Omega=Omega, Theta=Matrix(Theta, sparse = TRUE), id=id,
             add.counter=add.counter)
  return(out)
}

modifyCov <- function(Omega.start, id, corr.in.start, corr.out.start,
                      corr.in=0.3, corr.out=0.1, corr.out.off=NA, diag.add=.05,
                      min.diag.add=0){
  # Handling off-diagonals (option added 2/6)
  if(is.na(corr.out.off)){ corr.out.off <- corr.out }
  
  if(corr.in.start!=0){
    sc <- -corr.in.start/Omega.start[1,2]
  } else if (corr.out.start!=0){
    tmp <- Omega.start
    diag(tmp) <-0
    sc <- -corr.out.start/tmp[tmp!=0][1]
  } else {
    Theta.start <- (Omega.start==0)
    diag(Theta.start) <- 0
    Sigma.start <- solve(Omega.start)
    return(list(Sigma=Sigma.start, Omega=Omega.start, 
                Theta=Matrix(Theta.start, sparse = TRUE), id=id, add.counter=0))
  }
  Omega.sc <- Omega.start*sc
  
  d <- length(id) # d = total number of variables
  g <- length(unique(id)) # g = number of joint nodes
  g.list <- as.numeric(table(id)[unique(id)]) 
  g.ind <- rep(c(1:g), g.list)
  
  
  # Templates for adjacency matrices
  Theta.in <- Theta.out <- matrix(0, d, d)

  # Revise corr.in entries
  for (i in 1:g) {
    tmp = which(g.ind == i)
    Theta.in[tmp, tmp] <- 1
    rm(tmp)
    gc()
  }  
  
  # Revise corr.out entries
  Theta.out <- Omega.start!=0
  Theta.out[which(Theta.in!=0)] <- 0
  Theta <- Theta.in + Theta.out
  diag(Theta) <- 0
  Theta.out.off <- Theta.out.on <- matrix(0, nrow(Theta), ncol(Theta))
  dup <- duplicated(id)
  Theta.out.on[dup,dup] <- Theta.out[dup,dup]
  Theta.out.on[!dup,!dup] <- Theta.out[!dup,!dup]
  Theta.out.off[dup,!dup ] <- Theta.out[dup,!dup]
  Theta.out.off[!dup,dup] <- Theta.out[!dup,dup]
  Omega <- -1*(Theta.in*corr.in + Theta.out.on*corr.out +
                 Theta.out.off*corr.out.off)
  diag(Omega) <- diag(Omega.sc)
  mineval <- min(eigen(Omega)$values)
  add.counter <- 0
  while( mineval < 0.5 | add.counter < min.diag.add){
    diag(Omega) =   diag(Omega) + rexp(nrow(Omega)) * diag.add
#     print(diag(Omega)[1])
    mineval <- min(eigen(Omega)$values)
    add.counter <- add.counter + 1
  }
  if(add.counter>0){print('warning: adding to diag')}
  Omega <- Omega/mean(diag(Omega))
  Sigma <- solve(Omega)
#   Sigma <- solve(Omega)
#   sc2 <- mean(diag(Sigma))
#   Sigma <- Sigma/sc2
#   Omega <- Omega*sc2
  
  out = list(Sigma=Sigma, Omega=Omega, Theta=Matrix(Theta, sparse = TRUE), id=id,
             add.counter=add.counter)
  return(out)
}


generateData <- function(Sigma, n=100, snr=NA, perturb.loc=NA, vis=FALSE, 
                         returnCov=TRUE, png=NULL){
  # Number of variables
  p <- nrow(Sigma)
  
  # Error-checking
  if((!is.na(sum(perturb.loc)) & mean(perturb.loc %in% 1:p) < 1) | 
      (length(snr) != length(perturb.loc)) ){
    stop('Error: invalid perturbation/SNR specified') 
  } else if (is.na(sum(perturb.loc))) {
    print('No perturbation or SNR specified, will simulate control data')
    perturb.loc <- 1
    snr <- 0
  }
  
  # Set up perturbation vector
  mu <- rep(0, p)
  mu[perturb.loc] <- snr
  
  # Generate data
  Y <- rmvnorm(n, as.numeric(Sigma %*% mu), Sigma)
  Sigma.hat <- crossprod(Y)/n
  
  
  if (vis == TRUE) {
    fullfig = par(mfrow = c(1, 2), pty = "s", omi = c(0.3, 0.3, 0.3, 0.3),
                  mai = c(0.3, 0.3, 0.3, 0.3))
    fullfig[1] = image(Sigma, col=gray.colors(256), main="Covariance Matrix")
    fullfig[2] = image(Sigma.hat, col=gray.colors(256), main="Empirical Matrix")
    rm(fullfig)
    gc()
   }
  if(returnCov){
   out <- list(S=crossprod(Y)/n, n=n, mu=mu) 
  } else {
    out <- list(Y=Y, mu=mu)
  }
  return(out)
}


dropCovEdge <- function(Omega.truth, id){
  Theta.truth <- blockNorms(M=Omega.truth, id=id)!=0
  diag(Theta.truth) <- 0
  Theta.truth[!upper.tri(Theta.truth)] <- 0
  print(length(which(Theta.truth==1)))
  if(length(which(Theta.truth==1))>1){
    killEdge <- sample(which(Theta.truth==1),1)
    Theta.truth[killEdge] <- 0
  } else if (length(which(Theta.truth==1))==1){
    killEdge <- which(Theta.truth==1)
    Theta.truth[killEdge] <- 0
  }
  Theta.truth <- Theta.truth + t(Theta.truth)
  diag(Theta.truth) <- 1
  bigTheta.truth <- Theta.truth[id,id]
  newOmega <- Omega.truth * bigTheta.truth
}
