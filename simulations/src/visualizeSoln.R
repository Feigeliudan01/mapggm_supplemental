# visualizeSoln.R -- vizualize the solution from optimization & get distance
library(igraph)

vizSoln <- function(ans, true.adj, id){
  if(!is.null(ans)){
    Omega.found <- ans$Omega
    adj.found <- abs(Omega.found)>0
    diag(adj.found) <- 0
    graph.found <- graph.adjacency(adj.found, mode='undirected')
    graph.true <- graph.adjacency(true.adj, mode='undirected')
    ly <- layout.circle(graph.true)
    par(mfrow=c(2,2))
    plot(graph.true, layout=ly, edge.color = "gray50", 
         vertex.color = ifelse(duplicated(id), 'red','white'),
         vertex.size = 15, vertex.label = id, 
         main = "True Pattern")
    plot(graph.found, layout=ly, edge.color = "gray50", 
         vertex.color = ifelse(duplicated(id), 'red','white'),
         vertex.size = 15, vertex.label = id, 
         main = sprintf("Found Pattern, l=%.2f", ans$lambda))
    image(as.matrix(true.adj), main = "True Adjacency Matrix")
    image(adj.found, main = "Found Adjacency Matrix")
  } else { print('Solution not found') }
}

vizSoln3 <- function(ans1, ans2, ans3, true.adj, true.omega, id, plot=FALSE){
  omega.list <- list(ans1$Omega, ans2$Omega, ans3$Omega)
  lambda.list <- c(ans1$lambda, ans2$lambda, ans3$lambda)
  labels <- c('Structured', 'Unstructured','Separated')
  
  true.adj <- as.matrix(true.adj==1)
  graph.true <- graph.adjacency(true.adj, mode='undirected')
#   ly <- layout.kamada.kawai(graph.true)
  ly <- layout.circle(graph.true)
  
  if(plot){
    par(mfcol=c(2,4), cex.main=1)
    image(true.omega, main = "True precision",
        col=gray.colors(256))
    plot(graph.true, layout=ly, edge.color = "black", 
       vertex.color = ifelse(duplicated(id), 'SkyBlue','white'),
       vertex.size = 15, vertex.label = id, 
       main = "True Pattern")
  }

  findings <- data.frame(ind=c(0,rep(c(1:3),3)), 
                         type=c('edges', rep('found',3), rep('added',3), rep('lost',3)),
                                              val=NA)
  findings[findings$type=='edges','val'] <- sum(true.adj)/2
  for(i in 1:length(labels)){
    omega <- omega.list[[i]]
    adj <- abs(omega)>0
    diag(adj) <- FALSE
    
    found <- adj & true.adj
    added <- adj & !true.adj
    lost <- !adj& true.adj
    
    findings[findings$ind==i & findings$type=='found','val'] <- sum(found)/2
    findings[findings$ind==i & findings$type=='added','val'] <- sum(added)/2
    findings[findings$ind==i & findings$type=='lost','val'] <- sum(lost)/2

    if(plot){
    precis.text <- sprintf('%s \n Lambda= %.1f', labels[i], lambda.list[i])
    graph.text <- sprintf('Found: %d/%d \n Added: %d, Lost: %d', sum(found)/2, 
                    sum(true.adj)/2, sum(added)/2, sum(lost)/2)
    
   graph.found <- graph.adjacency(adj & true.adj, mode='undirected')
   graph.added <- graph.adjacency(adj & !true.adj, mode='undirected')   
   graph.lost <- graph.adjacency(!adj & true.adj, mode='undirected')
   
   
   image(omega,  col=gray.colors(256), main=precis.text)
   plot(graph.found, edge.color='black', edge.width=1, layout=ly,      
     vertex.color = ifelse(duplicated(id), 'SkyBlue','white'),
     vertex.size = 15, vertex.label = id, main=graph.text )
   plot(graph.added, edge.color='green', edge.width=1, add=TRUE, layout=ly,      
        vertex.color = ifelse(duplicated(id), 'SkyBlue','white'),
        vertex.size = 15, vertex.label = id)
   plot(graph.lost, edge.color='red', edge.width=1, edge.lty=3, add=TRUE, layout=ly,      
        vertex.color = ifelse(duplicated(id), 'SkyBlue','white'),
        vertex.size = 15, vertex.label = id)
    }
  }
  return(findings)
}


recoverMu <- function(Omega, Y){
  mu  <- Omega %*% colMeans(Y)
  return(mu[,1,drop=TRUE])
}

vizRecovery <- function(mu.true, mu.hat, id){
#   par(mfrow=c(1,1))
  rng <- range(c(mu.true, mu.hat))
  plot(mu.true, type='p', col=ifelse(duplicated(id), 'red','black'), pch=1, ylim=rng)
  lines(mu.hat, type='o', col=ifelse(duplicated(id), 'red','black'), pch=16)
  abline(h=0, lty=2)
}






