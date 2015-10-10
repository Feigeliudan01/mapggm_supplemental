# 2_EstimatePrecisions.R
# Estimates the network/covariance/precision built by 1_SimulateCovariances.R
#
# This is designed to be run from the command line using something like 
#   Rscript 2_EstimatePrecisions.R --args BATCH i
# where BATCH is the batch indicator, and i is the index in setup.csv currently
# being run.  Can set testing=TRUE if performing a test run not from terminal.
# This will also reduce the number of lambda's attempted from 10 to 2.
#
# This file depends on the mapggm package available on github. The first time
# you run this script, set first.run=TRUE to download the package.  After that,
# set to FALSE so that the installed package is simply loaded.

code.dir <- '~/Dropbox/phd/kolaczyk/mapggm/github/mapggm_supplemental/simulations/' # YOUR FOLDER HERE
testing <- TRUE
first.run <- FALSE

# Packages & source
if(first.run){
  library(devtools)
  install_github('paulajgriffin/mapggm')  
}

library(mapggm) 
library(Matrix)
source('src/visualizeSoln.R')

# Command line arguments
if(testing){
  args <- c('demo', 1)
  setwd(code.dir)
} else{ 
  args <- commandArgs(TRUE)
} 
batch <- args[1]
index <- as.numeric(args[2])
setup <- read.csv(sprintf('results_%s/setup/setup.csv',args[1]), as.is=TRUE)
setup <- setup[setup$index==index,]

if(setup$estimate | setup$index==1 ){
  # Labels
  setup$lambda3 <- setup$lambda2 <- setup$lambda1 <- NA
  setup$n.comp3 <- setup$n.comp2 <- setup$n.comp1 <- NA
  setup$edges0 <- NA
  setup$found3 <- setup$found2 <- setup$found1 <- NA
  setup$added3 <- setup$added2 <- setup$added1 <- NA
  setup$lost3 <- setup$lost2 <- setup$lost1 <- NA
  setup$h3.combo <- setup$h3.all <- setup$h2.combo <- setup$h2.all <- NA
  setup$h1.combo <- setup$h1.all <- NA
  setup$run.hours <- NA
  setup$sigmaHat1 <- sprintf('results_%s/storage/sigmaHat1_%d.txt', batch, setup$index)
  setup$sigmaHat2 <- sprintf('results_%s/storage/sigmaHat2_%d.txt', batch, setup$index)
  setup$sigmaHat3 <- sprintf('results_%s/storage/sigmaHat3_%d.txt', batch, setup$index)
  setup$omegaHat1 <- sprintf('results_%s/storage/omegaHat1_%d.txt', batch, setup$index)
  setup$omegaHat2 <- sprintf('results_%s/storage/omegaHat2_%d.txt', batch, setup$index)
  setup$omegaHat3 <- sprintf('results_%s/storage/omegaHat3_%d.txt', batch, setup$index)
  
  # Start the thing
  if(!testing){
    sink(sprintf('results_%s/log/log_%s.txt', batch, setup$index))
  } 
  start.time <- Sys.time()
  set.seed(setup$index)
  
  # Read in previously defined variables
  p <- 2*setup$n.paired + setup$n.unpaired
  mu.true <- diag(p)
  nullY <- as.matrix(read.table(setup$nullY))
  nullS  <- crossprod(nullY)
  id <- rep(1:(setup$n.paired + setup$n.unpaired),
            times=c(rep(2,setup$n.paired), rep(1, setup$n.unpaired)))
  trueOmega <- as.matrix(read.table(setup$omegaTrue))
  trueTheta <- Matrix(as.numeric(trueOmega!=0), sparse=TRUE, 
                      nrow=nrow(trueOmega), ncol=ncol(trueOmega))
  diag(trueTheta) <- 0
  
  # Run opt
  lambda.range1 <- getLambdaRange(S=nullS, id=id, length.out=setup$n.lambdas)
  lambda.range2 <- getLambdaRange(S=nullS, id=1:length(id), 
                                  length.out=setup$n.lambdas)
  lambda.range3 <- getLambdaRange(S=nullS[!duplicated(id),!duplicated(id)],
                                  id=unique(id), length.out=setup$n.lambdas)
  
  if(testing ){
    lambda.range1 <- lambda.range1[1:2] 
    lambda.range2 <- lambda.range2[1:2]
    lambda.range3 <- lambda.range3[1:2]
  }
  
  # EBIC-based selection
  ans1 <- multiAttSelect(S=nullS, n=setup$n, id=id, 
                         lambda.range=lambda.range1,  mode=1, 
                         max.gap=setup$max.gap, max.iter=setup$max.iter, 
                         min.t=setup$min.t, method='EBIC', 
                         Theta.true=NULL, plot=setup$label, update=5)   
  ans2 <- multiAttSelect(S=nullS, n=setup$n, id=id, 
                         lambda.range=lambda.range1,  mode=2, 
                         max.gap=setup$max.gap, max.iter=setup$max.iter, 
                         min.t=setup$min.t, method='EBIC', 
                         Theta.true=NULL, plot=setup$label, update=5)   
  ans3 <- multiAttSelect(S=nullS, n=setup$n, id=id, 
                         lambda.range=lambda.range1,  mode=3, 
                         max.gap=setup$max.gap, max.iter=setup$max.iter, 
                         min.t=setup$min.t, method='EBIC', 
                         Theta.true=NULL, plot=setup$label, update=5)  
  
  # Write out results
  write.table(ans1$Omega, setup$omegaHat1, col.names=FALSE, row.names=FALSE)
  write.table(ans2$Omega, setup$omegaHat2, col.names=FALSE, row.names=FALSE)
  write.table(ans3$Omega, setup$omegaHat3, col.names=FALSE, row.names=FALSE)
  
  write.table(ans1$Sigma, setup$sigmaHat1, col.names=FALSE, row.names=FALSE)
  write.table(ans2$Sigma, setup$sigmaHat2, col.names=FALSE, row.names=FALSE)
  write.table(ans3$Sigma, setup$sigmaHat3, col.names=FALSE, row.names=FALSE)
  
  # Track info
  setup$lambda1 <- ans1$lambda
  setup$lambda2 <- ans2$lambda
  setup$lambda3 <- ans3$lambda
  
  setup$n.comp1 <- ans1$n.comp
  setup$n.comp2 <- ans2$n.comp
  setup$n.comp3 <- ans3$n.comp
  
  graph.stats <- vizSoln3(ans1=ans1, ans2=ans2, ans3=ans3, true.adj=trueTheta,
           true.omega=trueOmega, id=id, plot=FALSE)
  
  # Optimization stats
  graph.stats$match <- paste(graph.stats$type, graph.stats$ind, sep='')
  setup[,graph.stats$match] <- graph.stats$val
  
  # Hamming distance
  hamcols <- sprintf('h%d.%s', rep(1:3,each=2), rep(c('all','combo'),3))
  setup[,hamcols] <- c(getHammingDist(ans1$Omega, trueTheta, id), 
                       getHammingDist(ans2$Omega, trueTheta, id), 
                       getHammingDist(ans3$Omega, trueTheta, id))
  
  # Write out summary file
  setup$run.hours <- difftime(Sys.time(), start.time, units='hours')
  files.there <- system(sprintf('ls results_%s/summary', batch), intern=TRUE)
  first.file <- (length(files.there)==0)
  write.table(setup, 
            sprintf('results_%s/summary/esimationSummary.csv', batch), 
            append=!first.file, col.names=first.file, sep=',', row.names=FALSE)
  save.image(sprintf('results_%s/Rdata/outdata_%s.Rdata', batch, setup$index))
  sink()
  
  print('DONE')
} else {
  print('DONE: no estimation performed')
}
