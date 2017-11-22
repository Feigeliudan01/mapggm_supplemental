# 1_SimulateCovariances.R
# Simulates networks and control data based on the design set up in 
# the 0_DesignSims.R script.  
#
# This is designed to be run from the command line (simulations directory) using 
#   Rscript 1_SimulateCovariances.R --args BATCH j 
# where BATCH is the batch indicator, and j is the replicate number being run.
# This is designed to split the simulation into manageable chunks, but there's
# no particular importance in that design.
# Generated data and networks are sent to the 'storage' subdirectory in the
# batch-specific results folder.

code.dir <- '/restricted/projectnb/johnsonlab/yuqingz/gaussian_network/DC_SBM_yz/mapggm_supplemental/DCSBM/'  # YOUR FOLDER HERE
testing <- FALSE # set to TRUE if you're not running this from the command line
source('src/simFunctions.R')
source('src/simFunctions_dcsbm.R')

if(testing){
  args <- c('demo', 1)
  setwd(code.dir)
} else{ 
  args <- commandArgs(TRUE)
} 
batch <- args[1]
rep <- as.numeric(args[2])

# Read in the setup file created by 0_DesignSims.R
setup <- read.csv(sprintf('results_%s/setup/setup.csv', batch), as.is=TRUE)
setup$ab <- paste(abs(setup$corr.in), abs(setup$corr.out), setup$sb.probin, setup$sb.probout, setup$type, setup$n.paired)
setup <- setup[(setup$each == rep) , ]
sink(sprintf('results_%s/log/logSimCov_rep%d.log', batch, rep))

# Start generating networks
start.time <- Sys.time()
cat('\nGenerating covariances with similar diagonal terms')
for(rep in unique(setup$each)){
  print(rep)
  for(ab in unique(setup$ab[setup$each==rep])){
    # Fill in for (+,+) line
    rin.start <- max(setup$corr.in[setup$each==rep & setup$ab==ab])
    rout.start <- max(setup$corr.out[setup$each==rep & setup$ab==ab])
    line <- which(setup$each==rep & setup$ab==ab &
                    setup$corr.in==rin.start & setup$corr.out==rout.start)
    if(length(line)==1){
      simcov.start <- generateCov_new(n.paired=setup$n.paired[line], 
                                      n.unpaired=setup$n.unpaired[line], 
                                      corr.in=setup$corr.in[line], 
                                      corr.out=setup$corr.out[line], 
                                      corr.out.off=setup$corr.out.off[line], 
                                      type=setup$type[line], 
                                      n.clust=setup$n.clust[line], 
                                      er.prob=setup$er.prob[line], 
                                      sb.probin=setup$sb.probin[line],
                                      sb.probout=setup$sb.probout[line],
                                      nn.s=setup$s.nn[line],
                                      dcsbm.theta=setup$dcsbm.theta[line],
                                      min.diag.add=1, viz=FALSE)
    } else { stop(sprintf('error: rep %s ab %s', rep, ab)) }
    
#     # Allow perturbed control data -- only use if you're specifically testing!
#     if(!is.na(setup$n.perturbed)){
#       if(setup$n.perturbed>0){
#         perturb.loc <- which(simcov$id<setup$n.perturbed) 
#       } else{
#         perturb.loc <- NA
#       }
#     }
  
    perturb.loc <- NA
    nulldata <- generateData(Sigma=simcov.start$Sigma, 
                             n=setup$n[line], snr=0, 
                             perturb.loc=perturb.loc, 
                             vis=FALSE, returnCov=FALSE, png=NULL)
    
    write.table(simcov.start$Sigma, file=setup$sigmaTrue[line], row.names=FALSE, col.names=FALSE)
    write.table(simcov.start$Omega, file=setup$omegaTrue[line], row.names=FALSE, col.names=FALSE)
    write.table(nulldata$Y, file=setup$nullY[line], row.names=FALSE, col.names=FALSE)
    
    # Fill in for (+,-) line
    rin.now   <- rin.start
    rout.now  <- -rout.start
    line <- which(setup$each==rep & setup$ab==ab &
                    setup$corr.in==rin.now & setup$corr.out==rout.now)
    if(length(line)==1){
      simcov.new <- modifyCov(Omega.start=simcov.start$Omega, id=simcov.start$id, 
                              corr.in.start=rin.start, corr.out.start=rout.start,
                              corr.in = rin.now, corr.out=rout.now, corr.out.off=NA)
      nulldata <- generateData(Sigma=simcov.new$Sigma, 
                             n=setup$n[line], snr=0, perturb.loc=NA, 
                             vis=FALSE, returnCov=FALSE, png=NULL)
    
      write.table(simcov.new$Sigma, file=setup$sigmaTrue[line], row.names=FALSE, col.names=FALSE)
      write.table(simcov.new$Omega, file=setup$omegaTrue[line], row.names=FALSE, col.names=FALSE)
      write.table(nulldata$Y, file=setup$nullY[line], row.names=FALSE, col.names=FALSE)
    } else if (length(line)>1){ stop(sprintf('error: rep %s ab %s', rep, ab)) }


    # Fill in for (-,+) line
    rin.now <- -rin.start
    rout.now <- rout.start
    line <- which(setup$each==rep & setup$ab==ab & 
                    setup$corr.in==rin.now & setup$corr.out==rout.now)
    if(length(line)==1){
      simcov.new <- modifyCov(Omega.start=simcov.start$Omega, id=simcov.start$id, 
                              corr.in.start=rin.start, corr.out.start=rout.start,
                              corr.in = rin.now, corr.out=rout.now, corr.out.off=NA)
      nulldata <- generateData(Sigma=simcov.new$Sigma, 
                       n=setup$n[line], snr=0, perturb.loc=NA, 
                       vis=FALSE, returnCov=FALSE, png=NULL)
      
      write.table(simcov.new$Sigma, file=setup$sigmaTrue[line], row.names=FALSE, col.names=FALSE)
      write.table(simcov.new$Omega, file=setup$omegaTrue[line], row.names=FALSE, col.names=FALSE)
      write.table(nulldata$Y, file=setup$nullY[line], row.names=FALSE, col.names=FALSE)
    } else if (length(line)>1){ stop(sprintf('error: rep %s ab %s', rep, ab)) }
    
    # Fill in for (-,-) line
    rin.now <- -rin.start
    rout.now <- -rout.start
    line <- which(setup$each==rep & setup$ab==ab & 
                    setup$corr.in==rin.now & setup$corr.out==rout.now)
    if(length(line)==1){
      simcov.new <- modifyCov(Omega.start=simcov.start$Omega, id=simcov.start$id, 
                              corr.in.start=rin.start, corr.out.start=rout.start,
                              corr.in = rin.now, corr.out=rout.now, corr.out.off=NA)
      nulldata <- generateData(Sigma=simcov.new$Sigma, 
                               n=setup$n[line], snr=0, perturb.loc=NA, 
                               vis=FALSE, returnCov=FALSE, png=NULL)
      
      write.table(simcov.new$Sigma, file=setup$sigmaTrue[line], row.names=FALSE, col.names=FALSE)
      write.table(simcov.new$Omega, file=setup$omegaTrue[line], row.names=FALSE, col.names=FALSE)
      write.table(nulldata$Y, file=setup$nullY[line], row.names=FALSE, col.names=FALSE)
    } else if (length(line)>1){ stop(sprintf('error: rep %s ab %s', rep, ab)) }
  } # close ab within rep
} # close rep
end.time <- Sys.time()
dt <- difftime(Sys.time(), start.time, units='secs')
print(dt)
sink()
