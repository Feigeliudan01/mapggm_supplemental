# 3_EvaluateSims.R
# Calculates performance statistics and generates plots. 
#
# testing=TRUE is just for script editing purposes -- it will only evaluate a 
# subest of runs, and it won't use sink() to generate a log file. 

# Run settings
code.dir <- '~/Dropbox/phd/kolaczyk/mapggm/github/mapggm_supplemental/simulations/'  # YOUR FOLDER HERE
batch <- 'demo'
testing <- FALSE
useIdeal <- FALSE
perturb <- 'node' # either 'node' or 'single'
test <- 'node' # either 'node' or 'single'
r.in.use <- c(.8, .6)#seq(-.8,.8, .2)#c(.8,.6)#c(-.8,.8)##
r.out.use <- c(.2, .4, .6)#c(-.4,0,.4)#
n.perturb <- 1 #2
snr <- 0.25 # c(.3,.2)

edgeDrop <- FALSE # deprecated option, do not change
condenseOR <- FALSE # deprecated option, do not change
snr.print <- paste(snr, collapse='-')

# Source calls
setwd(code.dir)
source('src/simFunctions.R')
source('src/performance.R')
library(mapggm)

setwd('~/Dropbox/phd/kolaczyk/mapggm/github/mapggm_supplemental/simulations')
dir.create(sprintf('results_%s/post', batch), showWarnings=FALSE)


if(useIdeal){
  setup <- read.csv(sprintf('results_%s/setup/setup.csv', batch), as.is=TRUE) 
  storageList <- system(sprintf('ls results_%s/storage/omegaTrue*', batch), intern=TRUE)
  keep <- setup$omegaTrue %in% storageList
  setup <- setup[keep,]
} else {
  setup <- read.csv(sprintf('results_%s/summary/estimationSummary.csv',batch), as.is=TRUE)
  setup <- setup[setup$estimate,]
}
setup <- setup[setup$corr.in %in% r.in.use,]
setup <- setup[setup$corr.out %in% r.out.use,]
setup$facet <- sprintf('%s_n%d_npair%d_sbpo%d',setup$type, setup$n, setup$n.paired, 100*setup$sb.probout)
descr <- sprintf('%s%s_snr%s_perturb%s_test%s_rank%s_cov%s_%d', ifelse(testing, "TEST_", ""),
                 ifelse(useIdeal, 'IDEAL', 'SIM'), snr.print, toupper(perturb), toupper(test), 
                 ifelse(condenseOR, 'OR', 'AND'), ifelse(edgeDrop, 'DROP','SAME'), as.integer(Sys.time()))
print(descr)

newdir <- sprintf('results_%s/post/%s', batch, descr)
dir.create(newdir, showWarnings=FALSE)

if(!testing){
  sink(sprintf('%s/3log_EvaluateSims.log', newdir))
}

fac <- setup$facet[1]
if(testing) { setup <- setup[setup$facet==fac,] } 
for (fac in unique(setup$facet)){
  print(sprintf('facet %s at %s', fac, Sys.time()))
  r.ins <- sort(unique(setup$corr.in[setup$facet==fac]), decreasing=TRUE)
  r.outs <- sort(unique(setup$corr.out[setup$facet==fac]), decreasing=FALSE)
  if(testing){
   r.ins <- r.ins[1:2]
   r.outs <- r.outs[1:2]
  }
  res.list   <- truth.list <- as.list(1:(length(r.ins)*length(r.outs)))
  auc.df <- data.frame(ind=1:(length(r.ins)*length(r.outs)),
                       r.in=rep(r.ins, each=length(r.outs)),
                       r.out=rep(r.outs, length(r.ins)),
                       lrt=NA, lrt.sep=NA, lrt1=NA, 
                       det2=NA, de1=NA, de.sep=NA, ssem1=NA, cond=NA)
                      
  plt.ind <- 0
  rin <- r.ins[1] ; rout <- r.outs[1]; rep <- 1
  for (rin in r.ins){
    for(rout in r.outs){
      print(sprintf('rout %.2f rin %.2f at %s', rout, rin, Sys.time()))
      plt.ind <- plt.ind + 1
      reps <- sort(setup$each[setup$corr.in==rin & 
                                setup$corr.out==rout & setup$facet==fac ])
      if(testing) { reps <- reps[1:2]; rep <- reps[1]} 
      print(reps)
      for(rep in reps){
        print(sprintf('rep %d  at %s', rep, Sys.time()))
        line <- which(setup$corr.in==rin & setup$corr.out==rout & 
                        setup$each==rep & setup$facet==fac )
        id <- rep(1:(setup$n.paired[line] + setup$n.unpaired[line]),
            times=c(rep(2,setup$n.paired[line]), rep(1, setup$n.unpaired[line])))
        ind1 <- which(!duplicated(id))
        ind2 <- which(duplicated(id))
        # Pull from external files
        if(useIdeal){
          # Read in cov, precision, and data
          Omega.truth <- as.matrix(read.table(setup$omegaTrue[line]))
          Sigma.truth <- as.matrix(read.table(setup$sigmaTrue[line]))
          Sigma.sim <- Sigma.truth
          nullY <- as.matrix(read.table(setup$nullY[line]))
          
          # Generate "separated"
          dup <- duplicated(id)
          Omega.sep <- Omega.truth
          Omega.sep[dup,!dup] <- 0
          Omega.sep[!dup, dup] <- 0
          Sigma.sep <- solve(Omega.sep)
          Sigma.sep <- (Sigma.sep + t(Sigma.sep))/2
          
        } else {
          # Read in cov, precision, and data
          Omega.truth <- as.matrix(read.table(setup$omegaHat1[line]))
          Sigma.truth <- as.matrix(read.table(setup$sigmaHat1[line]))
          Sigma.sim <- as.matrix(read.table(setup$sigmaTrue[line]))
          nullY <- as.matrix(read.table(setup$nullY[line]))
          
          Omega.sep <- as.matrix(read.table(setup$omegaHat2[line]))
          Sigma.sep <- as.matrix(read.table(setup$sigmaHat2[line]))          
        } # end if/else useIdeal
        
        # Set up perturbation and testing matrices
        mu.test1 <- sapply(id, function(i)(id==i))
        mu.test1 <- mu.test1[,!duplicated(t(mu.test1))]
        mu.true1 <- mu.test1
        stat.truth <- diag(length(ind1))
        
        # If n.perturb>1, need to select additional genes to perturb
        if(n.perturb>1){
          l <- 1
          for(l in 1:ncol(mu.true1)){
            perturbed.now <- unique(id[mu.true1[,l]==1])
            perturbed.other <- max(id)+1 - perturbed.now #sample(unique(id[mu.true1[,l]!=1]), n.perturb-1)
            mu.true1[id==perturbed.other,l] <- TRUE
          }
        }
        stat.truth <- (1+((mu.true1[ind1,]==1)-1))
          
        
        # Set up empty matrices for results 
        ssem1 <- de1 <- det2 <- ma.lrt <- cond.lrt <- 
          sep.lrt <- lrt1 <- matrix(NA, nrow=ncol(mu.test1), ncol=ncol(mu.test1))
        if (!condenseOR){
          sep.lrt <- matrix(NA, nrow=ncol(Omega.truth), ncol=ncol(mu.test1)) 
        }
        de.sep <- matrix(NA, nrow=ncol(Omega.truth), ncol=ncol(mu.test1)) 

        j <- 1
        for( j in 1:ncol(mu.true1)){
          # Generate perturbed data
          simdata <- generateData(Sigma=Sigma.sim, n=setup$n[line], 
                                  snr=rep(sample(snr,n.perturb), each=2), 
                                  perturb.loc=which(mu.true1[,j]==1), 
                                  vis=FALSE, returnCov=FALSE)
          # Store test statistics for all methods
          ma.lrt[,j] <- perturbNodeTests(Y=simdata$Y, Omega=Omega.truth,
                                  Sigma=Sigma.truth, id=id, sequential=FALSE, 
                                  return.value='stat')  
          cond.lrt[,j] <- perturbNodeTests(Y=simdata$Y, Omega=Omega.truth,
                                  Sigma=Sigma.truth, id=id, sequential=TRUE, 
                                  return.value='stat')                
          sep.lrt[,j] <- perturbTests(Y=simdata$Y, Omega=Omega.sep,
                                  Sigma=Sigma.sep, perturb.mat=diag(length(id)),
                                  return.value='stat')
          lrt1[,j] <- perturbTests(Y=simdata$Y[,ind1], Omega=Omega.sep[ind1,ind1],
                                Sigma=Sigma.sep[ind1,ind1], 
                                perturb.mat=mu.test1[ind1,],
                                return.value='stat')
          ssem1[,j] <- abs(rowMeans(Omega.sep[ind1,ind1] %*% t(simdata$Y[,ind1])))
          de1[,j] <- -log10(testDE(Y.perturb=simdata$Y[,ind1], Y.null=nullY[,ind1]))
          det2[,j] <- -log10(testDET2(Y.perturb=simdata$Y, Y.null=nullY, perturb.mat=mu.test1))
          de.sep[,j] <- -log10(testDE(Y.perturb=simdata$Y, Y.null=nullY))

        } # end for j in 1:ncol(mu.true1)
        
        # Fill into final results
        if(rep==reps[1]){
          # If first, create new mu vectors 
          true.full <- mu.true1
          test.full <- mu.test1
          truestat.full <- stat.truth
          pick.full <- diag(ncol(mu.test1))
          ma.lrt.full <- ma.lrt
          cond.lrt.full <- cond.lrt
          sep.lrt.full <- sep.lrt
          lrt1.full <- lrt1
          de1.full <- de1
          det2.full <- det2
          ssem1.full <- ssem1
          de.sep.full <- de.sep
        } else {
          # Otherwise, append
          true.full <- cbind(true.full, mu.true1)
          test.full <- cbind(test.full, mu.test1)
          truestat.full <- cbind(truestat.full, stat.truth)
          pick.full <- cbind(pick.full, diag(ncol(mu.test1)))
          ma.lrt.full <- cbind(ma.lrt.full, ma.lrt)
          cond.lrt.full <- cbind(cond.lrt.full, cond.lrt)
          sep.lrt.full <- cbind(sep.lrt.full, sep.lrt)
          lrt1.full <- cbind(lrt1.full, lrt1)
          de1.full <- cbind(de1.full, de1)
          det2.full <- cbind(det2.full, det2)
          ssem1.full <- cbind(ssem1.full, ssem1)
          de.sep.full <- cbind(de.sep.full, de.sep)
        } # end if/else for i==1
      } # end for rep in reps
    
      # Build list for storing answers
      plt.ind <- which(auc.df$r.in==rin & auc.df$r.out==rout)
      if(condenseOR){
        truth.list[[plt.ind]] <-list(truestat.full, truestat.full, truestat.full, truestat.full,
                                      truestat.full, truestat.full, truestat.full , truestat.full)
        res.list[[plt.ind]] <- list(ma.lrt.full, condenseResultsMatrix(sep.lrt.full, id),  lrt1.full,  det2.full, 
                                    de1.full, condenseResultsMatrix(de.sep.full, id), ssem1.full, cond.lrt.full)
      } else {
        sep.truth <- matrix(as.numeric(test.full), nrow=nrow(test.full), ncol=ncol(test.full))
        truth.list[[plt.ind]] <- list(truestat.full, sep.truth, truestat.full, truestat.full,
                                      truestat.full, sep.truth, truestat.full, truestat.full)
        res.list[[plt.ind]] <- list(ma.lrt.full, sep.lrt.full,  lrt1.full,  det2.full, 
                                  de1.full, de.sep.full, ssem1.full, cond.lrt.full)
      }
    } # end for rin in r.ins
  } # end for rout in r.outs

  save(res.list, truth.list, auc.df, snr,
       file=sprintf('%s/rd_%s_%s.Rdata', newdir, fac, ifelse(useIdeal, 'ideal','est')))

} # end for fac in unique(setup$facets)

print('DONE')
sink()
print(descr)

