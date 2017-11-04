# O_DesignSims.R
# Takes arguments for how many simulation settings and writes out a csv 
# to direct the rest of the pipeline.  When more than one value is specified
# for each of these, all possible combinations are output.  
#
# Setup parameters:
# batch = label to identify this set of simulations (recommended: date)
# n.replicates = how many replicates to run of each setting
# n.perclust = number of nodes in each stochastic blockmodel cluster
# n.paired = number of nodes with k=2 attributes
# n.unpaired = number of nodes with k=1 attributes
# n = sample size for control
# n.perturbed = sample size for perturbed data
# n.clust = number of clusters in stochastic blockmodel
# er.prob = probability of cross-node links in Erdos-Renyi model
# sb.probin = within-cluster probability of cross-node links in stochastic blockmodel
# sb.probout = cross-cluster probability of cross-node links in stochastic blockmodel
# s.nn = parameters for sample size based Kolar et al (2014), section 5
# corr.in = negative partial correlation for within-node entries in Omega
# corr.in = negative partial correlation for cross-node entries in Omega
# type = stochastic blockmodel ('sb'), Erdos-Renyi ('er'),  or KNN-chain ('chain')
# n.lambdas = number of lambdas to try in estimation of Omega
# max.iter = maximum number of iterations per component in estimation of Omega
# min.t = minimum step size for optimization in estimation of Omega
# max.gap = maximum primal-dual gap optimzation in estimation of Omega

code.dir <- '~/Dropbox/phd/kolaczyk/mapggm/github/mapggm_supplemental/simulations/'  # YOUR FOLDER HERE
batch <- 'demo'
n.replicates <- 2
n.perclust <- 10

setwd(code.dir)
build <- expand.grid(index=1,
                     each=1:n.replicates,
                     theta=NA,
                     n.paired=20,
                     n.unpaired=0,
                     n=50,
                     n.perturbed=50,
                     n.clust=1,
                     er.prob=.1,
                     sb.probin=.4,
                     sb.probout=.2,
                     s.nn=4,
                     corr.in=c(.6,.8),
                     corr.out=c(.2,.4,.6),
                     corr.out.off=NA,
                     snr=c(.25),
                     type=c('sb'),
                     n.lambdas=10,
                     max.iter=200,
                     min.t=2e-30,
                     max.gap=0.01)

# Corrections that didn't need to go factorial
build$n.clust <- ceiling(build$n.paired/n.perclust)
build$corr.out <- round(build$corr.out, 2)
build$corr.in <- round(build$corr.in, 2)

# Set below to FALSE if the pipeline should be run with true (not estimated) Omega
build$estimate <- TRUE 

# Removing rows when no network is present
build<- build[!(build$type!='er' & build$corr.out==0),] 
build$index <- 1:nrow(build)

# Calibrating sample size based on number of nearest neighbors (Kolar section 5)
build$s.nn[build$type=='er'] <- build$er.prob[build$type=='er']*10
build$s.nn[build$type=='chain'] <- 2

# Add columns to maintain a record of where results are stored
build$parentdir <- sprintf('results_%s/setup', batch)
build$sigmaTrue <- sprintf('results_%s/storage/sigmaTrue_%s.txt', batch, build$index)
build$omegaTrue <- sprintf('results_%s/storage/omegaTrue_%s.txt', batch, build$index)
build$nullY <- sprintf('results_%s/storage/nullY_%s.txt', batch, build$index)

# Create directory for results
dir.create(sprintf('results_%s/', batch))

# Structures within results directory
dir.create(sprintf('results_%s/setup', batch))
dir.create(sprintf('results_%s/summary', batch))
dir.create(sprintf('results_%s/log', batch))
dir.create(sprintf('results_%s/plot', batch))
dir.create(sprintf('results_%s/Rdata', batch))
dir.create(sprintf('results_%s/storage', batch))

# Write out build info csv
write.csv(build, sprintf('results_%s/setup/setup.csv', batch), row.names=FALSE)