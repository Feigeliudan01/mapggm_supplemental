require(ROCR)
#require(RColorBrewer)
#require(ggplot2)
#require(gridExtra)
require(Hotelling)

TidyPerformance <- function(perf, label = NULL) {
  thresholds <- unlist(perf@alpha.values)
  thresholdRange <- range(thresholds[is.finite(thresholds)])
  thresholds <- rev(seq(thresholdRange[1], thresholdRange[2],
                        length = max(sapply(perf@alpha.values, length))))
  metrics <- list()
  for (i in 1:length(perf@x.values)) {
    x <- approxfun(perf@alpha.values[[i]], perf@x.values[[i]],
                   rule=2, ties=mean)(thresholds)
    y <- approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                   rule=2, ties=mean)(thresholds)
    metrics[[i]] <- data.table(x = x, y = y, threshold = thresholds)
  }
  metrics <- rbindlist(metrics)
  metrics <- metrics[, lapply(.SD, mean), by = threshold]
  metrics[, label := label]
  return(metrics)
}


getPerformance <- function(mu.list, truth.list){
  roclist <- list(1:length(mu.list))
  for(i in 1:length(mu.list)){
    # Performance calculations
      perf <- prediction(as.vector(apply(abs(mu.list[[i]]),2,rank)), 
                         as.vector(truth.list[[i]]))
    roclist[[i]] <- performance(perf, 'tpr', 'fpr')
  }
  return(roclist)
}


aucList <- function(mu.list, truth.list, method.labels, r.in=NA, r.out=NA, 
                    use.stats=FALSE){
  roclist <- list()
  aucs <- rep(NA, length(mu.list))  
  for(i in 1:length(mu.list)){
    # Performance calculations
    if(use.stats){
      perf <- prediction( as.vector(mu.list[[i]]), as.vector(truth.list[[i]]))
    } else {
      perf <- prediction(as.vector(apply(abs(mu.list[[i]]),2,rank)), 
                         as.vector(truth.list[[i]]))
    }
    roclist[[i]] <- performance(perf, 'tpr', 'fpr')
    aucs[i] <-  mean(unlist(performance(perf , 'auc')@y.values))
  }
  return(aucs)
}


testDE <- function(Y.perturb, Y.null){
  pval <- sapply(1:ncol(Y.perturb), 
                 function(i){ 
                   t.res <- t.test(Y.perturb[,i], Y.null[,i])
                   t.res$p.value
                 })
  return(pval) 
}


testDET2 <- function(Y.perturb, Y.null, perturb.mat, ret='pval'){
  out <- apply(perturb.mat, 2, 
                function(perturb.vec){
                   inds <- which(perturb.vec)
                   t.res <- hotelling.test(x=Y.perturb[,inds], y=Y.null[,inds])
                   if(ret=='pval'){
                      return(t.res$pval)
                   } else {
                     return(t.res$stats$statistic) 
                   }
                 })
  return(out) 
}


topNProb <- function(mu.list, truth.list1, method.labels, n){

  props <- rep(NA, length(mu.list))
  
  for (i in 1:length(mu.list)){ 
    ranks <- as.vector(apply(-abs(mu.list[[i]]),2,rank))
    truth <- as.logical(as.vector(truth.list1[[i]]))
    ranks.keep <- ranks[which(truth)]
    props[i] <- round(mean(ranks.keep<=n),2)
  }
  
  names(props) <- method.labels
  return(props)
}