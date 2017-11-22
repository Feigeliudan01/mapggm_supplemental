# 4_PlotSims.R
# Plots the results of 3_EvaluateSims.R.
# Supply a batch indicator,  case data tag (the subdirectory in post),
# and the name of the Rdata file that you want to use.

# Run settings
code.dir <- '/restricted/projectnb/johnsonlab/yuqingz/gaussian_network/DC_SBM_yz/mapggm_supplemental/DCSBM/'  # YOUR FOLDER HERE
batch <- 'sim_DCSBM'
tag <- "SIM_snr0.2_perturbNODE_testNODE_rankAND_covSAME_1510028066"
rdfile <- 'rd_dcsbm_n50_npair20_sbpo20_est.Rdata'

# Source calls
setwd(code.dir)
library(ggplot2)
library(data.table)
library(xtable)
source('src/performance.R')

# Getting ready
post_dir <- sprintf('results_%s/post/%s/', batch, tag)
load(sprintf('results_%s/post/%s/%s', batch, tag, rdfile))
n.perturb <- 1
fac <- strsplit(strsplit(rdfile,'rd_')[[1]][2],'.Rdata')[[1]][1]
ideal <- strsplit(tag,'_')[[1]][1]


# Plotting settings
method.labels <- c('Multi-att. NF', 'Separated NF', 'Single-att. NF', 
                   'Multi-att. diff. expr.', 'Single-att. diff. expr.',  
                   'Separated diff. expr.', 'Single-att. SSEM-lasso', 
                   'Sequential multi-att. NF')
method.colors <- c('black', "darkorange1", "deepskyblue", "purple" )[c(1,1,1,2,2,2,3,4)]
method.lines <-c(1,2,3,1,3,2,2,4)
method.pch <- c(16,10,1,17,2,14,4,15)

# Table for p(first) -- no exports currently specified.
prob.tables1 <- sapply(1:length(res.list), function(i){
                        topNProb(res.list[[i]], truth.list[[i]], method.labels, n=n.perturb*1)})
prob.tables1[2,] <- sapply(1:length(res.list), function(i){
                        topNProb(res.list[[i]][2], truth.list[[i]][2], method.labels[2], n=n.perturb*2)})
prob.tables1[6,] <- sapply(1:length(res.list), function(i){
                        topNProb(res.list[[i]][6], truth.list[[i]][6], method.labels[6], n=n.perturb*2)})
prob.tab.all1<- data.frame(auc.df$r.in, auc.df$r.out, t(prob.tables1))
names(prob.tab.all1) <- c('rin','rout', method.labels)

# Table for AUC
auc.tables1 <- sapply(1:length(res.list), function(i){
  aucList(res.list[[i]], truth.list[[i]], method.labels, r.in=auc.df$r.in, r.out=auc.df$r.out)})
auc.tables1[2,] <- sapply(1:length(res.list), function(i){
  aucList(res.list[[i]][2], truth.list[[i]][2], method.labels[2], r.in=auc.df$r.in, r.out=auc.df$r.out)})
auc.tables1[6,] <- sapply(1:length(res.list), function(i){
  aucList(res.list[[i]][6], truth.list[[i]][6], method.labels[6], r.in=auc.df$r.in, r.out=auc.df$r.out)})
auc.tab.all1<- data.frame(auc.df$r.in, auc.df$r.out, t(auc.tables1))
names(auc.tab.all1) <- c('rin','rout', method.labels)


first <- TRUE
for(i in 1:nrow(auc.df)){
  perflist <- getPerformance(res.list[[i]], truth.list[[i]])
  l <- 1
  for(l in 1:length(perflist)){
    tmp <- TidyPerformance(perflist[[l]], label=method.labels[l])
    tmp$r.in <- auc.df$r.in[i]
    tmp$r.out <- auc.df$r.out[i]
    tmp$pch <- method.pch[l]
    tmp$lty <- method.lines[l]
    
    tmp0 <- data.table(threshold=0, x=0, y=0, label=method.labels[l],
                       r.in=auc.df$r.in[i], 
                       r.out=auc.df$r.out[i],
                       pch=method.pch[l],
                       lty=method.lines[l])
    if(first) {
      longres <- rbindlist(list(tmp, tmp0))
    } else {
      longres <- rbindlist(list(longres, tmp, tmp0))
    }
    first <- FALSE
  }
}


names(method.pch) <- method.labels
names(method.colors) <- method.labels
names(method.lines) <- method.labels

label_parseall <- function(variable, value) {
  if(variable=='r.in'){
    plyr::llply(value, function(x) {
      parse(text = paste('rho["in"]', x, sep = "=="))})
  } else if(variable=='r.out'){
    plyr::llply(value, function(x) {
      parse(text = paste('rho["out"]', x, sep = "=="))})
  }
}

# Subsets for plotting/
inc.list <- list( c(8,1,2,3,4,6,5,7),
                  c(1,3,4,5,7),
                  c(1,2,3),
                  c(8,1,4))
inc.labels <- c('all', 'method', 'lrt', 'cond')

longres$lty <- as.factor(longres$lty)
longres$pch <- as.factor(longres$pch)
longres$label <- factor(longres$label, levels=method.labels[inc.list[[1]]], ordered=TRUE)

for(inc in 1:length(inc.labels)){
  l.inc <- inc.list[[inc]]
  tmp <- longres[longres$label %in% method.labels[l.inc],]
  tmp$label <- factor(tmp$label, levels=method.labels[l.inc], ordered=TRUE)
  plt <- ggplot(data = tmp, aes_string(x = "x", y = "y", color = "label", 
                                           linetype="label")) +
      geom_line(size=.7) + 
      geom_abline(intercept = 0, slope = 1,  color = "gray", linetype = 1) +
      scale_color_manual(values=method.colors[l.inc]) +
      scale_linetype_manual(values=method.lines[l.inc]) +
      ylab("Proportion of perturbations found") + 
      xlab("Proportion of sites considered") +
      theme_bw() + 
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            strip.text.x = element_text(size = 12),
            strip.text.y = element_text(size = 12),
            legend.text = element_text(size=12), 
            legend.key.width = unit(1, "cm")) +
      facet_grid(sprintf("%s ~ %s", 'r.in', 'r.out'), labeller=label_parseall) 
  plt 
   ggsave(filename=sprintf('%s/ggroc_%s.pdf', post_dir, inc.labels[inc]),
         width=9.5, height=6.8)
  rm(plt)
  
  c.rep <- paste(rep('c', length(l.inc)), collapse='')
  auc.tab <- auc.tab.all1[,c(1:2,(l.inc+2))]
  #names(auc.tab) <- c('r.in','r.out', method.labels[l.inc])
  print(xtable(auc.tab, align=sprintf('ccc|%s', c.rep), 
               digits=c(rep(1,3), rep(2,length(l.inc)))),
        file=sprintf('%s/aucs_%s_%s_%s.tex', post_dir, inc.labels[inc], fac, 
                     ideal, table.placement='H', include.rownames=FALSE, 
                     include.colnames=TRUE ))

  print.tab1 <- prob.tab.all1[,c(1:2,(l.inc+2))]
  print(xtable(print.tab1, align=sprintf('ccc|%s', c.rep), 
               digits=c(rep(1,3), rep(2,length(l.inc)))),
        file=sprintf('%s/prob1_%s_%s_%s.tex', post_dir, inc.labels[inc], fac, 
                     ideal, table.placement='H', include.rownames=FALSE, 
                     include.colnames=TRUE ))
}