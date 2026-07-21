# mgus2 全 complete-case (N=1338), A={age,hgb,creat}, z={sex,mspike} の λ_X を CVPL で選択（査読 r1 Major#2）。
#   nmf.cox3.cv の正しい引数: smooth.grid（=λ_X グリッド）, rank.grid, criterion。返り値 $cvpl/$se/$best.smooth/$best.smooth.1se。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))
args<-commandArgs(trailingOnly=TRUE); MC<-if(length(args)>=1) as.integer(args[1]) else 1L
cc<-c("futime","death","sex","age","hgb","creat","mspike"); d<-survival::mgus2; d<-d[complete.cases(d[,cc]),]; d$y<-Surv(d$futime,d$death)
cat(sprintf("N=%d events=%d\n", nrow(d), sum(d$death)))
gr <- c(100,300,1000,3000,10000)
cv <- nmf.cox.cv(y~sex+mspike, data=d, A=~age+hgb+creat, rank=2, X.L2.smooth=gr,
                  criterion="cvpl", nfolds=5, ties="breslow", seed=1, x.update="newton", mc.cores=MC, verbose=TRUE)
cat("\n--- extracted ---\n")
cat("CVPL:\n"); print(round(cv$cvpl,3))
cat("SE:\n");   print(round(cv$se,3))
cat(sprintf("best.smooth (CVPL argmax) = %g\n", cv$X.L2.smooth.best))
cat(sprintf("best.smooth.1se (most parsimonious within 1-SE) = %g\n", cv$X.L2.smooth.best.1se))
