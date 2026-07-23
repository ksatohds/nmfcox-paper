# Monte Carlo: 二段階クロスフィット γ̂（nmf.cox3.cf）の定理1検証（オラクル同値・カバレッジ・安定性）。
#  設定ごとに R 反復を mclapply 並列。各反復で cross-fit / 層別オラクル / プールドオラクル / 全標本 / naive Cox の γ̂。
#  実行(サーバー): OPENBLAS_NUM_THREADS=1 Rscript verification/mc_cf_oracle_20260712.R <settings> <R> <cores> <outtag>
#    <settings>: 実行する設定番号のカンマ区切り（例 1,2）。<R>:反復数。<cores>:mc.cores。<outtag>:出力接尾辞。
#  進捗は verification/mc_progress_<outtag>.txt に追記（pollable）。結果は verification/mc_cf_oracle_<outtag>.rds。
suppressMessages({library(survival); library(parallel)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))

args <- commandArgs(trailingOnly=TRUE)
sets  <- if(length(args)>=1) as.integer(strsplit(args[1],",")[[1]]) else 1L
Rrep  <- if(length(args)>=2) as.integer(args[2]) else 50L
cores <- if(length(args)>=3) as.integer(args[3]) else 4L
outtag<- if(length(args)>=4) args[4] else "local"
progf <- sprintf("verification/mc_progress_%s.txt", outtag)
cat(sprintf("MC start: settings=%s R=%d cores=%d tag=%s\n", paste(sets,collapse=","), Rrep, cores, outtag), file=progf)

## 設定: (n, beta_A, rank)。beta_A=非PH共変量 a の時間変動効果強度、rank=cf の作業ランク。
SETTINGS <- list(
  `1`=list(n=300, betaA=0.4, rank=1, tag="n300_bA0.4_r1"),
  `2`=list(n=600, betaA=0.4, rank=1, tag="n600_bA0.4_r1"),
  `3`=list(n=300, betaA=0.8, rank=1, tag="n300_bA0.8_r1"),
  `4`=list(n=300, betaA=0.4, rank=2, tag="n300_bA0.4_r2"))
GAMMA0 <- 0.5; TAU <- 3; SM <- 300; NFOLD <- 5

sim_data <- function(n, betaA, seed){ set.seed(seed)
  z<-rnorm(n); a<-rnorm(n); K<-2000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K
  H<-exp(outer(z*GAMMA0,rep(1,K+1))+outer(a*betaA,grid/TAU)); cumH<-t(apply(H[,-1,drop=FALSE]*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx+1],Inf); Cc<-pmin(rexp(n,rate=0.6),TAU)
  data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z,a=a) }

## 真オフセット rho=a*betaA*(t/TAU) を既知としたオラクル: strat=strata(fold), pool=全標本。
oracle <- function(d, betaA, fold=NULL){
  time<-d$time; status<-d$status; et<-sort(unique(time[status==1])); K<-length(et); s_prev<-c(0,et[-K])
  rid<-rj<-rstart<-rstop<-rev<-integer(0)
  for(j in seq_len(K)){ risk<-which(time>=et[j]); m<-length(risk)
    rid<-c(rid,risk); rj<-c(rj,rep(j,m)); rstart<-c(rstart,rep(s_prev[j],m))
    rstop<-c(rstop,rep(et[j],m)); rev<-c(rev,as.integer(status[risk]==1 & time[risk]==et[j])) }
  off<-d$a[rid]*betaA*(et[rj]/TAU)
  dd<-data.frame(start=as.numeric(rstart),stop=as.numeric(rstop),ev=rev,id=rid,off=off,z=d$z[rid])
  if(is.null(fold)) coef(coxph(Surv(start,stop,ev)~z+offset(off)+cluster(id), data=dd))["z"]
  else { dd$fold<-fold[rid]; coef(coxph(Surv(start,stop,ev)~z+offset(off)+strata(fold)+cluster(id), data=dd))["z"] }
}

one_rep <- function(r, S){
  d <- sim_data(S$n, S$betaA, seed=r); d$y <- Surv(d$time, d$status)
  out <- c(cf=NA, se=NA, orc_fold=NA, orc_pool=NA, full=NA, naive=NA, nconv=NA)
  cf <- tryCatch(nmf.cox.cf(y~z, data=d, A=~a, rank=S$rank, X.L2.smooth=SM, nfolds=NFOLD, seed=r),
                 error=function(e) NULL)
  if(!is.null(cf)){ out["cf"]<-cf$gamma["z"]; out["se"]<-cf$se["z"]; out["nconv"]<-mean(cf$nuisance.converged)
    out["orc_fold"]<-oracle(d, S$betaA, cf$fold); out["orc_pool"]<-oracle(d, S$betaA, NULL) }
  fs <- tryCatch(nmf.cox(y~z, data=d, A=~a, rank=S$rank, X.L2.smooth=SM, ties="breslow",
                          maxit=30, verbose=FALSE, inference=FALSE, seed=r), error=function(e) NULL)
  if(!is.null(fs)) out["full"]<-fs$gamma["z"]
  out["naive"]<-coef(coxph(Surv(time,status)~z, data=d))["z"]
  out
}

results <- list()
for(si in sets){ S <- SETTINGS[[as.character(si)]]
  t0 <- Sys.time()
  cat(sprintf("[%s] setting %d (%s) start R=%d\n", format(Sys.time(),"%H:%M:%S"), si, S$tag, Rrep), file=progf, append=TRUE)
  mm <- mclapply(seq_len(Rrep), function(r){
    v <- tryCatch(one_rep(r, S), error=function(e) rep(NA,7))
    if(r %% 25 == 0) cat(sprintf("[%s] setting %d: %d/%d done\n", format(Sys.time(),"%H:%M:%S"), si, r, Rrep), file=progf, append=TRUE)
    v }, mc.cores=cores)
  M <- do.call(rbind, mm); colnames(M) <- c("cf","se","orc_fold","orc_pool","full","naive","nconv")
  results[[S$tag]] <- list(setting=S, M=M)
  cov <- mean(abs(M[,"cf"]-GAMMA0) <= 1.96*M[,"se"], na.rm=TRUE)
  cat(sprintf("[%s] setting %d DONE (%.1f min). cf: mean=%.3f sd=%.3f bias=%+.3f | cov95=%.3f | |cf-orcfold| mean=%.4f | fails=%d\n",
      format(Sys.time(),"%H:%M:%S"), si, as.numeric(difftime(Sys.time(),t0,units="mins")),
      mean(M[,"cf"],na.rm=TRUE), sd(M[,"cf"],na.rm=TRUE), mean(M[,"cf"],na.rm=TRUE)-GAMMA0,
      cov, mean(abs(M[,"cf"]-M[,"orc_fold"]),na.rm=TRUE), sum(is.na(M[,"cf"]))), file=progf, append=TRUE)
}
saveRDS(results, sprintf("verification/mc_cf_oracle_%s.rds", outtag))
cat(sprintf("[%s] ALL DONE tag=%s\n", format(Sys.time(),"%H:%M:%S"), outtag), file=progf, append=TRUE)
