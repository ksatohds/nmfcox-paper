# Monte Carlo B: 査読対応の追加測定。γ推論(bias/empSD/meanSE/coverage/oracle gap)＋
#   β(t)推定精度(ISE/RMSE)＋β(t)のpointwise 95%被覆、を複数 n で（漸近同値のn依存を見る）。
#   真: λ_i(t)=exp(z·γ0 + a·βA·t/τ)（rank-1 非PH）、γ0=0.5、z⊥a。β_0(t)=βA·t/τ（raw a スケール）。
#   実行(サーバー): OPENBLAS_NUM_THREADS=1 Rscript verification/mc_beta_20260712.R <settings> <R> <cores> <tag>
suppressMessages({library(survival); library(parallel)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))
args<-commandArgs(trailingOnly=TRUE)
sets  <- if(length(args)>=1) as.integer(strsplit(args[1],",")[[1]]) else 1L
Rrep  <- if(length(args)>=2) as.integer(args[2]) else 30L
cores <- if(length(args)>=3) as.integer(args[3]) else 1L
outtag<- if(length(args)>=4) args[4] else "local"
progf <- sprintf("verification/mc_beta_progress_%s.txt", outtag)
cat(sprintf("MC-beta start: settings=%s R=%d cores=%d tag=%s\n", paste(sets,collapse=","),Rrep,cores,outtag), file=progf)

GAMMA0<-0.5; TAU<-3; SM<-300; NFOLD<-5
## 複数 n（漸近同値・oracle gap の n 依存）。betaA=0.4, rank=1 固定。
SETTINGS <- list(
  `1`=list(n=150,  betaA=0.4, rank=1, tag="n150"),
  `2`=list(n=300,  betaA=0.4, rank=1, tag="n300"),
  `3`=list(n=600,  betaA=0.4, rank=1, tag="n600"),
  `4`=list(n=1200, betaA=0.4, rank=1, tag="n1200"))

sim_data <- function(n, betaA, seed){ set.seed(seed)
  z<-rnorm(n); a<-rnorm(n); K<-2000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K
  H<-exp(outer(z*GAMMA0,rep(1,K+1))+outer(a*betaA,grid/TAU)); cumH<-t(apply(H[,-1,drop=FALSE]*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx+1],Inf); Cc<-pmin(rexp(n,rate=0.6),TAU)          # 非交差=Inf→正しく打切り
  data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z,a=a) }

oracle <- function(d, betaA, fold=NULL){
  time<-d$time; status<-d$status; et<-sort(unique(time[status==1])); K<-length(et); sp<-c(0,et[-K])
  rid<-rj<-rs<-rt<-rv<-integer(0)
  for(j in seq_len(K)){ rk<-which(time>=et[j]); m<-length(rk)
    rid<-c(rid,rk); rj<-c(rj,rep(j,m)); rs<-c(rs,rep(sp[j],m)); rt<-c(rt,rep(et[j],m)); rv<-c(rv,as.integer(status[rk]==1&time[rk]==et[j])) }
  off<-d$a[rid]*betaA*(et[rj]/TAU); dd<-data.frame(start=as.numeric(rs),stop=as.numeric(rt),ev=rv,id=rid,off=off,z=d$z[rid])
  if(is.null(fold)) coef(coxph(Surv(start,stop,ev)~z+offset(off)+cluster(id),data=dd))["z"]
  else { dd$fold<-fold[rid]; coef(coxph(Surv(start,stop,ev)~z+offset(off)+strata(fold)+cluster(id),data=dd))["z"] } }

one_rep <- function(r, S){
  d<-sim_data(S$n,S$betaA,seed=r); d$y<-Surv(d$time,d$status)
  out<-c(cf=NA,se=NA,orc_fold=NA,orc_pool=NA,naive=NA,beta_ise=NA,beta_rmse=NA,beta_cov=NA,nconv=NA)
  cf<-tryCatch(nmf.cox.cf(y~z,data=d,A=~a,rank=S$rank,X.L2.smooth=SM,nfolds=NFOLD,seed=r),error=function(e)NULL)
  if(!is.null(cf)){ out["cf"]<-cf$gamma["z"]; out["se"]<-cf$se["z"]; out["nconv"]<-mean(cf$nuisance.converged)
    out["orc_fold"]<-tryCatch(oracle(d,S$betaA,cf$fold),error=function(e)NA)
    out["orc_pool"]<-tryCatch(oracle(d,S$betaA,NULL),error=function(e)NA) }
  ## β(t) 推定精度＋pointwise被覆（全標本 fit, inference=TRUE）。raw a スケールで真 β0(t)=βA·t/τ と比較。
  fsI<-tryCatch(nmf.cox(y~z,data=d,A=~a,rank=S$rank,X.L2.smooth=SM,ties="breslow",
                         maxit=30,verbose=FALSE,inference=TRUE,seed=r),error=function(e)NULL)
  if(!is.null(fsI)){ et<-fsI$event.times; sc<-fsI$A.scale[1]
    bhat<-fsI$beta.t[,1]/sc; se<-fsI$se.beta.t[,1]/sc; b0<-S$betaA*(et/TAU)
    out["beta_ise"]<-mean((bhat-b0)^2); out["beta_rmse"]<-sqrt(mean((bhat-b0)^2))
    ok<-is.finite(se)&se>0; out["beta_cov"]<-if(any(ok)) mean(abs(bhat[ok]-b0[ok])<=1.96*se[ok]) else NA }
  out["naive"]<-coef(coxph(Surv(time,status)~z,data=d))["z"]
  out }

results<-list()
for(si in sets){ S<-SETTINGS[[as.character(si)]]; t0<-Sys.time()
  cat(sprintf("[%s] setting %d (%s) start R=%d\n",format(Sys.time(),"%H:%M:%S"),si,S$tag,Rrep),file=progf,append=TRUE)
  mm<-mclapply(seq_len(Rrep),function(r){ v<-tryCatch(one_rep(r,S),error=function(e)rep(NA,9))
    if(r%%25==0) cat(sprintf("[%s] setting %d: %d/%d\n",format(Sys.time(),"%H:%M:%S"),si,r,Rrep),file=progf,append=TRUE); v },mc.cores=cores)
  M<-do.call(rbind,mm); colnames(M)<-c("cf","se","orc_fold","orc_pool","naive","beta_ise","beta_rmse","beta_cov","nconv")
  results[[S$tag]]<-list(setting=S,M=M)
  cov<-mean(abs(M[,"cf"]-GAMMA0)<=1.96*M[,"se"],na.rm=TRUE)
  cat(sprintf("[%s] setting %d DONE (%.1fm). cf bias=%+.3f empSD=%.3f meanSE=%.3f cov95=%.3f | oracgap=%.4f | betaRMSE=%.3f betaCov=%.3f | fails=%d\n",
    format(Sys.time(),"%H:%M:%S"),si,as.numeric(difftime(Sys.time(),t0,units="mins")),
    mean(M[,"cf"],na.rm=TRUE)-GAMMA0, sd(M[,"cf"],na.rm=TRUE), mean(M[,"se"],na.rm=TRUE), cov,
    mean(abs(M[,"cf"]-M[,"orc_fold"]),na.rm=TRUE), mean(M[,"beta_rmse"],na.rm=TRUE), mean(M[,"beta_cov"],na.rm=TRUE),
    sum(is.na(M[,"cf"]))),file=progf,append=TRUE) }
saveRDS(results, sprintf("verification/mc_beta_%s.rds", outtag))
cat(sprintf("[%s] ALL DONE tag=%s\n",format(Sys.time(),"%H:%M:%S"),outtag),file=progf,append=TRUE)
