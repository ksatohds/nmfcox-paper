# 【20260724 監査】必須-6: ランク過小指定が未検証なのに、実データ3件はすべて Q=2 < R=4/5/3 の過小指定領域で運用されている。
#   真ランク2 の DGP（a1 は u=t/tau、a2 は 1-u に載る）を、作業ランク Q=1（過小）と Q=2（正しい）で同一データに当て、
#   gamma の bias / 95%cov / 層別oracle 差を比較する。過小指定の代償を直接測る。
# tab:mcsim の「over-rank」設定を R=2 に再設計（査読 Major#1: 旧設定は R=1 で作業rank2＝Q>R 非同定）。
#   R=2 共変量 (a1,a2)、真 rank 1（両者が1つの ramp 時間関数に載る: β_l(t)=θ_l·(t/τ), θ=(0.4,0.4)）、作業 rank Q=2（=R, Q≤R の妥当な過大指定）。
#   z⊥a, γ0=0.5, n=300, 500反復。cross-fit の bias/95%cov/oracle差/fails と 層別oracle mean±sd を出力。
#   実行(サーバー): OPENBLAS_NUM_THREADS=1 Rscript verification/mc_overrank_R2_20260714.R <R> <cores> <tag>
suppressMessages({library(survival); library(parallel)})
# Working directory = repository root. If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))
args<-commandArgs(trailingOnly=TRUE)
Rrep  <- if(length(args)>=1) as.integer(args[1]) else 20L
cores <- if(length(args)>=2) as.integer(args[2]) else 1L
outtag<- if(length(args)>=3) args[3] else "local"
OFFS  <- if(length(args)>=4) as.integer(args[4]) else 0L   # seed オフセット（2台分割用）
GAMMA0<-0.5; TAU<-3; SM<-300; NFOLD<-5; N<-300; TH<-c(0.5,0.5)   # θ=(a1,a2 の負荷); a1->u, a2->(1-u) の異なる時間関数＝真 rank 2
progf <- sprintf("verification/mc_underrank_progress_%s.txt", outtag)
cat(sprintf("MC-underrank(true rank2) start: R=%d cores=%d tag=%s\n",Rrep,cores,outtag),file=progf)

sim_data <- function(n, seed){ set.seed(seed)
  z<-rnorm(n); a1<-rnorm(n); a2<-rnorm(n); K<-2000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K
  u<-grid/TAU
  lp<-outer(z*GAMMA0,rep(1,K+1)) + outer(a1*TH[1],u) + outer(a2*TH[2],1-u)   # 異なる時間関数＝真 rank 2
  H<-exp(lp); cumH<-t(apply(H[,-1,drop=FALSE]*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx+1],Inf); Cc<-pmin(rexp(n,rate=0.6),TAU)
  data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z,a1=a1,a2=a2) }

oracle <- function(d, fold=NULL){
  time<-d$time; status<-d$status; et<-sort(unique(time[status==1])); K<-length(et); s_prev<-c(0,et[-K])
  rid<-rj<-rstart<-rstop<-rev<-integer(0)
  for(j in seq_len(K)){ risk<-which(time>=et[j]); m<-length(risk)
    rid<-c(rid,risk); rj<-c(rj,rep(j,m)); rstart<-c(rstart,rep(s_prev[j],m)); rstop<-c(rstop,rep(et[j],m))
    rev<-c(rev,as.integer(status[risk]==1 & time[risk]==et[j])) }
  uu<-et[rj]/TAU
  off<-d$a1[rid]*TH[1]*uu + d$a2[rid]*TH[2]*(1-uu)                        # 真オフセット(rank2)
  dd<-data.frame(start=as.numeric(rstart),stop=as.numeric(rstop),ev=rev,id=rid,off=off,z=d$z[rid])
  if(is.null(fold)) coef(coxph(Surv(start,stop,ev)~z+offset(off)+cluster(id), data=dd))["z"]
  else { dd$fold<-fold[rid]; coef(coxph(Surv(start,stop,ev)~z+offset(off)+strata(fold)+cluster(id), data=dd))["z"] } }

one_rep <- function(r){
  d <- sim_data(N, seed=r); d$y <- Surv(d$time, d$status)
  out <- c(cf1=NA, se1=NA, cf2=NA, se2=NA, orc_fold=NA, nconv1=NA, nconv2=NA)
  c1 <- tryCatch(nmf.cox.cf(y~z, data=d, A=~a1+a2, rank=1, X.L2.smooth=SM, nfolds=NFOLD, seed=r), error=function(e) NULL)  # 過小 Q=1
  c2 <- tryCatch(nmf.cox.cf(y~z, data=d, A=~a1+a2, rank=2, X.L2.smooth=SM, nfolds=NFOLD, seed=r), error=function(e) NULL)  # 正しい Q=2
  if(!is.null(c1)){ out["cf1"]<-c1$gamma["z"]; out["se1"]<-c1$se["z"]; out["nconv1"]<-mean(c1$nuisance.converged); out["orc_fold"]<-oracle(d, c1$fold) }
  if(!is.null(c2)){ out["cf2"]<-c2$gamma["z"]; out["se2"]<-c2$se["z"]; out["nconv2"]<-mean(c2$nuisance.converged) }
  out }

t0<-Sys.time()
mm<-mclapply(seq_len(Rrep), function(i){ r<-OFFS+i; v<-tryCatch(one_rep(r),error=function(e)c(cf1=NA,se1=NA,cf2=NA,se2=NA,orc_fold=NA,nconv1=NA,nconv2=NA))
  if(r%%50==0) cat(sprintf("[%s] %d/%d\n",format(Sys.time(),"%H:%M:%S"),r,Rrep),file=progf,append=TRUE); v }, mc.cores=cores)
M<-do.call(rbind,mm); saveRDS(M, sprintf("verification/mc_underrank_%s.rds", outtag))
ok<-is.finite(M[,"cf"])&is.finite(M[,"se"]); cf<-M[ok,"cf"]; se<-M[ok,"se"]; of<-M[ok,"orc_fold"]
cat("=== done ===\n")
