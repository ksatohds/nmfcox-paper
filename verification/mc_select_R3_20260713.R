# Monte Carlo: rank Q / 平滑化 λ_X の選択率。【20260713 再設計】非PH共変量数 R=3 に固定し、
#   評価ランク Q=1,2,3 を常に Q≤R=3 に収める（旧 mc_select_20260712.R は R=1 で Q=2,3 を評価＝Q>R 非同定で内部矛盾）。
#   2 DGP（いずれも z⊥a, γ0=0.5, R=3 共変量 a1,a2,a3）:
#     (A) 真 rank 1：3 共変量が 1 つの共有時間形 t/τ に載る  β_l(t)=θ_l·(t/τ), θ=(0.6,0.45,0.30)。
#     (B) 真 rank 2：3 共変量が 2 つの共有時間形 {t/τ, 1-t/τ} に載る
#          β1=0.7(t/τ), β2=0.7(1-t/τ), β3=0.5(t/τ)+0.3(1-t/τ)（3 本は 2 次元 span＝rank2）。
#   各反復で rank.grid=1:3, smooth.grid で規準最良の (rank,smooth) を記録→選択率。criterion で AIC/BIC/CVPL 切替。
#   実行: OPENBLAS_NUM_THREADS=1 Rscript verification/mc_select_R3_20260713.R <settings> <R> <cores> <tag> <criterion>
suppressMessages({library(survival); library(parallel)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))
args<-commandArgs(trailingOnly=TRUE)
sets  <- if(length(args)>=1) as.integer(strsplit(args[1],",")[[1]]) else 1L
Rrep  <- if(length(args)>=2) as.integer(args[2]) else 20L
cores <- if(length(args)>=3) as.integer(args[3]) else 1L
outtag<- if(length(args)>=4) args[4] else "local"
CRIT  <- if(length(args)>=5) args[5] else "aic"   # "aic" / "bic" / "cvpl"(1-SE 則で最簡素)
progf <- sprintf("verification/mc_selR3_progress_%s.txt", outtag)
cat(sprintf("MC-select(R3) start: settings=%s R=%d cores=%d tag=%s crit=%s\n",paste(sets,collapse=","),Rrep,cores,outtag,CRIT),file=progf)

GAMMA0<-0.5; TAU<-3; SMGRID<-c(30,100,300,1000); RGRID<-1:3
SETTINGS <- list(
  `1`=list(n=400, betas=list(function(u)0.6*u, function(u)0.45*u, function(u)0.30*u),
           true.rank=1, Af=~a1+a2+a3, tag="rank1_R3_n400"),
  `2`=list(n=500, betas=list(function(u)0.7*u, function(u)0.7*(1-u), function(u)0.5*u+0.3*(1-u)),
           true.rank=2, Af=~a1+a2+a3, tag="rank2_R3_n500"))

sim_multi <- function(n, betafuns, seed){ set.seed(seed)
  L<-length(betafuns); z<-rnorm(n); A<-matrix(rnorm(n*L),n,L)
  K<-2000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K; u<-grid/TAU
  lp<-outer(z*GAMMA0,rep(1,K+1)); for(l in 1:L) lp<-lp+outer(A[,l],betafuns[[l]](u))
  H<-exp(lp); cumH<-t(apply(H[,-1,drop=FALSE]*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx+1],Inf); Cc<-pmin(rexp(n,rate=0.6),TAU)
  d<-data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z)
  for(l in 1:L) d[[paste0("a",l)]]<-A[,l]; d$y<-Surv(d$time,d$status); d }

one_rep <- function(r, S){
  d<-sim_multi(S$n,S$betas,seed=r)
  cv<-tryCatch(nmf.cox.cv(y~z, data=d, A=S$Af, rank=RGRID, X.L2.smooth=SMGRID,
                           criterion=CRIT, nfolds=5, verbose=FALSE, seed=r), error=function(e)NULL)
  if(is.null(cv)) c(rank=NA, smooth=NA, rankmin=NA)
  else if(CRIT=="cvpl") c(rank=cv$rank.best.1se, smooth=cv$X.L2.smooth.best.1se, rankmin=cv$rank.best)  # 1-SE 則 と 素の最大
  else c(rank=cv$rank.best, smooth=cv$X.L2.smooth.best, rankmin=cv$rank.best) }

results<-list()
for(si in sets){ S<-SETTINGS[[as.character(si)]]; t0<-Sys.time()
  cat(sprintf("[%s] setting %d (%s) start R=%d crit=%s\n",format(Sys.time(),"%H:%M:%S"),si,S$tag,Rrep,CRIT),file=progf,append=TRUE)
  mm<-mclapply(seq_len(Rrep),function(r){ v<-tryCatch(one_rep(r,S),error=function(e)c(rank=NA,smooth=NA,rankmin=NA))
    if(r%%25==0) cat(sprintf("[%s] setting %d: %d/%d\n",format(Sys.time(),"%H:%M:%S"),si,r,Rrep),file=progf,append=TRUE); v },mc.cores=cores)
  M<-do.call(rbind,mm); colnames(M)<-c("rank","smooth","rankmin"); results[[S$tag]]<-list(setting=S,M=M)
  rt<-table(factor(M[,"rank"],levels=RGRID)); rtm<-table(factor(M[,"rankmin"],levels=RGRID))
  cat(sprintf("[%s] setting %d DONE (%.1fm). true.rank=%d | 1SE sel=%.3f dist=%s | MAX sel=%.3f dist=%s | median smooth=%g | fails=%d\n",
    format(Sys.time(),"%H:%M:%S"),si,as.numeric(difftime(Sys.time(),t0,units="mins")),S$true.rank,
    mean(M[,"rank"]==S$true.rank,na.rm=TRUE), paste(sprintf("Q%d:%.2f",RGRID,as.numeric(rt)/sum(rt)),collapse=" "),
    mean(M[,"rankmin"]==S$true.rank,na.rm=TRUE), paste(sprintf("Q%d:%.2f",RGRID,as.numeric(rtm)/sum(rtm)),collapse=" "),
    median(M[,"smooth"],na.rm=TRUE), sum(is.na(M[,"rank"]))),file=progf,append=TRUE) }
saveRDS(results, sprintf("verification/mc_selR3_%s.rds", outtag))
cat(sprintf("[%s] ALL DONE tag=%s crit=%s\n",format(Sys.time(),"%H:%M:%S"),outtag,CRIT),file=progf,append=TRUE)
