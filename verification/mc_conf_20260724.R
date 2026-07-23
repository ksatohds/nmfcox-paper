# Monte Carlo: 交絡＋非PH シミュレーション。z と a を相関(rho)させ、a を強い非PH(betaA(t)=cA*(t/TAU-0.5))に。
#  主張の実証: ① z-only Cox と ② full-PH Cox(a を PH で誤特定)は主効果 γ が偏り、③ NMF-COX(cross-fit)と
#  オラクル(真オフセット既知)が真値 γ0=0.5 を回復する。→「非PHを正しく推定すると主効果の偏りが補正される」。
#  実行(サーバー): OPENBLAS_NUM_THREADS=1 Rscript verification/mc_conf_20260712.R <settings> <R> <cores> <outtag>
suppressMessages({library(survival); library(parallel); library(splines)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))

args <- commandArgs(trailingOnly=TRUE)
sets  <- if(length(args)>=1) as.integer(strsplit(args[1],",")[[1]]) else 1L
Rrep  <- if(length(args)>=2) as.integer(args[2]) else 30L
cores <- if(length(args)>=3) as.integer(args[3]) else 1L
outtag<- if(length(args)>=4) args[4] else "local"
progf <- sprintf("verification/mc_conf_progress_%s.txt", outtag)
cat(sprintf("MC-conf start: settings=%s R=%d cores=%d tag=%s\n", paste(sets,collapse=","), Rrep, cores, outtag), file=progf)

GAMMA0<-0.5; TAU<-3; NFOLD<-5; H0<-0.4   # H0=一定ベースラインハザード（イベントを[0,TAU]に分散させ非PHを露出）
## 設定: (n, rho=corr(z,a), cA=非PH強度, rank, SM, censor.rate)
SETTINGS <- list(
  `1`=list(n=300, rho=0.5, cA=3, rank=1, SM=100, cens=0.25, tag="n300_rho.5_cA3_Q1"),  # 主デモ(交絡)
  `2`=list(n=600, rho=0.5, cA=3, rank=1, SM=100, cens=0.25, tag="n600_rho.5_cA3_Q1"),  # 大標本
  `3`=list(n=300, rho=0.0, cA=3, rank=1, SM=100, cens=0.25, tag="n300_rho.0_cA3_Q1"),  # z⊥a(②の非PH誤特定のみ)
  `4`=list(n=300, rho=0.7, cA=3, rank=1, SM=100, cens=0.25, tag="n300_rho.7_cA3_Q1"))  # 強交絡

betaA <- function(u, cA) cA*u                    # ramp 型 非PH（u=t/TAU∈[0,1]; 単調増加）

sim_conf <- function(n, rho, cA, cens, seed){ set.seed(seed)
  z<-rnorm(n); e<-rnorm(n); a<-rho*z + sqrt(1-rho^2)*e          # corr(z,a)=rho, 周辺 N(0,1)
  K<-3000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K; bg<-betaA(grid/TAU,cA)
  Hmid<-H0*exp(outer(z*GAMMA0,rep(1,K)) + outer(a,(bg[-1]+bg[-(K+1)])/2))  # 中点則; λ=H0 exp(zγ0+a βA(t))
  cumH<-t(apply(Hmid*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx+1],Inf); Cc<-pmin(rexp(n,rate=cens),TAU)     # 非交差=Inf→正しく打切り
  data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z,a=a) }

oracle_pooled <- function(d, cA){                # 真オフセット a*βA(t) 既知
  time<-d$time; status<-d$status; et<-sort(unique(time[status==1])); K<-length(et); s_prev<-c(0,et[-K])
  rid<-rj<-rstart<-rstop<-rev<-integer(0)
  for(j in seq_len(K)){ risk<-which(time>=et[j]); m<-length(risk)
    rid<-c(rid,risk); rj<-c(rj,rep(j,m)); rstart<-c(rstart,rep(s_prev[j],m))
    rstop<-c(rstop,rep(et[j],m)); rev<-c(rev,as.integer(status[risk]==1 & time[risk]==et[j])) }
  off<-d$a[rid]*betaA(et[rj]/TAU,cA)
  dd<-data.frame(start=as.numeric(rstart),stop=as.numeric(rstop),ev=rev,off=off,z=d$z[rid])
  coef(coxph(Surv(start,stop,ev)~z+offset(off), data=dd))["z"] }

one_rep <- function(r, S){
  d<-sim_conf(S$n,S$rho,S$cA,S$cens,seed=r); d$y<-Surv(d$time,d$status)
  out<-c(cox_z=NA,se0=NA,full=NA,se1=NA,nmf=NA,se3=NA,tvc=NA,se_tvc=NA,oracle=NA,evt=mean(d$status))
  c0<-tryCatch(coxph(Surv(time,status)~z,data=d),error=function(e)NULL)
  if(!is.null(c0)){ out["cox_z"]<-coef(c0)["z"]; out["se0"]<-sqrt(vcov(c0)["z","z"]) }
  c1<-tryCatch(coxph(Surv(time,status)~z+a,data=d),error=function(e)NULL)
  if(!is.null(c1)){ out["full"]<-coef(c1)["z"]; out["se1"]<-sqrt(vcov(c1)["z","z"]) }
  cf<-tryCatch(nmf.cox.cf(y~z,data=d,A=~a,rank=S$rank,X.L2.smooth=S$SM,nfolds=NFOLD,seed=r),error=function(e)NULL)
  if(!is.null(cf)){ out["nmf"]<-cf$gamma["z"]; out["se3"]<-cf$se["z"] }
  # tvc-Cox: z を PH, a を自然スプライン(df=3)の時間変動係数で（比較法）。
  ctvc<-tryCatch(coxph(Surv(time,status)~z+tt(a),data=d,tt=function(x,t,...) x*splines::ns(t,df=3)),error=function(e)NULL)
  if(!is.null(ctvc)){ out["tvc"]<-coef(ctvc)["z"]; out["se_tvc"]<-sqrt(vcov(ctvc)["z","z"]) }
  out["oracle"]<-tryCatch(oracle_pooled(d,S$cA),error=function(e)NA)
  out }

results<-list()
for(si in sets){ S<-SETTINGS[[as.character(si)]]; t0<-Sys.time()
  cat(sprintf("[%s] setting %d (%s) start R=%d\n", format(Sys.time(),"%H:%M:%S"), si, S$tag, Rrep), file=progf, append=TRUE)
  mm<-mclapply(seq_len(Rrep), function(r){ v<-tryCatch(one_rep(r,S),error=function(e)rep(NA,10))
    if(r%%25==0) cat(sprintf("[%s] setting %d: %d/%d\n",format(Sys.time(),"%H:%M:%S"),si,r,Rrep),file=progf,append=TRUE); v }, mc.cores=cores)
  M<-do.call(rbind,mm); colnames(M)<-c("cox_z","se0","full","se1","nmf","se3","tvc","se_tvc","oracle","evt")
  results[[S$tag]]<-list(setting=S,M=M)
  bias<-function(x)mean(x,na.rm=TRUE)-GAMMA0
  cvg<-function(est,se)mean(abs(est-GAMMA0)<=1.96*se,na.rm=TRUE)
  cat(sprintf("[%s] setting %d DONE (%.1fm). evt=%.2f | bias ①=%+.3f ②=%+.3f ③=%+.3f tvc=%+.3f orac=%+.3f | cov ②=%.2f ③=%.2f tvc=%.2f | fails③=%d\n",
    format(Sys.time(),"%H:%M:%S"), si, as.numeric(difftime(Sys.time(),t0,units="mins")), mean(M[,"evt"],na.rm=TRUE),
    bias(M[,"cox_z"]),bias(M[,"full"]),bias(M[,"nmf"]),bias(M[,"tvc"]),bias(M[,"oracle"]),
    cvg(M[,"full"],M[,"se1"]),cvg(M[,"nmf"],M[,"se3"]),cvg(M[,"tvc"],M[,"se_tvc"]), sum(is.na(M[,"nmf"]))), file=progf, append=TRUE)
}
saveRDS(results, sprintf("verification/mc_conf_%s.rds", outtag))
cat(sprintf("[%s] ALL DONE tag=%s\n", format(Sys.time(),"%H:%M:%S"), outtag), file=progf, append=TRUE)
