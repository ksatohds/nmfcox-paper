# (A) newton ソルバ堅牢性の全データ検証: nmf.cox3(newton) vs nmf.cox3(optim) vs nmf.cox2 を
#   sim + 論文6データ（veteran/gbsg/pharma/lung/pbc/mgus2）で比較。
#   収束・故障の有無、α̂ 一致、β(keyA)(t) 相関、penalized loglik（newton≥optim を期待）、速度。
# 実行(サーバー): OPENBLAS_NUM_THREADS=1 Rscript verification/verify_newton_alldata_20260711b.R
suppressMessages({library(survival); has.asaur <- requireNamespace("asaur", quietly=TRUE); if(has.asaur) library(asaur)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
source("nmf_cox2.R"); suppressMessages(library(nmfkc))
ALPHA0<-0.5; BETA_A<-0.4; TAU<-3

sim_data <- function(n=300, seed=1){ set.seed(seed)
  z<-rnorm(n); a<-rnorm(n); K<-2000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K
  H<-exp(outer(z*ALPHA0,rep(1,K+1))+outer(a*BETA_A,grid/TAU)); cumH<-t(apply(H[,-1,drop=FALSE]*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx],TAU); Cc<-pmin(rexp(n,rate=0.6),TAU)
  data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z,a=a) }

prep <- list()
{ d<-sim_data(300,1); d$y<-Surv(d$time,d$status); prep$sim<-list(d=d, zf="z", Af=~a, keyA="a", SM=300) }
{ d<-survival::veteran; d$y<-Surv(d$time,d$status); prep$veteran<-list(d=d, zf="trt+age+diagtime+prior", Af=~karno+celltype, keyA="karno", SM=1000) }
{ d<-survival::gbsg; d$y<-Surv(d$rfstime,d$status); prep$gbsg<-list(d=d, zf="hormon+age+meno", Af=~grade+nodes+size+pgr+er, keyA="pgr", SM=3000) }
if(has.asaur){ d<-subset(asaur::pharmacoSmoking,ttr>0); d$grp<-relevel(d$grp,ref="patchOnly"); d$y<-Surv(d$ttr,d$relapse)
  prep$pharma<-list(d=d, zf="grp", Af=~age+gender+yearsSmoking+levelSmoking+priorAttempts+longestNoSmoke, keyA="yearsSmoking", SM=1000) }
{ d<-survival::lung; d<-d[complete.cases(d[,c("time","status","sex","age","ph.ecog","ph.karno","wt.loss")]),]; d$y<-Surv(d$time,d$status==2)
  prep$lung<-list(d=d, zf="sex+age", Af=~ph.ecog+ph.karno+wt.loss, keyA="ph.karno", SM=3000) }
{ d<-survival::pbc; d<-d[!is.na(d$trt),]; d$logbili<-log(d$bili)
  d<-d[complete.cases(d[,c("time","status","trt","age","sex","logbili","albumin","edema","stage","protime")]),]; d$y<-Surv(d$time,d$status==2)
  prep$pbc<-list(d=d, zf="trt+age+sex", Af=~logbili+albumin+edema+stage+protime, keyA="protime", SM=1000) }
{ d<-survival::mgus2; d<-d[complete.cases(d[,c("futime","death","sex","age","hgb","creat","mspike")]),]
  d$y<-Surv(d$futime,d$death); prep$mgus2<-list(d=d, zf="sex+mspike", Af=~age+hgb+creat, keyA="age", SM=3000) }  # 全1338, cox.zph 準拠 split（査読 r1 Major#2）

cat(sprintf("%-9s %6s %5s | %-22s | %-18s | %-14s | %-12s\n","data","n","K","alpha(interest): cox2/opt/new","beta cor new:cox2/opt","pen.loglik c2/o/n","time o/n (conv)"))
for(nm in names(prep)){ p<-prep[[nm]]; fz<-as.formula(paste0("y~",p$zf)); ik<-strsplit(p$zf,"[+]")[[1]][1]
  ff<-function(upd,use3) if(use3) nmf.cox(fz,data=p$d,A=p$Af,rank=2,X.L2.smooth=p$SM,ties="breslow",maxit=30,verbose=FALSE,inference=FALSE,X.init="kmeans++",seed=1,x.update=upd) else nmf.cox2(fz,data=p$d,A=p$Af,rank=2,X.smooth=p$SM,ties="breslow",maxit=30,verbose=FALSE,inference=FALSE)
  r2<-tryCatch(system.time(f2<<-ff("optim",FALSE)),error=function(e){cat(nm,"cox2 ERR:",conditionMessage(e),"\n");NA})
  ro<-tryCatch(system.time(fo<<-ff("optim",TRUE)),error=function(e){cat(nm,"optim ERR:",conditionMessage(e),"\n");NA})
  rn<-tryCatch(system.time(fn<<-ff("newton",TRUE)),error=function(e){cat(nm,"NEWTON ERR:",conditionMessage(e),"\n");NA})
  if(any(is.na(c(r2[1],ro[1],rn[1])))) next
  et<-f2$event.times; al<-function(f) approx(f$event.times,f$beta.t[,p$keyA],xout=et,rule=2)$y
  b2<-f2$beta.t[,p$keyA]; bo<-al(fo); bn<-al(fn)
  cat(sprintf("%-9s %6d %5d | %+.3f/%+.3f/%+.3f | %.3f / %.3f | %.1f/%.1f/%.1f | %.1f/%.1fs (%s/%s)\n",
    nm, nrow(p$d), length(et), unname(f2$alpha[ik]),unname(fo$gamma[ik]),unname(fn$gamma[ik]),  # nmfkc の nmf.cox は $gamma のみ返す（$alpha 別名廃止）; f2=nmf.cox2 は $alpha 維持
    cor(b2,bn),cor(bo,bn), f2$best.obj,fo$best.obj,fn$best.obj, ro[3],rn[3], fo$converged,fn$converged)) }
cat("\n=== done ===\n")
