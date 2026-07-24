# 【20260724 監査】E: cox.zph による z/a 分割の感度（transform × 閾値）／D: リスク集合内相関 |Cor_{π_t}(z,a)| の β_A 依存。
#   E → 必須-5「分割が cox.zph の既定 transform 1つで反転しうる」の検証。
#   D → 必須-2「Assumption 3(直交性)は非PH強度とともに破れるはずだが β_A=3 での実測が無い」の検証。
suppressMessages(library(survival))
# Working directory = repository root. If needed: setwd("/path/to/nmfcox-paper")

## ---------- E: z/a split sensitivity ----------
cat("=========== E: cox.zph z/a split sensitivity (p-values) ===========\n")
DS <- list()
{ d<-survival::veteran; d$y<-Surv(d$time,d$status); DS$veteran<-list(d=d, f=y~trt+age+diagtime+prior+karno+celltype) }
{ d<-survival::gbsg; d$y<-Surv(d$rfstime,d$status); DS$gbsg<-list(d=d, f=y~hormon+nodes+size+age+meno+grade+pgr+er) }
{ d<-survival::mgus2; d<-d[complete.cases(d[,c("futime","death","sex","age","hgb","creat","mspike")]),]
  d$y<-Surv(d$futime,d$death); DS$mgus2<-list(d=d, f=y~sex+age+hgb+creat+mspike) }

for(nm in names(DS)){ P<-DS[[nm]]; fit<-coxph(P$f, data=P$d)
  cat(sprintf("\n--- %s ---\n", nm))
  tabs <- lapply(c("km","rank","identity","log"), function(tr)
    tryCatch(cox.zph(fit, transform=tr)$table[,"p"], error=function(e) NULL))
  names(tabs) <- c("km","rank","identity","log")
  vars <- setdiff(rownames(cox.zph(fit)$table), "GLOBAL")
  cat(sprintf("%-12s %8s %8s %8s %8s | split@.05 (km->identity)\n","covariate","km","rank","identity","log"))
  for(v in vars){
    p <- sapply(tabs, function(z) if(is.null(z)) NA else unname(z[v]))
    flip <- if(!is.na(p["km"])&&!is.na(p["identity"]) && ((p["km"]<.05) != (p["identity"]<.05))) "  <== FLIPS" else ""
    cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f | %s%s\n", v, p["km"],p["rank"],p["identity"],p["log"],
        ifelse(p["km"]<.05,"a","z"), flip)) }
}

## ---------- D: within-risk-set correlation vs non-PH strength ----------
cat("\n\n=========== D: mean |Cor_{pi_t}(z,a)| vs beta_A (z indep a at baseline) ===========\n")
GAMMA0<-0.5; TAU<-3
sim_data <- function(n, betaA, seed){ set.seed(seed)
  z<-rnorm(n); a<-rnorm(n); K<-2000; grid<-seq(0,TAU,length.out=K+1); ds<-TAU/K
  H<-exp(outer(z*GAMMA0,rep(1,K+1))+outer(a*betaA,grid/TAU))
  cumH<-t(apply(H[,-1,drop=FALSE]*ds,1,cumsum)); E<-rexp(n)
  idx<-max.col(-(cumH<E),ties.method="first"); crossed<-cumH[cbind(1:n,idx)]>=E
  Tt<-ifelse(crossed,grid[idx+1],Inf); Cc<-pmin(rexp(n,rate=0.6),TAU)
  data.frame(time=pmin(Tt,Cc),status=as.integer(Tt<=Cc),z=z,a=a) }
wcor <- function(z,a,w){ w<-w/sum(w); mz<-sum(w*z); ma<-sum(w*a)
  vz<-sum(w*(z-mz)^2); va<-sum(w*(a-ma)^2); if(vz<=0||va<=0) return(NA)
  sum(w*(z-mz)*(a-ma))/sqrt(vz*va) }
cat(sprintf("%6s %6s %10s %10s %10s\n","n","beta_A","mean|Cor|","max|Cor|","events"))
for(n in c(300,600)) for(bA in c(0.4,0.8,1,2,3,4)){
  vals<-c(); ev<-c()
  for(r in 1:40){ d<-sim_data(n,bA,seed=r); et<-sort(unique(d$time[d$status==1])); ev<-c(ev,length(et))
    cs<-sapply(et, function(t){ rk<-which(d$time>=t)
      if(length(rk)<5) return(NA)
      w<-exp(d$z[rk]*GAMMA0 + d$a[rk]*bA*(t/TAU)); wcor(d$z[rk],d$a[rk],w) })
    vals<-c(vals, mean(abs(cs),na.rm=TRUE)) }
  cat(sprintf("%6d %6.1f %10.4f %10.4f %10.1f\n", n, bA, mean(vals,na.rm=TRUE), max(vals,na.rm=TRUE), mean(ev))) }
cat("\n=== done ===\n")
