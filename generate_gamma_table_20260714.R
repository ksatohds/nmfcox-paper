# 【2026-07-14 査読 r1 Major#2】3解析例の PH 変数 γ を3手法で比較（Fig8 fig:gamma の元データ）。
#   mgus2 を全 complete-case (N=1338) に変更し、cox.zph 準拠で z=sex+mspike, a=age+hgb+creat（旧: 先頭700, z=sex+hgb+creat, a=age+mspike）。
#  ① z-only Cox ② full-PH Cox ③ NMF-COX cross-fit（nmf.cox3.cf）。ストーリー: γ は保存（③≈②）。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))

DS <- list(
  veteran = list(surv="Surv(time,status)",  z="trt+age+diagtime+prior", a="karno+celltype",       SM=1000),
  gbsg    = list(surv="Surv(rfstime,status)",z="hormon+nodes+size",      a="age+meno+grade+pgr+er",SM=3000),
  mgus2   = list(surv="Surv(futime,death)",  z="sex+mspike",             a="age+hgb+creat",        SM=3000))
getdata <- function(nm){
  if(nm=="veteran") d<-survival::veteran
  if(nm=="gbsg")    d<-survival::gbsg
  if(nm=="mgus2"){ d<-survival::mgus2; d<-d[complete.cases(d[,c("futime","death","sex","age","hgb","creat","mspike")]),] }  # 全 complete-case
  d
}
ci <- function(est, se) c(est, est-1.96*se, est+1.96*se)
fmt <- function(v) sprintf("%+.3f (%+.3f, %+.3f)", v[1], v[2], v[3])

rows <- list()
for(nm in names(DS)){ p<-DS[[nm]]; d<-getdata(nm)
  mr <- model.frame(as.formula(paste0(p$surv,"~1")),d)[[1]]; d$y<-Surv(mr[,1], mr[,2])
  c0 <- coxph(as.formula(paste0(p$surv,"~",p$z)), data=d); s0<-sqrt(diag(vcov(c0)))
  c1 <- coxph(as.formula(paste0(p$surv,"~",p$z,"+",p$a)), data=d); s1<-sqrt(diag(vcov(c1)))
  cf <- nmf.cox.cf(as.formula(paste0("y~",p$z)), data=d, A=as.formula(paste0("~",p$a)),
                    rank=2, X.L2.smooth=p$SM, nfolds=5, seed=1)
  zc <- names(coef(c0))
  cat(sprintf("\n===== %s (n=%d) =====\n", nm, nrow(d)))
  for(v in zc){
    g0<-ci(coef(c0)[v], s0[v]); g1<-ci(coef(c1)[v], s1[v]); g3<-ci(cf$gamma[v], cf$se[v])
    cat(sprintf("%-10s | %-26s | %-26s | %-26s\n", v, fmt(g0), fmt(g1), fmt(g3)))
    rows[[length(rows)+1]] <- data.frame(data=nm, var=v,
       cox_z=fmt(g0), cox_full=fmt(g1), nmfcox=fmt(g3),
       e0=g0[1], e1=g1[1], e3=g3[1], stringsAsFactors=FALSE)
  }
}
tab <- do.call(rbind, rows)
saveRDS(tab, "verification/gamma_table_20260714.rds")
write.csv(tab[,1:5], "verification/gamma_table_20260714.csv", row.names=FALSE)
cat(sprintf("\n② full-PH と ③ NMF-COX の γ̂ 差: 平均|Δ|=%.3f 最大|Δ|=%.3f\n",
    mean(abs(tab$e1-tab$e3)), max(abs(tab$e1-tab$e3))))
cat("saved verification/gamma_table_20260714.{rds,csv}\n=== done ===\n")
