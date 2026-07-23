# 【20260724】偽イベント修正後の再実行結果を集計（tab:mcsim baseline 行 と fig:confound/§3.1 の数値）。
#   入力: mc_cf_oracle_cf724.rds（設定1-3=baseline/larger N/strong nuisance）, mc_conf_conf724.rds（rank=1 に修正）
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
GAMMA0 <- 0.5

cat("=================== tab:mcsim (cross-fit oracle, fixed DGP) ===================\n")
res <- readRDS("verification/mc_cf_oracle_cf724.rds")
summ <- function(x) c(mean=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE))
cat(sprintf("%-18s %5s %4s %5s | %8s %7s %9s | %-14s | %s\n",
    "setting","n","bA","rank","bias","cov95","|cf-orcF|","strat.oracle","fails/R"))
for(tag in names(res)){ S<-res[[tag]]$setting; M<-res[[tag]]$M
  cf<-summ(M[,"cf"]); orf<-summ(M[,"orc_fold"])
  bias<-cf["mean"]-GAMMA0
  cov <-mean(abs(M[,"cf"]-GAMMA0)<=1.96*M[,"se"],na.rm=TRUE)
  d_of<-mean(abs(M[,"cf"]-M[,"orc_fold"]),na.rm=TRUE)
  fails<-sum(is.na(M[,"cf"]))
  cat(sprintf("%-18s %5d %4.1f %5d | %+8.3f %7.3f %9.4f | %+.3f±%.3f | %d/%d\n",
      tag,S$n,S$betaA,S$rank,bias,cov,d_of,orf["mean"],orf["sd"],fails,nrow(M)))
}

cat("\n=================== fig:confound / §3.1 (mc_conf, rank=1) ===================\n")
cf2 <- readRDS("verification/mc_conf_conf724.rds")
cat("--- structure of one setting ---\n"); print(colnames(cf2[[1]]$M))
est <- colnames(cf2[[1]]$M)
cat(sprintf("\n%-22s %s\n","setting", paste(sprintf("%-16s",est), collapse="")))
for(tag in names(cf2)){ M<-cf2[[tag]]$M
  cat(sprintf("%-22s %s\n", tag,
      paste(sprintf("%-16s", sprintf("%+.4f", colMeans(M,na.rm=TRUE)-GAMMA0)), collapse=""))) }
cat("\n(above = bias of each column vs GAMMA0=0.5; coverage columns interpreted below)\n")
for(tag in names(cf2)){ M<-cf2[[tag]]$M
  cat(sprintf("\n[%s] n=%d rho=%.1f rank=%d\n", tag, cf2[[tag]]$setting$n, cf2[[tag]]$setting$rho, cf2[[tag]]$setting$rank))
  print(round(colMeans(M,na.rm=TRUE),4)) }
cat("\n=== done ===\n")
