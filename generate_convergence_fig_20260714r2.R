# 【r2 2026-07-14】mgus2 を全1338・cox.zph 準拠 split（z=sex+mspike, a=age+hgb+creat）に修正（他図表と整合）。
# (B) block-coordinate の単調収束図: veteran/gbsg/mgus2 の反復ごとの罰則付き部分対数尤度（(B) のみ、3パネル）。
#   objfunc.iter は最小化量（負の罰則付き部分対数尤度）→ 罰則付き loglik = -objfunc.iter を y 軸に（大きいほど良い）。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))   # (B): Θ 罰則付き block-coordinate

DS <- list(
  veteran = list(get=function() survival::veteran, tv="time", ev="status", zf="trt+age+diagtime+prior", af=~karno+celltype, SM=1000),
  gbsg    = list(get=function() survival::gbsg, tv="rfstime", ev="status", zf="hormon+nodes+size", af=~age+meno+grade+pgr+er, SM=3000),
  mgus2   = list(get=function(){ d<-survival::mgus2; d<-d[complete.cases(d[,c("futime","death","sex","age","hgb","creat","mspike")]),]; d},   # 全1338（査読 r1 Major#2）
                 tv="futime", ev="death", zf="sex+mspike", af=~age+hgb+creat, SM=3000))   # cox.zph 準拠 split
fit_iter <- function(P){ d<-P$get(); d$y <- Surv(d[[P$tv]], d[[P$ev]])
  f <- nmf.cox(as.formula(paste0("y~",P$zf)), data=d, A=P$af, rank=2, X.L2.smooth=P$SM, ties="breslow",
                maxit=30, verbose=FALSE, inference=FALSE, x.update="newton", seed=1)
  o <- f$objfunc.iter; -o[is.finite(o)] }               # 罰則付き部分対数尤度
res <- lapply(names(DS), function(nm) fit_iter(DS[[nm]])); names(res)<-names(DS)
mono <- function(v) if(length(v)<2) NA else round(mean(diff(v)>=-1e-8),2)
draw <- function(){
  par(mfrow=c(1,3), mar=c(4.2,4.8,3,1), mgp=c(2.6,0.8,0))
  for(nm in names(res)){ p<-res[[nm]]
    plot(seq_along(p), p, type="b", pch=19, lwd=2.2, cex=1.2, col="black",
         xlab="iteration", ylab="penalised partial log-likelihood", main=nm)
    grid(nx=NA, ny=NULL, col="gray85", lty=1) }
}
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
png("v3_fig/convergence.png", width=2050, height=680, res=150, type=ptype); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/convergence.pdf", width=2050/150, height=680/150); draw(); dev.off()
cat("=== monotone fraction (non-decreasing penalised loglik) ===\n")
for(nm in names(res)) cat(sprintf("  %-8s (B) mono=%.2f  iters=%d  final=%.2f\n", nm, mono(res[[nm]]), length(res[[nm]]), tail(res[[nm]],1)))
cat("saved v3_fig/convergence.{png,pdf}\n=== done ===\n")
