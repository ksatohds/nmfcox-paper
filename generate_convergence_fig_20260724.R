# 【20260724】フォント2倍版: Fig.1 収束図の全文字（軸ラベル・目盛・タイトル）を 2x に拡大。
#   r2 からの変更は draw() の cex.lab/cex.axis/cex.main=2 と、はみ出し防止の余白(mar/mgp)・出力高さのみ。
#   モデル・データ・数値は r2 と同一（mgus2 全1338・cox.zph 準拠 split）。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))

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
  # フォント2倍 → 目盛数字・軸タイトルのはみ出し防止に余白と mgp を拡大
  par(mfrow=c(1,3), mar=c(6.2,7.6,4.4,1.4), mgp=c(4.4,1.1,0))
  for(nm in names(res)){ p<-res[[nm]]
    plot(seq_along(p), p, type="b", pch=19, lwd=2.4, cex=1.3, col="black",
         xlab="iteration", ylab="penalised partial log-likelihood", main=sprintf("(%s)", nm),
         cex.lab=2, cex.axis=2, cex.main=2)
    grid(nx=NA, ny=NULL, col="gray85", lty=1) }
}
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
# 高さを 680→900 に拡大（2x フォントの縦方向の余裕を確保）
png("v3_fig/convergence.png", width=2050, height=900, res=150, type=ptype); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/convergence.pdf", width=2050/150, height=900/150); draw(); dev.off()
cat("=== monotone fraction (non-decreasing penalised loglik) ===\n")
for(nm in names(res)) cat(sprintf("  %-8s (B) mono=%.2f  iters=%d  final=%.2f\n", nm, mono(res[[nm]]), length(res[[nm]]), tail(res[[nm]],1)))
cat("saved v3_fig/convergence.{png,pdf} (fonts 2x)\n=== done ===\n")
