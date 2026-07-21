# 【20260714】Fig4 fig:v3basis（上=共有非負基底X, 下=非PH共変量の β(t)）を再描画。フォント1.5倍（pointsize 12→18）。
#   従来は generate_mono_figs_20260712.R が PNG を描き make_pdf 経由で PDF 化していたが、本版は基底図のみを
#   cairo_pdf で直接出力（自己完結）。計算は reassigned_basis_beta_results.rds から（再フィット不要）。
#   フォント拡大に合わせラベル自動配置の間隔 adj と右余白 pad を拡張（0.085→0.12, 0.13→0.17）。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
res <- readRDS("v3_fig/reassigned_basis_beta_results.rds")
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"

sty <- function(n){ list(
  col = c("black","black","gray45","black","gray45","gray60")[1:n],
  lty = c(1,      2,      1,       3,      2,       4)[1:n],
  lwd = c(2.0,    2.0,    2.3,     2.6,    2.3,     2.3)[1:n]) }

## 右端で線の脇に直接ラベル（縦の重なりを自動回避）。フォント拡大に合わせ adj を拡大。
rlab <- function(et, M, labels, cex=1.05, adj=0.12){
  k <- length(et); x0 <- et[k]; yv <- as.numeric(M[k, ])
  usr <- par("usr"); logy <- par("ylog"); ylo<-usr[3]; yhi<-usr[4]
  y <- if(logy) log10(yv) else yv; sep <- (yhi-ylo)*adj
  ord <- order(y); ys <- y[ord]
  for(i in seq_along(ys)[-1]) if(ys[i]-ys[i-1] < sep) ys[i] <- ys[i-1] + sep
  over <- max(ys)-(yhi-sep*0.4); if(over>0) ys <- ys-over
  under <- (ylo+sep*0.4)-min(ys); if(under>0) ys <- ys+under
  yo <- numeric(length(ys)); yo[ord] <- ys; yy <- if(logy) 10^yo else yo
  text(x0, yy, labels, pos=4, cex=cex, xpd=NA)
}
## 曲線脇ラベルは各プロット枠内の右側（データが無い余白領域）に収める → 右マージンを小さくでき、
## パネル間の余白が減り、かつラベルが枠外へはみ出さない。pad は最長ラベル（smallcell）が収まる幅。
xr <- function(et, pad=0.42) c(min(et), max(et)+pad*diff(range(et)))

draw <- function(){
  par(mfrow=c(2,3), mar=c(4.2,4.5,3,1.2))                              # 右マージンを詰める（ラベルは枠内右側に配置）
  for(nm in names(res)){ r<-res[[nm]]; s<-sty(2)
    matplot(r$event.times, r$X, type="l", lty=s$lty, lwd=s$lwd, col=s$col,
            xlab="event time", ylab="basis value (colSums=1)", xlim=xr(r$event.times),
            main=sprintf("(%s) shared non-negative basis X", nm), ylim=range(0,r$X))
    rlab(r$event.times, r$X, c("basis 1","basis 2")) }
  for(nm in names(res)){ r<-res[[nm]]; L<-ncol(r$beta.t); s<-sty(L)
    lab <- sub("^celltype","",colnames(r$beta.t))
    matplot(r$event.times, r$beta.t, type="l", lty=s$lty, lwd=s$lwd, col=s$col,
            xlab="event time", ylab=expression(paste(beta[l](t)," (per-SD log-HR)")), xlim=xr(r$event.times),
            main=sprintf("(%s) beta(t) of non-PH A (L=%d)", nm, L))
    abline(h=0, lty=3, col="gray70"); rlab(r$event.times, r$beta.t, lab) }
}
png("v3_fig/reassigned_basis_beta.png", width=2300, height=1450, res=170, type=ptype, pointsize=18); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/reassigned_basis_beta.pdf", width=2300/170, height=1450/170, pointsize=18); draw(); dev.off()
cat("saved v3_fig/reassigned_basis_beta.{png,pdf} (font 1.5x)\n=== done ===\n")
