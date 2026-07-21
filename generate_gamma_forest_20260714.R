# 【20260714】Fig8 fig:gamma（PH変数 γ の3手法比較 forest）を再描画。フォント1.5倍（pointsize 12→18）。
#   従来は PNG のみ→make_pdf 経由で PDF 化していたが、本版は cairo_pdf を直接出力（自己完結）。計算は rds から（再フィット不要）。
#   ① z-only Cox（灰・開円）② full-PH Cox（黒・開方）③ NMF-COX cross-fit（黒・塗円, 強調）。
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
tab <- readRDS("verification/gamma_table_20260714.rds")   # 全1338例・mgus2 z={sex,mspike}（査読 r1 Major#2）
pnum <- function(s) as.numeric(regmatches(s, gregexpr("[-+][0-9]+\\.[0-9]+", s))[[1]])  # c(point, low, high)

meth <- c("cox_z","cox_full","nmfcox")
mlab <- c("(1) Cox (z-only)","(2) full-PH Cox","(3) NMF-COX (cross-fit)")
mpch <- c(1, 0, 19); mcol <- c("gray45","black","black"); mlwd <- c(1.5,1.5,2.3); mcex <- c(1.1,1.1,1.25)
moff <- c(0.26, 0, -0.26)     # 各共変量スロット内の縦オフセット（①上 ②中 ③下）
datasets <- unique(tab$data)

draw <- function(){
  par(mfrow=c(1,length(datasets)), mar=c(4.4,6.2,3,1))
  for(di in seq_along(datasets)){ nm<-datasets[di]; sub<-tab[tab$data==nm,,drop=FALSE]
    covs<-sub$var; nc<-length(covs); ybase<-nc:1     # 上から下へ（covs[1]が上）
    allx<-unlist(lapply(seq_len(nc), function(i) unlist(lapply(meth, function(m) pnum(sub[i,m])[2:3]))))
    xr<-range(c(allx,0)); xr<-xr+c(-1,1)*0.06*diff(xr)
    plot(NA, xlim=xr, ylim=c(0.5,nc+0.6), yaxt="n", xlab=expression(hat(gamma)~"(log-HR)"), ylab="",
         main=sprintf("(%s)", nm))
    abline(v=0, lty=3, col="gray55")
    axis(2, at=ybase, labels=covs, las=1, cex.axis=.9)
    for(i in seq_len(nc)) for(j in seq_along(meth)){ v<-pnum(sub[i,meth[j]]); y<-ybase[i]+moff[j]
      arrows(v[2], y, v[3], y, angle=90, code=3, length=0.02, col=mcol[j], lwd=mlwd[j])
      points(v[1], y, pch=mpch[j], col=mcol[j], lwd=mlwd[j], cex=mcex[j], bg="white") }
    if(di==1) legend("bottomright", mlab, pch=mpch, col=mcol, pt.lwd=mlwd, bty="n", cex=.78)
  }
}
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
png("v3_fig/gamma_forest.png", width=2300, height=820, res=175, type=ptype, pointsize=18); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/gamma_forest.pdf", width=2300/175, height=820/175, pointsize=18); draw(); dev.off()
cat("saved v3_fig/gamma_forest.{png,pdf} (font 1.5x)\n=== done ===\n")
