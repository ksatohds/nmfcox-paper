# 【20260714】fig:rrrfit の再描画。フォント1.5倍（pointsize 12→18）。ラベルは M=（RRR辞書サイズ, r=2固定）。
# fig:rrrfit（学習基底 vs 与える基底, fit vs 有効df）の再描画（モノクロ PNG + ベクター PDF）。
#   【20260713】表記統一：RRR の点ラベルを「Q=」→「M=」に変更（M=辞書サイズ, 縮約ランク r=2 固定）。
#   計算は不変（`compare_rrr_basis_count_results.rds` をそのまま使用；再フィット不要）。generate_mono_figs_20260712.R の
#   Figure 4 ブロックと同一だが点ラベルのみ M= に統一。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
fd <- readRDS("v3_fig/compare_rrr_basis_count_results.rds")
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"

draw <- function(){
  par(mfrow=c(1,length(fd)), mar=c(4.4,4.6,3,1))
  for(nm in names(fd)){ tb<-fd[[nm]]; rr<-tb[tb$model=="RRR",]; nn<-tb[tb$model!="RRR",]
    xrp<-range(c(rr$df,nn$df)); xrp<-xrp+c(-1,1)*0.11*diff(xrp)          # 端ラベルのクリップ回避（フォント拡大対応）
    yrp<-range(c(rr$loglik,nn$loglik)); yrp<-yrp+c(-1,1)*0.07*diff(yrp)
    plot(rr$df, rr$loglik, type="b", pch=19, col="gray60", lwd=2, cex=1.1,
         xlab="effective df", ylab="Breslow partial loglik", main=sprintf("(%s) fit vs df", nm),
         xlim=xrp, ylim=yrp)
    text(rr$df, rr$loglik, labels=paste0("M=",rr$q), pos=1, cex=1.1, font=2, col="black")   # M=辞書サイズ（r=2 固定）
    points(nn$df, nn$loglik, pch=4, col="black", cex=1.7, lwd=2.6)
    text(nn$df, nn$loglik, "NMF-COX(Q=2)", pos=3, cex=1.05, font=2, col="black")
    legend("topleft", c("RRR (spline, rank r=2): gray filled circle","NMF-COX (learned basis): cross"),
           pch=c(19,4), col=c("gray60","black"), bty="n", cex=0.95) }
}

png("v3_fig/compare_rrr_basis_count.png", width=2300, height=780, res=170, type=ptype, pointsize=18); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/compare_rrr_basis_count.pdf", width=2300/170, height=780/170, pointsize=18); draw(); dev.off()
cat("saved v3_fig/compare_rrr_basis_count.{png,pdf} (labels: M= for RRR dictionary size, r=2 fixed)\n=== done ===\n")
