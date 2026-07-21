# 【20260714】Fig9 rank_compare を再描画。フォント1.5倍（pointsize 12→18）。計算は rds から（再フィット不要）。
# 論文図: gbsg で NMF-COX(Q=2,3,4) vs RRR(coxvcf, 固定辞書 M=5 アンカ付き（=max rank, L=5） → 縮約ランク r=2,3,4) の β(t) 比較（モノクロ）。
#   【20260713r2】RRR を真の reduced rank に修正（M=4 固定スプライン辞書に対し縮約ランク r=2,3,4（全て真の縮約）。
#   上段=NMF-COX（罰則＋非負学習基底, Q=2,3,4）、下段=RRR（M=5 固定辞書 → 縮約 r=2,3,4, 全て真の縮約）。列=非PH共変量。
#   rank(Q/r) は線種で区別（=2 実線 / =3 破線 / =4 点線, 全黒）。結果 rds から再描画（再フィット不要）。
#   PNG＋ベクター PDF（cairo_pdf, フォント埋め込み）を出力。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
R <- readRDS("v3_fig/rank_compare_gbsg_results.rds")
et<-R$et; acov<-R$acov; Qs<-c(2,3,4)
bN<-R$beta_nmf; bR<-R$beta_rrr; rok<-as.integer(names(bR))
rng_nmf<-R$rng_nmf; rng_rrr<-R$rng_rrr; Mdim<-if(!is.null(R$rrr_Mdim)) R$rrr_Mdim else 4

## モノクロ rank スタイル（実線/破線/点線, 全て黒。色差を用いず線種のみで区別）: =2 実線 / =3 破線 / =4 点線
rlab<-as.character(c(2,3,4)); rcol<-c("black","black","black"); rlty<-c(1,2,3); rlwd<-c(2.0,2.2,2.4)
names(rcol)<-names(rlty)<-names(rlwd)<-rlab

draw <- function(){
  par(mfrow=c(2,length(acov)), mar=c(4,4.3,3,0.8), mgp=c(2.4,0.8,0))
  ## 各共変量（列）で NMF-COX と RRR の y 軸を共通（両者を合わせた最小・最大）にそろえる
  ysh <- lapply(acov, function(l) range(c(sapply(Qs, function(q) bN[[as.character(q)]][,l]),
                                          sapply(rok,function(r) bR[[as.character(r)]][l,]))))
  names(ysh) <- acov
  # 上段 NMF-COX（Q=2,3,4 の学習非負基底; K→Q 直接）
  for(l in acov){
    plot(NA,xlim=range(et),ylim=ysh[[l]],xlab="event time",ylab=expression(beta[l](t)~"(per-SD log-HR)"),
         main=sprintf("NMF-COX: %s", l)); abline(h=0,lty=3,col="gray70")
    for(q in Qs){ k<-as.character(q); lines(et, bN[[k]][,l], col=rcol[k], lty=rlty[k], lwd=rlwd[k]) }
    if(l==acov[1]) legend("topleft", sprintf("Q=%d",Qs), col=rcol, lty=rlty, lwd=rlwd, bty="n", cex=.85, seg.len=4) }
  # 下段 RRR（固定辞書 M=Mdim → 縮約ランク r=2,3,4; K→M→r）
  for(l in acov){
    plot(NA,xlim=range(et),ylim=ysh[[l]],xlab="event time",ylab=expression(beta[l](t)~"(per-SD log-HR)"),
         main=sprintf("RRR (M=%d anch.): %s", Mdim, l)); abline(h=0,lty=3,col="gray70")
    for(r in rok){ k<-as.character(r); lines(et, bR[[k]][l,], col=rcol[k], lty=rlty[k], lwd=rlwd[k]) }
    if(l==acov[1]) legend("topleft", sprintf("r=%d",rok), col=rcol, lty=rlty, lwd=rlwd, bty="n", cex=.85, seg.len=4) }
}

ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
# PNG（モノクロ）
png("v3_fig/rank_compare_gbsg.png", width=520*length(acov), height=1060, res=140, type=ptype, pointsize=18); draw(); dev.off()
# PDF（ベクター, フォント埋め込み）
grDevices::cairo_pdf("v3_fig/rank_compare_gbsg.pdf", width=520*length(acov)/140, height=1060/140, pointsize=18); draw(); dev.off()

cat("=== rank instability table (max pointwise range across ranks; small=stable) ===\n")
print(round(rbind(NMF_COX=rng_nmf, RRR=rng_rrr),3))
cat(sprintf("most stable (NMF-COX): %s (%.3f);  NMF/RRR mean range = %.3f / %.3f\n",
    names(which.min(rng_nmf)), min(rng_nmf), mean(rng_nmf), mean(rng_rrr)))
cat("saved v3_fig/rank_compare_gbsg.{png,pdf} (mono)\n=== done ===\n")
