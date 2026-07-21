# 【2026-07-14 査読 r1 Major#2】Fig5 fig:mgus2hr（mgus2 連続非PH共変量の値別 HR(t)）を新設定で再描画。
#   全 complete-case (N=1338)・A={age,hgb,creat}（旧: 先頭700・{age,mspike}）。3パネル。フォント1.5倍（pointsize 18）。
#   reassigned_basis_beta_results.rds（generate_reassigned_basis_beta_20260714.R が生成）から HR を再計算（再フィット不要）。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
res <- readRDS("v3_fig/reassigned_basis_beta_results.rds")
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"

## 右端で曲線脇に直接ラベル（縦重なり自動回避）。枠内右側に収める。
rlab <- function(et, M, labels, cex=0.9, adj=0.11){
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
xr <- function(et, pad=0.42) c(min(et), max(et)+pad*diff(range(et)))   # 幅広ラベル（"10.60(p10)" 等）を枠内に収める

dm <- survival::mgus2; dm <- dm[complete.cases(dm[,c("futime","death","sex","age","hgb","creat","mspike")]),]
r <- res$mgus2; et <- r$event.times; qs <- c(.10,.25,.50,.75,.90)
ramp <- gray(seq(0.72,0,length=length(qs)))     # 薄い(p10)->濃い(p90)
covs <- c("age","hgb","creat")

draw <- function(){
  par(mfrow=c(1,3), mar=c(4.2,4.8,3,2.6))
  for(l in covs){
    m<-r$A.center[l]; sdl<-r$A.scale[l]; vals<-quantile(dm[[l]],qs,names=FALSE); z<-(vals-m)/sdl
    HR<-sapply(z, function(zz) exp(zz*r$beta.t[,l]))
    matplot(et, HR, type="l", lty=1, lwd=2.2, col=ramp, log="y", xlim=xr(et),
            xlab="event time", ylab=sprintf("HR(t) vs mean patient [%s]",l),
            main=sprintf("mgus2: time-varying HR by %s",l))
    abline(h=1, lty=3, col="gray55")
    rlab(et, HR, sprintf("%.2f(p%d)", vals, as.integer(qs*100))) }
}
png("v3_fig/mgus2_HR_by_value.png", width=2300, height=820, res=165, type=ptype, pointsize=18); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/mgus2_HR_by_value.pdf", width=2300/165, height=820/165, pointsize=18); draw(); dev.off()
cat("saved v3_fig/mgus2_HR_by_value.{png,pdf} (age/hgb/creat, N=1338, font 1.5x)\n=== done ===\n")
