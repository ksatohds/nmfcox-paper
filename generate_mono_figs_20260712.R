# Figure 1-3 をモノクロ化（黒・灰＋線種）。凡例をまとめず、線の脇に直接ラベル（右端、縦方向は自動 de-collision）。
#   結果 rds から再計算（再フィット不要）。上書き先は既存 PNG（カラー版は *_color.png にバックアップ済）。
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
res <- readRDS("v3_fig/reassigned_basis_beta_results.rds")
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"

## モノクロ線スタイル（n 本まで区別）: 黒/灰 × 実線/破線/点線/一点鎖線
sty <- function(n){ list(
  col = c("black","black","gray45","black","gray45","gray60")[1:n],
  lty = c(1,      2,      1,       3,      2,       4)[1:n],
  lwd = c(2.0,    2.0,    2.3,     2.6,    2.3,     2.3)[1:n]) }

## 右端で線の脇に直接ラベル（縦の重なりを自動回避）。xlim は右に余白を確保しておくこと。
rlab <- function(et, M, labels, cex=1.05, adj=0.085){
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
xr <- function(et, pad=0.13) c(min(et), max(et)+pad*diff(range(et)))   # 右余白付き xlim（ラベル拡大に伴い縮小）

## ============ Figure 1: 上=基底X, 下=β(t) （veteran/gbsg/mgus2）============
png("v3_fig/reassigned_basis_beta.png", width=2300, height=1450, res=170, type=ptype)
par(mfrow=c(2,3), mar=c(4.2,4.4,3,3.6))
for(nm in names(res)){ r<-res[[nm]]; s<-sty(2)
  matplot(r$event.times, r$X, type="l", lty=s$lty, lwd=s$lwd, col=s$col,
          xlab="event time", ylab="basis value (colSums=1)", xlim=xr(r$event.times,0.13),
          main=sprintf("(%s) shared non-negative basis X", nm), ylim=range(0,r$X))
  rlab(r$event.times, r$X, c("basis 1","basis 2")) }
for(nm in names(res)){ r<-res[[nm]]; L<-ncol(r$beta.t); s<-sty(L)
  lab <- sub("^celltype","",colnames(r$beta.t))
  matplot(r$event.times, r$beta.t, type="l", lty=s$lty, lwd=s$lwd, col=s$col,
          xlab="event time", ylab=expression(paste(beta[l](t)," (per-SD log-HR)")), xlim=xr(r$event.times),
          main=sprintf("(%s) beta(t) of non-PH A (L=%d)", nm, L))
  abline(h=0, lty=3, col="gray70"); rlab(r$event.times, r$beta.t, lab) }
dev.off()

## ============ Figure 2: mgus2 連続量の値別 HR(t)（灰の濃淡ランプ p10->p90）============
dm <- survival::mgus2; dm <- dm[complete.cases(dm[,c("futime","death","sex","age","hgb","creat","mspike")]),]
if(nrow(dm)>700) dm <- dm[1:700,]
r <- res$mgus2; et <- r$event.times; qs <- c(.10,.25,.50,.75,.90)
ramp <- gray(seq(0.72,0,length=length(qs)))     # 薄い(p10)->濃い(p90)
png("v3_fig/mgus2_HR_by_value.png", width=2000, height=900, res=165, type=ptype)
par(mfrow=c(1,2), mar=c(4.2,4.6,3,5.6))
for(l in c("age","mspike")){
  m<-r$A.center[l]; sdl<-r$A.scale[l]; vals<-quantile(dm[[l]],qs,names=FALSE); z<-(vals-m)/sdl
  HR<-sapply(z, function(zz) exp(zz*r$beta.t[,l]))
  matplot(et, HR, type="l", lty=1, lwd=2.2, col=ramp, log="y", xlim=xr(et,0.20),
          xlab="event time", ylab=sprintf("HR(t) vs mean patient [%s]",l),
          main=sprintf("(mgus2) time-varying HR by %s",l))
  abline(h=1, lty=3, col="gray55")
  rlab(et, HR, sprintf("%.2f(p%d)", vals, as.integer(qs*100)), cex=.78) }
dev.off()

## ============ Figure 3: ダミー 1-vs-参照 HR(t)============
hr_dummy <- function(r,cols) sapply(cols, function(l) exp(r$beta.t[,l]/r$A.scale[l]))
vet<-res$veteran; vd<-grep("^celltype",colnames(vet$beta.t),value=TRUE); gb<-res$gbsg
png("v3_fig/dummy_hr_by_covariate.png", width=2000, height=880, res=165, type=ptype)
par(mfrow=c(1,2), mar=c(4.2,4.8,3,5.6))
HRv<-hr_dummy(vet,vd); s<-sty(length(vd))
matplot(vet$event.times, HRv, type="l", lty=s$lty, lwd=s$lwd, col=s$col, log="y", xlim=xr(vet$event.times,0.20),
        xlab="event time", ylab="HR(t) vs squamous", main="(veteran) celltype: HR(t) vs squamous")
abline(h=1, lty=3, col="gray55"); rlab(vet$event.times, HRv, sub("^celltype","",vd), cex=.82)
HRg<-hr_dummy(gb,"meno")
matplot(gb$event.times, HRg, type="l", lty=1, lwd=2.4, col="black", log="y", xlim=xr(gb$event.times,0.20),
        xlab="event time", ylab="HR(t): post vs pre", main="(gbsg) meno: HR(t) post vs pre")
abline(h=1, lty=3, col="gray55"); rlab(gb$event.times, cbind(HRg), "meno", cex=.82)
dev.off()

## ============ Figure 4: 学習基底 vs 与える基底 fit-vs-df（モノクロ）============
fd <- readRDS("v3_fig/compare_rrr_basis_count_results.rds")
png("v3_fig/compare_rrr_basis_count.png", width=2300, height=780, res=170, type=ptype)
par(mfrow=c(1,length(fd)), mar=c(4.4,4.6,3,1))
for(nm in names(fd)){ tb<-fd[[nm]]; rr<-tb[tb$model=="RRR",]; nn<-tb[tb$model!="RRR",]
  plot(rr$df, rr$loglik, type="b", pch=19, col="gray60", lwd=2, cex=1.1,   # RRR 折れ線＝灰色（テキストを前面に）
       xlab="effective df", ylab="Breslow partial loglik", main=sprintf("(%s) fit vs df", nm),
       xlim=range(c(rr$df,nn$df)), ylim=range(c(rr$loglik,nn$loglik)))
  text(rr$df, rr$loglik, labels=paste0("Q=",rr$q), pos=1, cex=1.1, font=2, col="black")  # 灰線の上に黒テキスト
  points(nn$df, nn$loglik, pch=4, col="black", cex=1.7, lwd=2.6)
  text(nn$df, nn$loglik, "NMF-COX(Q=2)", pos=3, cex=1.05, font=2, col="black")
  legend("bottomright", c("RRR (spline, rank2): gray filled circle","NMF-COX (learned basis): cross"),
         pch=c(19,4), col=c("gray60","black"), bty="n", cex=1.0) }
dev.off()
cat("saved mono: reassigned_basis_beta.png, mgus2_HR_by_value.png, dummy_hr_by_covariate.png, compare_rrr_basis_count.png\n=== done ===\n")
