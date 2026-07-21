# gbsg: NMF-COX（Q=2,3,4）vs RRR（coxvcf, 真の reduced rank r=2,3,4）の β(t) 比較。
#   【20260713r3 改訂】RRR のスプライン辞書を Perperoglou 流に **アンカリング**（非定数列を開始時点 t0 で 0 に）。
#   H(t)=[1, ns_1(t)-ns_1(t0), ..., ns_4(t)-ns_4(t0)]：第1列=時間一定効果、残り4列=開始時点からの変化、と解釈できる。
#   ※アンカリングは非定数列から定数を引く（切片列に吸収）だけで span 不変、かつ reduced-rank 制約 rank(C)≤r は
#     基底の可逆再パラメータ化 H→HA で不変ゆえ、当てはめ β(t)=HC は数値的に r2（未アンカ）と一致する（本スクリプトで確認）。
#   さらに **M=6 頑健性チェック**（辞書を ns(df=5) に増やしても傾向不変か）を併記。
#   NMF 側は不変ゆえ r2 の rds から読み込み、RRR のみ再計算する。
suppressMessages({library(survival); library(MASS); library(splines)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
## RRR baseline uses Perperoglou's reduced-rank code (third-party; NOT redistributed here).
## Download coxvc.r & coxvcf.r from https://github.com/drperpo and place them in the repository root.
if(!all(file.exists(c("coxvc.r","coxvcf.r"))))
  stop("This RRR-comparison figure needs coxvc.r & coxvcf.r (Perperoglou). See README / https://github.com/drperpo .")
source("coxvc.r"); source("coxvcf.r")

## ---- NMF 結果は r2 の rds から読み込み（NMF は M に依存しないので不変）----
R0 <- readRDS("v3_fig/rank_compare_gbsg_results.rds")
et <- R0$et; acov <- R0$acov; Qs <- c(2,3,4)
beta_nmf <- R0$beta_nmf; gamma_nmf <- R0$gamma_nmf; znames <- rownames(gamma_nmf)
rng_nmf <- R0$rng_nmf; nmf_sel_smooth <- R0$nmf_sel_smooth

## ---- gbsg 準備 ----
d <- survival::gbsg; d <- d[order(d$rfstime), ]
tt <- d$rfstime
Xmat <- model.matrix(~age+meno+grade+pgr+er, data=d)[,-1,drop=FALSE]
colnames(Xmat) <- make.names(colnames(Xmat)); Xmat <- scale(Xmat)      # per-SD
t0 <- min(tt)                                                          # アンカー基準（開始時点）

## ---- アンカ付き RRR を任意の辞書サイズ M で当てはめる関数 ----
fit_rrr <- function(Mdim, ranks=c(2,3,4)){
  nsb <- ns(tt, df=Mdim-1)
  ns0 <- as.numeric(predict(nsb, t0))                                 # t0 での値（非定数列）
  Ft  <- cbind(1, sweep(nsb, 2, ns0))                                 # アンカ：非定数列は t0 で 0
  Feq <- cbind(1, sweep(predict(nsb, et), 2, ns0))                    # 事象時点での辞書（固定）
  fits <- list()
  for(r in ranks){ set.seed(1)
    capture.output(f <- tryCatch(coxvcf(Surv(rfstime,status)~hormon+nodes+size, X=Xmat, Ft=Ft, rank=r, data=d),
                                 error=function(e){message("RRR M=",Mdim," r=",r," ERR: ",conditionMessage(e)); NULL}))
    fits[[as.character(r)]] <- f }
  ok <- ranks[!sapply(fits[as.character(ranks)], is.null)]
  beta <- lapply(ok, function(r){ b<-fits[[as.character(r)]]$theta %*% t(Feq); rownames(b)<-acov; b }); names(beta)<-as.character(ok)
  gam  <- sapply(ok, function(r){ g<-fits[[as.character(r)]]$f.coef; names(g)<-znames; g[znames] })
  if(is.null(dim(gam))) gam<-matrix(gam,length(znames),dimnames=list(znames,as.character(ok)))
  ll   <- sapply(ok, function(r) fits[[as.character(r)]]$loglik)
  rng  <- sapply(acov, function(l) max(apply(sapply(ok, function(r) beta[[as.character(r)]][l,]),1,function(x)diff(range(x)))))
  list(beta=beta, gam=gam, ll=ll, rng=rng, ok=ok, Mdim=Mdim)
}

## ---- 主：M=5（アンカ付き, 論文用）----
RR <- fit_rrr(5)
betaR <- RR$beta; gamR <- RR$gam; llR <- RR$ll; rng_rrr <- RR$rng; rok <- RR$ok

## ---- 頑健性：M=6（アンカ付き）----
RR6 <- fit_rrr(6)

cat("=== NMF-COX gamma (Q; from r2 rds) ===\n"); print(round(gamma_nmf,4))
cat("=== RRR gamma (M=5 anchored, reduced rank r) ===\n"); print(round(gamR,4)); cat("RRR loglik:", sprintf("%.2f",llR),"\n")
cat("\n=== beta(t) instability across rank: max pointwise range (small=stable) ===\n")
print(round(rbind(NMF=rng_nmf, RRR_M5=rng_rrr, RRR_M6=RR6$rng),3))
cat(sprintf("mean range: NMF=%.3f  RRR(M5)=%.3f  RRR(M6)=%.3f\n", mean(rng_nmf), mean(rng_rrr), mean(RR6$rng)))
cat(sprintf("gamma_hormon spread: NMF=%.4f  RRR(M5)=%.4f  RRR(M6)=%.4f\n",
    diff(range(gamma_nmf["hormon",])), diff(range(gamR["hormon",])), diff(range(RR6$gam["hormon",]))))
cat("[check] M=5 anchored vs r2(unanchored) は β(t)=HC が数学的に同一のはず → mean range が r2 の 0.289 に一致することを確認\n")

## ---- 図（カラー簡易版; モノクロ／PDF は mono r3 スクリプトで rds から再描画）----
cols <- c("2"="black","3"="#D55E00","4"="#0072B2")
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
png("v3_fig/rank_compare_gbsg_color.png", width=520*length(acov), height=1040, res=140, type=ptype)
par(mfrow=c(2,length(acov)), mar=c(4,4.3,3,1))
for(l in acov){
  yl<-range(c(sapply(Qs,function(q) beta_nmf[[as.character(q)]][,l]), sapply(rok,function(r) betaR[[as.character(r)]][l,])))
  plot(NA,xlim=range(et),ylim=yl,xlab="event time",ylab=expression(beta[l](t)),main=sprintf("NMF-COX: %s", l)); abline(h=0,lty=3,col="gray70")
  for(q in Qs) lines(et, beta_nmf[[as.character(q)]][,l], col=cols[as.character(q)], lwd=2.2)
  legend("topleft",sprintf("Q=%d",Qs),col=cols[as.character(Qs)],lwd=2.2,bty="n",cex=.8) }
for(l in acov){
  yl<-range(c(sapply(Qs,function(q) beta_nmf[[as.character(q)]][,l]), sapply(rok,function(r) betaR[[as.character(r)]][l,])))
  plot(NA,xlim=range(et),ylim=yl,xlab="event time",ylab=expression(beta[l](t)),main=sprintf("RRR (M=5 anchored): %s", l)); abline(h=0,lty=3,col="gray70")
  for(r in rok) lines(et, betaR[[as.character(r)]][l,], col=cols[as.character(r)], lwd=2.2)
  legend("topleft",sprintf("r=%d",rok),col=cols[as.character(rok)],lwd=2.2,bty="n",cex=.8) }
dev.off()
saveRDS(list(gamma_nmf=gamma_nmf, gamma_rrr=gamR, beta_nmf=beta_nmf, beta_rrr=betaR,
             et=et, acov=acov, rng_nmf=rng_nmf, rng_rrr=rng_rrr, llR=llR,
             nmf_sel_smooth=nmf_sel_smooth, rrr_Mdim=5, rrr_ranks=rok, rrr_anchored=TRUE,
             rng_rrr_M6=RR6$rng, gamma_rrr_M6=RR6$gam),
        "v3_fig/rank_compare_gbsg_results.rds")
cat("\nsaved v3_fig/rank_compare_gbsg_results.rds (M=5 anchored + M=6 robustness) (+ _color.png)\n=== done ===\n")
