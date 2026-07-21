# cox.zph 再割り当て版（A=真に非PH）の解釈図: 上段=共有非負基底 X（2本）、
#   下段=A 各共変量の時間変動係数 β(t)=XΘ（標準化スケール=per-SD log-HR）。veteran/gbsg/mgus2。
#   （2026-07-13: 現行 nmf_cox3.R の基底列並べ替え〔時間重心 Σ(p/P)X 昇順, §6.12〕を反映して rds/図を再生成。
#    以前の rds/図は §6.12 の並べ替え導入前で、veteran/gbsg は basis1=後期側に付番されていた。
#    追加で (i) 再フィットの β(t) が旧 rds と一致すること（＝並べ替えは順序のみで β(t)=XΘ 不変）を検証、
#    (ii) seed を振った基底の安定性〔列マッチ後の相関〕を診断する。)
suppressMessages(library(survival))
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc)); dir.create("v3_fig", showWarnings=FALSE)

DS <- list()
{ d<-survival::veteran; d$y<-Surv(d$time,d$status)
  DS$veteran<-list(d=d, zf=y~trt+age+diagtime+prior, Af=~karno+celltype, SM=1000) }
{ d<-survival::gbsg; d$y<-Surv(d$rfstime,d$status)
  DS$gbsg<-list(d=d, zf=y~hormon+nodes+size, Af=~age+meno+grade+pgr+er, SM=3000) }
# mgus2 (血液内科: 単クローン性ガンマグロブリン血症)。【2026-07-14 査読 r1 Major#2】全 complete-case (N=1338) を使用（旧: 先頭700 は順序依存の恣意選択）。
#   全例で cox.zph 機械割当: 非PH(p<0.05)={age,hgb,creat}→A, PH={sex,mspike}→z（先頭700 の {age,mspike} は subset 依存だった）。
{ d<-survival::mgus2; d<-d[complete.cases(d[,c("futime","death","sex","age","hgb","creat","mspike")]),]
  d$y<-Surv(d$futime,d$death)
  DS$mgus2<-list(d=d, zf=y~sex+mspike, Af=~age+hgb+creat, SM=3000) }   # SM は CVPL で確認（2026-07-14）

fit_one <- function(p, seed) nmf.cox(p$zf, data=p$d, A=p$Af, rank=2, X.L2.smooth=p$SM, ties="breslow",
                                      maxit=30, verbose=FALSE, inference=FALSE, x.update="newton", seed=seed)
## 明示的な多スタート: seed=1..NSTART で当てはめ、罰則付き部分尤度（objfunc）最良の解を採用。
##   基底が高相関で弱識別のデータ（gbsg, basis.cor≈0.91）で稀に生じる劣位局所解を排除し、図を再現可能にする。
NSTART <- 8
fit_best <- function(p){
  best <- NULL
  for(sd in 1:NSTART){ f <- fit_one(p, sd)
    if(is.null(best) || (!is.na(f$objfunc) && f$objfunc < best$objfunc)) best <- f }   # objfunc は小さいほど良い（=負の罰則付き部分尤度）
  best
}

pal <- c("black","#D55E00","#0072B2","#009E73","#CC79A7","#E69F00")
fits <- list()
for(nm in names(DS)){ p<-DS[[nm]]
  f <- fit_best(p); fits[[nm]] <- f
  cen <- as.vector((seq_len(nrow(f$X))/nrow(f$X)) %*% f$X)   # 並べ替えキー（列ごとの時間重心, 小=早期）
  cat(sprintf("[%s] A covariates (L=%d): %s | objfunc=%.4f | basis cor=%+.3f | centroid(basis1,basis2)=(%.3f,%.3f) %s\n",
              nm, ncol(f$beta.t), paste(colnames(f$beta.t),collapse=", "), f$objfunc, cor(f$X[,1],f$X[,2]),
              cen[1], cen[2], if(cen[1]<=cen[2]) "[basis1=早期 OK]" else "[!! basis1 が後期]"))
}

## (i) 検証: 再フィット β(t) が旧 rds（並べ替え前）と一致するか（並べ替えは順序置換のみ→ β(t)=XΘ 不変のはず）
old <- tryCatch(readRDS("v3_fig/reassigned_basis_beta_results_20260712.rds"), error=function(e) NULL)
if(!is.null(old)){
  cat("\n--- (i) β(t) invariance vs pre-reorder rds (max abs diff over covariates) ---\n")
  for(nm in names(fits)){
    if(is.null(old[[nm]]) || !identical(dim(fits[[nm]]$beta.t), dim(old[[nm]]$beta.t))){
      cat(sprintf("  %-8s skipped (A config changed since 2026-07-12 rds; e.g. mgus2 A now {age,hgb,creat})\n", nm)); next }  # 診断のみ; 図生成には無関係
    d <- max(abs(fits[[nm]]$beta.t - old[[nm]]$beta.t))
    cat(sprintf("  %-8s max|Δβ(t)| = %.3e  %s\n", nm, d, if(d<1e-6) "[unchanged: order-only]" else "[!! fit changed]"))
  }
}

## (ii) seed 安定性: 各 seed の解を「多スタート最良解」と比較。
##   ・objfunc が最良に一致する再スタートの割合（＝支配的 basin の再現率）
##   ・基底 X の列相関（弱識別なら劣位 basin で低下しうる）
##   ・β(t)=XΘ の列相関（同定される積＝どの basin でも安定なはず）
cat("\n--- (ii) restart stability vs multi-start best (basis X and identified beta(t)) ---\n")
seeds <- 1:12
for(nm in names(DS)){ p<-DS[[nm]]
  best<-fits[[nm]]; Xref<-best$X; Bref<-best$beta.t; obest<-best$objfunc
  nbest<-0; xmin<-c(); bmin<-c()
  for(sd in seeds){ f2 <- fit_one(p, sd)
    if(!is.na(f2$objfunc) && abs(f2$objfunc-obest) < 1e-3) nbest<-nbest+1
    xmin <- c(xmin, min(cor(Xref[,1],f2$X[,1]), cor(Xref[,2],f2$X[,2])))
    bmin <- c(bmin, min(sapply(seq_len(ncol(Bref)), function(j) cor(Bref[,j],f2$beta.t[,j])))) }
  cat(sprintf("  %-8s best-basin %d/%d restarts | min basis-corr=%.3f | min beta(t)-corr=%.4f  %s\n",
              nm, nbest, length(seeds), min(xmin), min(bmin),
              if(min(bmin)>0.99) "[beta(t) stable]" else "[!! beta(t) varies]"))
}

# 軽量な結果のみ保存（両図の再描画に必要な成分：full fit を保存すると 22MB になるため）
saveRDS(lapply(fits, function(f) list(
          event.times = f$event.times, X = f$X, beta.t = f$beta.t,
          A.center = f$A.center, A.scale = f$A.scale,
          anames = colnames(f$beta.t), basis.cor = cor(f$X[,1], f$X[,2]))),
        "v3_fig/reassigned_basis_beta_results.rds")

ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"

## ===== Fig A（カラー版; 公表版は generate_mono_figs で上書き。ここでは点検用）=====
png("v3_fig/reassigned_basis_beta.png", width=2300, height=1450, res=170, type=ptype)
par(mfrow=c(2,3), mar=c(4.2,4.4,3,1))
for(nm in names(fits)){ r<-fits[[nm]]; et<-r$event.times
  matplot(et, r$X, type="l", lty=1, lwd=2.6, col=c("black","gray55"),
          xlab="event time", ylab="basis value (colSums=1)",
          main=sprintf("(%s) shared non-negative basis X", nm), ylim=range(0,r$X))
  legend("topright", c("basis 1","basis 2"), lty=1, lwd=2.6, col=c("black","gray55"), bty="n", cex=.85) }
for(nm in names(fits)){ r<-fits[[nm]]; et<-r$event.times; L<-ncol(r$beta.t); cc<-pal[((seq_len(L)-1)%%length(pal))+1]
  matplot(et, r$beta.t, type="l", lty=1, lwd=2.2, col=cc,
          xlab="event time", ylab=expression(paste(beta[l](t), "  (per-SD log-HR)")),
          main=sprintf("(%s) beta(t) of non-PH A  (L=%d)", nm, L))
  abline(h=0, lty=3, col="gray70")
  legend("topright", colnames(r$beta.t), lty=1, lwd=2.2, col=cc, bty="n", cex=.72) }
dev.off()

cat("\nsaved v3_fig/reassigned_basis_beta.png + reassigned_basis_beta_results.rds (reordered).",
    "\nNext: source generate_mono_figs_20260712.R then make_pdf_figs_20260712.R to refresh mono PNG/PDF.\n=== done ===\n")
