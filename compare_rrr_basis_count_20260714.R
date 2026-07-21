# 「基底を与える RRR vs 基底を学習する NMF-COX」の定量比較。
#   同じ Breslow 部分尤度で両者の当てはまりを評価（タイ処理の交絡を除去）。
#   RRR: rank=2 固定, 時間関数の数 q を 2..6 で増やす（q=2 が NMF-COX の基底数 Q=2 と同数）。
#   NMF-COX: Q=2（基底は自動学習）。df は AIC 機構（df.eff）を使用。
#   予想: q=2(同数)では RRR の当てはまりが劣り、NMF-COX に追いつくには大きい q（＝より多くの基底）が要る。
suppressMessages({library(survival); library(MASS); library(splines)})
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
suppressMessages(library(nmfkc))
## RRR baseline uses Perperoglou's reduced-rank code (third-party; NOT redistributed here).
## Download coxvc.r & coxvcf.r from https://github.com/drperpo and place them in the repository root.
if(!all(file.exists(c("coxvc.r","coxvcf.r"))))
  stop("This RRR-comparison figure needs coxvc.r & coxvcf.r (Perperoglou). See README / https://github.com/drperpo .")
source("coxvc.r"); source("coxvcf.r")

# 【2026-07-14 査読 r1 Major#2】mgus2 を全 complete-case (N=1338)・cox.zph 準拠 split（z=sex+mspike, a=age+hgb+creat）に更新。
DS <- list(
  veteran = list(tvar="time", evar="status", zf="trt+age+diagtime+prior", af="karno+celltype", SM=1000),
  gbsg    = list(tvar="rfstime", evar="status", zf="hormon+nodes+size", af="age+meno+grade+pgr+er", SM=3000),
  mgus2   = list(tvar="futime", evar="death", zf="sex+mspike", af="age+hgb+creat", SM=3000))
getdata <- function(nm){
  if(nm=="veteran") d<-survival::veteran
  if(nm=="gbsg")    d<-survival::gbsg
  if(nm=="mgus2"){ d<-survival::mgus2; d<-d[complete.cases(d[,c("futime","death","sex","age","hgb","creat","mspike")]),] }
  d
}
## 共通 Breslow 部分対数尤度: η_i(t_j)=z_i'α + a_i'β(t_j)。betaK=K×L(β at et), Amat=L×n(そのモデルの尺度)。
breslow_ll <- function(time, status, za, Amat, betaK, et){
  ll <- 0
  for(j in seq_along(et)){ tj<-et[j]
    risk <- which(time>=tj); fail <- which(status==1 & time==tj)
    if(!length(fail) || !length(risk)) next
    bt <- betaK[j,]
    eta <- za[risk] + as.vector(crossprod(Amat[,risk,drop=FALSE], bt))
    ef  <- za[fail] + as.vector(crossprod(Amat[,fail,drop=FALSE], bt))
    m <- max(eta); ll <- ll + sum(ef) - length(fail)*(m + log(sum(exp(eta-m))))
  }
  ll
}

RANK <- 2; QSEQ <- c(2,3,4,6)   # RRR の時間関数数 q（q=2 が NMF の基底数と同数）
out <- list()
for(nm in names(DS)){ p<-DS[[nm]]; d<-getdata(nm); d<-d[order(d[[p$tvar]]),]
  time<-d[[p$tvar]]; status<-as.integer(d[[p$evar]]); n<-nrow(d)
  et <- sort(unique(time[status==1])); K<-length(et)
  Zmat <- model.matrix(as.formula(paste0("~",p$zf)), data=d)[,-1,drop=FALSE]
  Amat0 <- t(model.matrix(as.formula(paste0("~",p$af)), data=d)[,-1,drop=FALSE])   # L×n (raw)
  L<-nrow(Amat0); pz<-ncol(Zmat)
  ## --- NMF-COX (Q=2, 基底自動) ---
  d$y<-Surv(time,status)
  fN <- nmf.cox(as.formula(paste0("y~",p$zf)), data=d, A=as.formula(paste0("~",p$af)),
                 rank=RANK, X.L2.smooth=p$SM, ties="breslow", maxit=30, verbose=FALSE,
                 inference=FALSE, aic=TRUE, x.update="newton", seed=1)
  aN <- fN$gamma[colnames(Zmat)]; zaN <- as.vector(Zmat %*% aN)
  llN <- breslow_ll(time, status, zaN, fN$A, fN$beta.t, et)     # A は標準化, beta.t は per-SD
  cat(sprintf("[%s] check: breslow_ll(NMF)=%.2f vs fit$loglik=%.2f (差=%.3f)\n", nm, llN, fN$loglik, llN-fN$loglik))
  rows <- list(data.frame(model="NMF-COX(Q=2)", q=NA, loglik=llN, df=fN$df.eff, aic=-2*llN+2*fN$df.eff))
  ## --- RRR (rank=2, q=2..6, 基底=固定スプライン) ---
  Acen <- Amat0 - rowMeans(Amat0)                                # coxvc は中心化のみ
  tm <- mean(time); ts <- sd(time)
  for(q in QSEQ){
    if(q>=3){ nsb <- ns(time, df=q-1); Ft <- cbind(1, nsb); Fe <- cbind(1, predict(nsb, et)) }
    else    { Ft <- cbind(1, (time-tm)/ts); Fe <- cbind(1, (et-tm)/ts) }   # q=2: 定数+線形
    set.seed(1)
    fR <- tryCatch(coxvcf(as.formula(paste0("Surv(",p$tvar,",",p$evar,")~",p$zf)),
                          X=t(Amat0), Ft=Ft, rank=RANK, data=d),
                   error=function(e){cat("  q=",q," ERR:",conditionMessage(e),"\n");NULL})
    if(is.null(fR)) next
    betaK <- t(fR$theta %*% t(Fe))                               # K×L
    aR <- fR$f.coef; zaR <- as.vector(Zmat %*% aR)
    llR <- breslow_ll(time, status, zaR, Acen, betaK, et)
    dfR <- RANK*(L + ncol(Ft) - RANK) + pz                       # r(p+q-r)+v
    rows[[length(rows)+1]] <- data.frame(model="RRR", q=q, loglik=llR, df=dfR, aic=-2*llR+2*dfR)
  }
  tab <- do.call(rbind, rows); out[[nm]] <- tab
  cat(sprintf("=== %s (n=%d, K=%d, L=%d) ===\n", nm, n, K, L)); print(format(tab, digits=5)); cat("\n")
}
saveRDS(out, "v3_fig/compare_rrr_basis_count_results.rds")

## 図: 各データ Breslow loglik vs df（RRR の q 系列 ● と NMF-COX ★）
ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
png("v3_fig/compare_rrr_basis_count.png", width=2300, height=780, res=170, type=ptype)
par(mfrow=c(1,length(out)), mar=c(4.4,4.6,3,1))
for(nm in names(out)){ tb<-out[[nm]]
  rr<-tb[tb$model=="RRR",]; nn<-tb[tb$model!="RRR",]
  plot(rr$df, rr$loglik, type="b", pch=19, col="#D55E00", lwd=2,
       xlab="effective df", ylab="Breslow partial loglik",
       main=sprintf("(%s) fit vs df", nm),
       xlim=range(c(rr$df,nn$df)), ylim=range(c(rr$loglik,nn$loglik)))
  text(rr$df, rr$loglik, labels=paste0("q=",rr$q), pos=1, cex=.75, col="#D55E00")
  points(nn$df, nn$loglik, pch=8, col="black", cex=1.8, lwd=2)
  text(nn$df, nn$loglik, "NMF-COX\nQ=2", pos=3, cex=.8, font=2)
  legend("bottomright", c("RRR (spline, rank2)","NMF-COX (learned basis)"),
         pch=c(19,8), col=c("#D55E00","black"), bty="n", cex=.8) }
dev.off()
cat("saved v3_fig/compare_rrr_basis_count.png + results.rds\n=== done ===\n")
