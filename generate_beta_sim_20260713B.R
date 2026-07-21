# 【20260713B】fig:betasim / tab:betasim を (B) 版 mc_beta_beta500B57.rds（500reps, n=150/300/600/1200）から再描画。
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")
GAMMA0<-0.5
res <- readRDS("verification/mc_beta_beta500B57.rds")
ord <- intersect(c("n150","n300","n600","n1200"), names(res))
tab <- do.call(rbind, lapply(ord, function(tg){ M<-res[[tg]]$M
  data.frame(n=res[[tg]]$setting$n, bias=mean(M[,"cf"],na.rm=T)-GAMMA0, empSD=sd(M[,"cf"],na.rm=T),
    meanSE=mean(M[,"se"],na.rm=T), cov=mean(abs(M[,"cf"]-GAMMA0)<=1.96*M[,"se"],na.rm=T),
    ogap=mean(abs(M[,"cf"]-M[,"orc_fold"]),na.rm=T), brmse=mean(M[,"beta_rmse"],na.rm=T),
    bcov=mean(M[,"beta_cov"],na.rm=T)) }))
cat("=== MC-beta (B) summary (500 reps/n) ===\n"); print(format(tab,digits=3), row.names=FALSE)
cat("\n--- LaTeX rows (tab:betasim) ---\n")
for(i in 1:nrow(tab)) with(tab[i,], cat(sprintf("$%d$ & $%+.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.4f$ & $%.3f$ & $%.2f$ \\\\\n",
  n,bias,empSD,meanSE,cov,ogap,brmse,bcov)))
grDevices::cairo_pdf("v3_fig/beta_sim.pdf", width=9, height=3.7)
par(mfrow=c(1,2), mar=c(4.3,4.8,3,1))
plot(tab$n, tab$ogap, log="xy", type="b", pch=19, lwd=2, cex=1.4, xlab="sample size n",
     ylab=expression("mean |"*hat(gamma)-hat(gamma)^{"orac"}*"|"), main="(A) oracle gap vs n (log-log)")
ref<-tab$ogap[1]*sqrt(tab$n[1]/tab$n); lines(tab$n, ref, lty=2, lwd=1.6, col="gray50")
legend("topright", c("oracle gap (empirical)", expression("reference "*propto~n^{-1/2})), lty=c(1,2), pch=c(19,NA),
       lwd=2, col=c("black","gray50"), bty="n", cex=.95)
plot(tab$n, tab$cov, log="x", type="b", pch=19, lwd=2, cex=1.4, ylim=c(0.3,1), xlab="sample size n",
     ylab="95% coverage", main="(B) coverage vs n (log-x)")
lines(tab$n, tab$bcov, type="b", pch=1, lwd=2, cex=1.4, lty=2, col="gray35")
abline(h=0.95, lty=3, col="gray55")
legend("bottomleft", c(expression(gamma*" (inferential target)"), expression(beta(t)*" pointwise")),
       lty=c(1,2), pch=c(19,1), lwd=2, col=c("black","gray35"), bty="n", cex=1.0)
dev.off()
cat("\nsaved v3_fig/beta_sim.pdf (B)\n=== done ===\n")
