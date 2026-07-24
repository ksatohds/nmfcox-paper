# 【20260724】査読 §1.3 対応: 非同定 Q=2 を Q=1（=R, 同定）で再実行した mc_conf_conf724.rds から再描画。
# 【20260714】fig:confound を (B) 版 mc_conf_conf724.rds から再描画（PNG＋PDF）。フォント1.5倍（pointsize 12→18）。
# 【20260713B】(B) 版 mc_conf_conf724.rds（500reps, 全4設定）から再描画。
setwd(ifelse(dir.exists("D:/My Documents/Data/NMF-COXv3"),"D:/My Documents/Data/NMF-COXv3","~/NMF-COX"))
GAMMA0<-0.5
res <- readRDS("verification/mc_conf_rrr724.rds")   # 【20260724r1】RRR を第5推定量として追加（6本）
ord <- c("n300_rho.0_cA3_Q1","n300_rho.5_cA3_Q1","n300_rho.7_cA3_Q1","n600_rho.5_cA3_Q1"); ord<-ord[ord %in% names(res)]
est  <- c("cox_z","full","nmf","tvc","rrr","oracle")
elab <- c("(1) Cox (z-only)","(2) full-PH Cox","(3) NMF-COX (CF)","tvc-Cox (splines)","RRR (M=4,r=1)","oracle")
cat(sprintf("%-16s %3s %4s | %-7s %-7s %-7s %-7s %-7s | cov2/3/tvc\n","setting","n","rho","1","2","3","tvc","orac"))
summ<-list()
for(tg in ord){ M<-res[[tg]]$M; S<-res[[tg]]$setting
  bias<-sapply(est,function(e) mean(M[,e],na.rm=TRUE)-GAMMA0)
  cov2<-mean(abs(M[,"full"]-GAMMA0)<=1.96*M[,"se1"],na.rm=TRUE)
  cov3<-mean(abs(M[,"nmf"]-GAMMA0)<=1.96*M[,"se3"],na.rm=TRUE)
  covt<-mean(abs(M[,"tvc"]-GAMMA0)<=1.96*M[,"se_tvc"],na.rm=TRUE)
  summ[[tg]]<-list(S=S,M=M)
  cat(sprintf("%-16s %3d %4.1f | %+.3f %+.3f %+.3f %+.3f %+.3f | %.2f/%.2f/%.2f\n",
    tg,S$n,S$rho,bias["cox_z"],bias["full"],bias["nmf"],bias["tvc"],bias["oracle"],cov2,cov3,covt)) }
draw <- function(){
  par(mfrow=c(1,length(ord)),mar=c(4.3,0.6,3,0.6),oma=c(0,7.4,0,0.5))
  ecol<-c("gray45","black","black","gray30","gray20","gray25"); epch<-c(1,0,19,5,2,4)
  elwd<-c(1.6,1.6,2.2,1.7,1.7,1.8); ecex<-c(1.15,1.15,1.3,1.2,1.2,1.2); yb<-length(est):1
  for(ti in seq_along(ord)){ tg<-ord[ti]; M<-summ[[tg]]$M; S<-summ[[tg]]$S
    qmat<-sapply(est,function(e) quantile(M[,e],c(.025,.5,.975),na.rm=TRUE)); mn<-sapply(est,function(e)mean(M[,e],na.rm=TRUE))
    xr<-range(c(qmat,GAMMA0,0)); xr<-xr+c(-1,1)*0.04*diff(xr)
    plot(NA,xlim=xr,ylim=c(0.5,length(est)+0.5),yaxt="n",xlab=expression(hat(gamma)),ylab="",
         main=bquote(rho==.(S$rho)*","~n==.(S$n)))
    abline(v=GAMMA0,lty=1,col="gray55"); abline(v=0,lty=3,col="gray70")
    for(k in seq_along(est)){ y<-yb[k]; qq<-qmat[,k]
      segments(qq[1],y,qq[3],y,col=ecol[k],lwd=elwd[k])
      points(mn[k],y,pch=epch[k],col=ecol[k],lwd=elwd[k],cex=ecex[k],bg="white") }
    if(ti==1){ axis(2,at=yb,labels=elab,las=1,cex.axis=.92,tick=FALSE)
      mtext(expression(gamma[0]==0.5),side=3,at=GAMMA0,cex=.62,line=-0.1,col="gray40") } }
}
ptype<-if(isTRUE(capabilities("cairo")))"cairo" else "Xlib"
png("v3_fig/confound_forest.png",width=2300,height=980,res=175,type=ptype,pointsize=18); draw(); dev.off()
grDevices::cairo_pdf("v3_fig/confound_forest.pdf", width=2300/175, height=980/175, pointsize=18); draw(); dev.off()
cat("saved v3_fig/confound_forest.{png,pdf} (B)\n=== done ===\n")
