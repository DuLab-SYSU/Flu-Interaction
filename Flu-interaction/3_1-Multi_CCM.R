
# @File    :   3_1-Multi_CCM.R
# @Time    :   2024/4/5 16:12:42
# @Author  :   DuLab
# @Desc    :   Causal analysis of influenza A and influenza B on a global level.



library("multispatialCCM")
library(rEDM)
filename = './data/'
files = read.csv(paste(filename,'3_1-Multi_CCM.csv',sep=''))
Accm<-files$AH3
Bccm<-files$B
maxE<-5
Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("AH3", "B")
for(E in 2:maxE) {Emat[E-1,"AH3"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
}
matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,xlab="E", ylab="rho", lwd=2)
legend("bottomleft", c("AH3", "B"), lty=1:2, col=1:2, lwd=2, bty="n")
E_A<-5
E_B<-2

signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,predsteplist=1:9)
signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,predsteplist=1:9)

#Does "A cause B"
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=1000)
#Does "B cause A" 
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=1000)
#test
(CCM_significance_test<-ccmtest(CCM_boot_A,CCM_boot_B))
#plot result
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))

#Plot “A causes B”
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, 
     lwd=4,xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),xlab="L", ylab="rho")
matlines(CCM_boot_A$Lobs,cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
                               CCM_boot_A$rho+CCM_boot_A$sdevrho),lty=4, col=1)

#Plot "B causes A"
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=4)
matlines(CCM_boot_B$Lobs,cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
                               CCM_boot_B$rho+CCM_boot_B$sdevrho),lty=4, col=2)
legend("topleft",c("AH3 causes B , P = 0.045","B causes AH3 , P = 0.135"),lty=c(1,2), col=c(1,2), lwd=2, bty="n")
