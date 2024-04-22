
# @File    :   3_2-Four-countries-subtype-CCM.R
# @Time    :   2024/4/5 16:12:42
# @Author  :   DuLab
# @Desc    :   Causal relationship between A/H3N2 and A/H1N1, and between A/H1N1 and B at the country level.



library(rEDM)
filename = './data/'
datas = read.csv(paste(filename,'1-Thailand.csv',sep=''))


h1 <- datas$AH1[13:95]
h3 <- datas$AH3[1:83]    ##AH1 xmap AH3
datas <- data.frame(AH3=h3,AH1=h1)

#CCM
ts<-datas$AH1     ##Select the best dimension for the result AH1
lib<-c(1,length(ts))
pred<-c(1,length(ts))
simplex_output<-simplex(ts,lib,pred,silent = TRUE)
smap_output<-s_map(ts,lib,pred,E=6,silent = TRUE)

n <- NROW(datas)
vars <- c('AH3','AH1')
ccm_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars))
for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    out_temp <- ccm(datas, E =6, lib_column = ccm_from, target_column = ccm_to, 
                    lib_sizes = n, replace = FALSE, silent = TRUE)
    ccm_matrix[ccm_from, ccm_to] <- out_temp$rho
  }
}

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,  vars))
for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    cf_temp <- ccf(datas[, ccm_from], datas[, ccm_to], type = "correlation", 
                   lag.max = n, plot = FALSE)$acf
    corr_matrix[ccm_from, ccm_to] <- max(abs(cf_temp))
  }
}


#thrips_xmap_maxT
paramecium_xmap_didinium <- ccm(datas, E =6, random_libs = TRUE, lib_column = "AH1", 
                                target_column = "AH3", lib_sizes = seq(5, 83, by = 5), num_samples = 300,silent = TRUE)

para_map_didi = ccm_means(paramecium_xmap_didinium)

par(mar = c(4, 4, 2, 2), mgp = c(2.5, 1, 0),cex.axis=1.25,cex.lab=1.25)  # set up margins for plotting

plot(para_map_didi$lib_size,pmax(0, para_map_didi$rho), type = "l", col = "blue", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), lwd = 4, cex.lab=1.3)

title("Thailand",cex.main=1.3)
legend(x = "topleft", col = c("blue"), lwd = 4, legend = c("AH1 xmap AH3 , P = 0.006"),
       inset = 0.02, bty = "n", cex = 1.2)
abline(h = 0, lty = 4)


#Hypothesis testing, testing for the cause AH3
num_surr <- 1000
surr_A3 <- make_surrogate_data(datas$AH3, method = "seasonal", 
                               T_period = 12, num_surr = num_surr)
rho_surr <- data.frame(A3 = numeric(num_surr))

for (i in 1:num_surr) {
  rho_surr$A3[i] <- ccm(cbind(datas$AH1, surr_A3[,  i]),
                        E =6, lib_column = 1, target_column = 2, lib_sizes = NROW(datas), 
                        replace = FALSE)$rho
}
(sum(ccm_matrix['AH1','AH3'] < rho_surr$A3) + 1) / 
  (length(rho_surr$A3) + 1)




#ZHIHOU
library(rEDM)
filename = './data/'
datas = read.csv(paste(filename,'1-Thailand.csv',sep=''))


h1 <- datas$AH1[1:83]      ##AH3 xmap AH1
h3 <- datas$AH3[13:95]   
datas <- data.frame(AH3=h3,AH1=h1)

#CCM
ts<-datas$AH3     ##Select the best dimension for the result AH3
lib<-c(1,length(ts))
pred<-c(1,length(ts))
simplex_output<-simplex(ts,lib,pred,silent = TRUE)
smap_output<-s_map(ts,lib,pred,E=3,silent = TRUE)


n <- NROW(datas)
vars <- c('AH3','AH1')
ccm_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars, vars))
for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    out_temp <- ccm(datas, E =3, lib_column = ccm_from, target_column = ccm_to, 
                    lib_sizes = n, replace = FALSE, silent = TRUE)
    ccm_matrix[ccm_from, ccm_to] <- out_temp$rho
  }
}

corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,  vars))
for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    cf_temp <- ccf(datas[, ccm_from], datas[, ccm_to], type = "correlation", 
                   lag.max = n, plot = FALSE)$acf
    corr_matrix[ccm_from, ccm_to] <- max(abs(cf_temp))
  }
}


#thrips_xmap_maxT
paramecium_xmap_didinium <- ccm(datas, E =3, random_libs = TRUE, lib_column = "AH3", 
                                target_column = "AH1", lib_sizes = seq(5, 83, by = 5), num_samples = 300,silent = TRUE)

para_map_didi = ccm_means(paramecium_xmap_didinium)

par(mar = c(4, 4, 2, 2), mgp = c(2.5, 1, 0),cex.axis=1.25,cex.lab=1.25)  # set up margins for plotting

plot(para_map_didi$lib_size,pmax(0, para_map_didi$rho), type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), lwd = 4, cex.lab=1.3)

title("Thailand",cex.main=1.3)
legend(x = "topleft", col = c("red"), lwd = 4, legend = c("AH3 xmap AH1, P = 0.919"),
       inset = 0.02, bty = "n", cex = 1.2)
abline(h = 0, lty = 4)


#Hypothesis testing, testing for the cause AH1
num_surr <- 1000
surr_A1 <- make_surrogate_data(datas$AH1, method = "seasonal", 
                               T_period = 12, num_surr = num_surr)
rho_surr <- data.frame(A1 = numeric(num_surr))

for (i in 1:num_surr) {
  rho_surr$A1[i] <- ccm(cbind(datas$AH3, surr_A1[,  i]),
                        E = 3, lib_column = 1, target_column = 2, lib_sizes = NROW(datas), 
                        replace = FALSE)$rho
}
(sum(ccm_matrix['AH3','AH1'] < rho_surr$A1) + 1) / 
  (length(rho_surr$A1) + 1)
