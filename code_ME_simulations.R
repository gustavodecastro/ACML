######################################
args   = commandArgs(TRUE)
TASKID = as.numeric(args[2])
njob   = TASKID
######################################

library(MASS)
library(glmmML)
library(survey)
require(rms)
require(Hmisc)
require(splines)
require(glmmML)
require(survey)
require(TwoPhaseReg)
library(sleev)
df = 500; Alpha=5; W=1

source("main_gen.R")
source("integrals_gen.R")

beta = c(1,.5,.5)
N = 4000; Nsim = 21
coef.Y = coef.wgt = coef.cml1 = coef.cml2 = coefs.cal = matrix(ncol=length(beta),nrow=Nsim)
ser.Y = ser.wgt = ser.cml1 = ser.cml2 = ser.cal = matrix(ncol=length(beta),nrow=Nsim)
coef.cml1b = coef.cml2b = matrix(ncol=length(beta),nrow=Nsim)
ser.cml1b = ser.cml2b = matrix(ncol=length(beta),nrow=Nsim)
error=0; nsimprint=10; Error = rep(0,Nsim); jj.track = rep(0,Nsim)
coefs_smle <- ses_smle <- coef.rak <- ser.rak <- coef.rakb <- ser.rakb <- matrix(NA, ncol=3, nrow=Nsim)

coef.wgt10 = coef.wgt1a = coef.wgt1b = coef.wgt1c = matrix(ncol=length(beta),nrow=Nsim)
coef.cml10 = coef.cml1a = coef.cml1b = coef.cml1c = matrix(ncol=length(beta),nrow=Nsim)
coef.cml20 = coef.cml2a = coef.cml2b = coef.cml2c = matrix(ncol=length(beta),nrow=Nsim)

ser.wgt10 = ser.wgt1a = ser.wgt1b = ser.wgt1c = matrix(ncol=length(beta),nrow=Nsim)
ser.cml10 = ser.cml1a = ser.cml1b = ser.cml1c = matrix(ncol=length(beta),nrow=Nsim)
ser.cml20 = ser.cml2a = ser.cml2b = ser.cml2c = matrix(ncol=length(beta),nrow=Nsim)

n = 400
rho <- 0
sig_errorprone <- sqrt(0.25)
sig_errorprone <- sqrt(1)
sig_errorprone <- sqrt(1.5)

N_SIEVE <- 30
degree <- 2

file_est <- paste0('coefs_SIM_Revision_error_prone_Nsim_Sieve_',N_SIEVE,'_sigmaEP_',sig_errorprone,'_rho_', rho*10, '_', N, n,'.txt')
file_var <- paste0('vars_SIM_Revision_error_prone_Nsim_Sieve_',N_SIEVE,'_sigmaEP_',sig_errorprone,'_rho_',  rho*10, '_', N, n, '.txt')

set.seed(2108 + njob)

for (ii in 1:Nsim) {
  x1 = rnorm(N)
  x3 = rnorm(N)
  x2 = rho*x1 + sqrt(1-rho)*x3
  x2_error = x2 + rnorm(N,0,sig_errorprone)
  
  e = rnorm(N)
  
  y = beta[1] + beta[2]*x1 + beta[3]*x2 + e
  
  cut1 = quantile(y, .25)
  cut2 = quantile(y, .75)
  Xcut1 = quantile(x1, .5)
  X2cut1= quantile(x2, .25)
  
  sy = sx1 = sx2 = rep(0,N)
  sy[y <= cut1] = 1
  sy[y >= cut2] = 2
  sx1[x1 <= Xcut1] = 1
  sx2[x2 <= X2cut1] = 1
  ind = 1:N; ymean = mean(y)
  
  R1y0 = sample(ind[sy==0],ceiling(n*0.4)); R1y1 = sample(ind[sy==1],n-2*ceiling(n*0.4)); R1y2 = sample(ind[sy==2],ceiling(n*0.4))
  R1 = rep(0,N); R1[R1y0] = 1; R1[R1y1] = 1; R1[R1y2] = 1
  
  counts = rep(1,N); use = (R1==1); always1 = rep(FALSE,N)
  data = data.frame(id = 1:N, R1=R1, y=y, dy=sy, x1=x1, x2=x2, x3=x2_error,
                    dx1=sx1, dx2=sx2, counts=counts,always1=always1,ycut=rep(c(cut1, cut2), length(y)/2))
  data$x2[data$R1==0] <- NA
  
  regformula = y ~ x1+x2
  desformula0 = R1 ~ dy
  desformula1 = R1 ~ dy + x3*x1*y
  desformula2 = R1 ~ dy + y*inf_naive
  desformula3 = R1 ~ dy + y*inf_naive_x1 + y*inf_naive_x2
  
  coef10 = f(regformula, desformula0, data, method="weighted")  
  coef1a = f(regformula, desformula1, data, method="weighted")
  bbeta <<- coef(lm(y ~ x1 + x2_error))
  inf_naive = (y - cbind(1,x1,x2_error) %*% bbeta[1:3])*x1
  data$inf_naive = inf_naive
  inf_naive_x1 = (y - cbind(1,x1,x2_error) %*% bbeta[1:3])*x1
  inf_naive_x2 = (y - cbind(1,x1,x2_error) %*% bbeta[1:3])*x2_error
  data$inf_naive_x1 = inf_naive_x1
  data$inf_naive_x2 = inf_naive_x2
  
  coef1b = f(regformula, desformula2, data, method="weighted")
  coef1c = f(regformula, desformula3, data, method="weighted")
  
  coef20 = f(regformula, desformula0, data=data, method="cml", niter=25, maxstepprop=.01)  
  coef2a = f(regformula, desformula1, data=data, method="cml", niter=25, maxstepprop=.01)
  coef2b = f(regformula, desformula2, data=data, method="cml", niter=25, maxstepprop=.01)
  coef2c = f(regformula, desformula3, data=data, method="cml", niter=25, maxstepprop=.01)
  
  coef30 = f(regformula, desformula0, Stilde = TRUE, data, method="cml", niter=25, maxstepprop=.01)
  coef3a = f(regformula, desformula1, Stilde = TRUE, data, method="cml", niter=25, maxstepprop=.01)
  coef3b = f(regformula, desformula2, Stilde = TRUE, data, method="cml", niter=25, maxstepprop=.01)
  coef3c = f(regformula, desformula3, Stilde = TRUE, data, method="cml", niter=25, maxstepprop=.01)
  
  coef.wgt10[ii,] = coef10$coef  
  coef.wgt1a[ii,] = coef1a$coef
  coef.wgt1b[ii,] = coef1b$coef
  coef.wgt1c[ii,] = coef1c$coef
  coef.cml10[ii,] = coef20$coef
  coef.cml1a[ii,] = coef2a$coef
  coef.cml1b[ii,] = coef2b$coef
  coef.cml1c[ii,] = coef2c$coef
  coef.cml20[ii,] = coef30$coef
  coef.cml2a[ii,] = coef3a$coef
  coef.cml2b[ii,] = coef3b$coef
  coef.cml2c[ii,] = coef3c$coef
  
  ### Raking
  data_rak <- data
  data_rak$id = 1:nrow(data_rak)
  inffun <- with(data_rak, c(y  - cbind(1, x1, x3) %*% t(t((bbeta))))*cbind(1, x1, x3))
  colnames(inffun) <- paste0('if', 1:3)
  data_rak <- cbind(data_rak, inffun)
  des2ph <- twophase(id = list(~1, ~1), strata = list(NULL, ~dy), data = data_rak, subset = R1==1)
  des2ph_cal <- survey::calibrate(des2ph, phase = 2, calfun = 'raking', formula = ~ if1 + if2 + if3)
  rak2ph <- svyglm(regformula, design = des2ph_cal)
  
  des2ph_calb <- survey::calibrate(des2ph, phase = 2, calfun = 'raking', formula = ~ if1*dy + if2*dy + if3*dy + x3*x1*y)
  rak2phb <- svyglm(regformula, design = des2ph_calb)
  
  coef.rak[ii,] = summary(rak2ph)$coef[,1]; ser.rak[ii,] = summary(rak2ph)$coef[,2]
  coef.rakb[ii,] = summary(rak2phb)$coef[,1]; ser.rakb[ii,] = summary(rak2phb)$coef[,2]
  
  ###
  N_SIEVE <- N_SIEVE
  Bspline <- bs(data$x3, df=N_SIEVE, degree=degree, Boundary.knots=range(data$x3), intercept=TRUE)
  colnames(Bspline) <- paste("bs", 1:N_SIEVE, sep="")
  dat = data.frame(Y=data$y, X=data$x2, Z=data$x1, X3=data$x3, R1=data$R1, Bspline)
  dat[dat$R1==0,"X"] = NA
  res_linear = sleev::linear2ph(Y="Y", X="X", Y_unval="Y", X_unval="X3",
                                Z="Z", Bspline=colnames(Bspline), data=dat)
  
  coefs_smle[ii,] <- res_linear$coefficients[c(1,3,2),1]
  ses_smle[ii,] <- res_linear$coefficients[c(1,3,2),2]
  
  ddest <- matrix(c(coef.wgt10[ii,], 
                    coef.wgt1a[ii,], 
                    coef.wgt1b[ii,], 
                    coef.wgt1c[ii,],
                    coef.cml10[ii,], 
                    coef.cml1a[ii,], 
                    coef.cml1b[ii,], 
                    coef.cml1c[ii,],
                    coef.cml20[ii,], 
                    coef.cml2a[ii,], 
                    coef.cml2b[ii,], 
                    coef.cml2c[ii,],
                    coef.rak[ii,], 
                    coef.rakb[ii,],
                    coefs_smle[ii,]), nrow=1)
  
  if (ii == 1) {
    save2 <- paste('N = ', N, 'n1 = ', n, 'beta = ', paste(beta, collapse=" "), ' njob = ', njob)
    if (!file.exists(file_est)) {
      file.create(file_est)
      write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
    }
  }
  
  write.table(ddest, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)
}
