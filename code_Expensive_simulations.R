######################################
args   = commandArgs(TRUE)
TASKID = as.numeric(args[2])
njob   = TASKID
######################################

require(rms)
require(Hmisc)
require(splines)
require(glmmML)
require(survey)
require(mvtnorm)
require(TwoPhaseReg)
library(sleev)
library(mice)
df = 500; Alpha=5; W=1

source("main.R")
source("integrals.R")

beta = c(1,.5,1)
N = 10000; Nsim = 11
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

n = 1000
rho <- 0.5

file_est <- paste0('coefs_expensive_REV_SIM_Nsim_',Nsim,'_rho_', rho*10, '_', N, n,'.txt')
file_var <- paste0('vars_expensive_REV_SIM_Nsim_',Nsim,'_rho_',  rho*10, '_', N, n, '.txt')

set.seed(2108 + njob)

for (ii in 1:Nsim) {
  x1 = rnorm(N)
  x3 = rnorm(N)
  x2 = rho*x1 + sqrt(1-rho)*x3
  
  e = rnorm(N)
  
  y = beta[1] + beta[2]*x1 + beta[3]*x2 + e
  
  cut1 = quantile(y, .33)
  cut2 = quantile(y, .67)
  
  sy = sx1 = sx2 = rep(0,N)
  sy[y <= cut1] = 1
  sy[y >= cut2] = 2
  sx1[x1 <= quantile(x1, .50)] = 1
  sx2[x2 <= quantile(x2, .50)] = 1
  ind = 1:N
  
  R1y0 = sample(ind[sy==0], ceiling(n*0.4)); R1y1 = sample(ind[sy==1], n-2*ceiling(n*0.4)); R1y2 = sample(ind[sy==2], ceiling(n*0.4))
  R1 = rep(0,N); R1[R1y0] = 1; R1[R1y1] = 1; R1[R1y2] = 1
  
  counts = rep(1,N); use = (R1==1); always1 = rep(FALSE,N)
  data = data.frame(id = 1:N, R1=R1, y=y, dy=sy, x1=x1, x2=x2,
                    dx1=sx1, dx2=sx2, counts=counts,always1=always1,ycut=rep(c(cut1, cut2), length(y)/2))
  data$x2[data$R1==0] <- NA
  
  regformula = y ~ x1+x2
  desformula0 = R1 ~ dy
  desformula1 = R1 ~ dy + dx2*x1*y
  desformula2 = R1 ~ dy + dx2*x1*y + y*inf_naive_x1
  desformula3 = R1 ~ dy + dx2*x1*y + y*inf_naive_x1 + y*inf_naive_x2
  
  coef10 = f(regformula, desformula0, data, method="weighted")  
  coef1a = f(regformula, desformula1, data, method="weighted")
  
  Nimpute <- 1
  mod_imp <- lm(x2 ~ y*x1+dx2, data=data, subset=R1==1)
  sig_imp <- sqrt(summary(mod_imp)$sigma^2*rchisq(Nimpute, summary(mod_imp)$df[2], ncp=0)/summary(mod_imp)$df[2])
  desmatrix <- rmvnorm(Nimpute, mean=mod_imp$coeff, sigma=summary(mod_imp)$sigma^2*summary(mod_imp)$cov)
  x2imp   <- as.vector(desmatrix[1,] %*% t(cbind(1, data$y, data$x1, data$dx2, data$y*data$x1)) + rnorm(nrow(data), 0, sig_imp[1]))
  data$x3 <- x2imp
  
  bbeta <<- coef(lm(y ~ x1 + x3, data = data))
  inf_naive_x0 = (y - cbind(1,data$x1,data$x3) %*% bbeta[1:3])
  inf_naive_x1 = (y - cbind(1,data$x1,data$x3) %*% bbeta[1:3])*data$x1
  inf_naive_x2 = (y - cbind(1,data$x1,data$x3) %*% bbeta[1:3])*data$x3
  data$inf_naive = inf_naive_x1
  data$inf_naive_x0 = inf_naive_x0
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
  
  coef.wgt10[ii,] = coef10$coef; ser.wgt10[ii,] = coef10$ser  
  coef.wgt1a[ii,] = coef1a$coef; ser.wgt1a[ii,] = coef1a$ser
  coef.wgt1b[ii,] = coef1b$coef; ser.wgt1b[ii,] = coef1b$ser
  coef.wgt1c[ii,] = coef1c$coef; ser.wgt1c[ii,] = coef1c$ser
  coef.cml10[ii,] = coef20$coef; ser.cml10[ii,] = coef20$ser
  coef.cml1a[ii,] = coef2a$coef; ser.cml1a[ii,] = coef2a$ser
  coef.cml1b[ii,] = coef2b$coef; ser.cml1b[ii,] = coef2b$ser
  coef.cml1c[ii,] = coef2c$coef; ser.cml1c[ii,] = coef2c$ser
  coef.cml20[ii,] = coef30$coef; ser.cml20[ii,] = coef30$ser
  coef.cml2a[ii,] = coef3a$coef; ser.cml2a[ii,] = coef3a$ser
  coef.cml2b[ii,] = coef3b$coef; ser.cml2b[ii,] = coef3b$ser
  coef.cml2c[ii,] = coef3c$coef; ser.cml2c[ii,] = coef3c$ser 
  
  #  print(ii)
  
  #  if (ii > 1){
  #    cat('\nWGT-1', round(apply(coef.wgt1a[1:ii,], 2, var)*1000, 3))
  #    cat('\nWGT-2', round(apply(coef.wgt1b[1:ii,], 2, var)*1000, 3))
  #    cat('\nWGT-3', round(apply(coef.wgt1c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nCML-1', round(apply(coef.cml1a[1:ii,], 2, var)*1000, 3))
  #    cat('\nCML-2', round(apply(coef.cml1b[1:ii,], 2, var)*1000, 3))
  #    cat('\nCML-3', round(apply(coef.cml1c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nACML-1', round(apply(coef.cml2a[1:ii,], 2, var)*1000, 3))
  #    cat('\nACML-2', round(apply(coef.cml2b[1:ii,], 2, var)*1000, 3))
  #    cat('\nACML-3', round(apply(coef.cml2c[1:ii,], 2, var)*1000, 3))
  #  }
  #}
  
  #coef3 = f(regformula, desformula, data, method="cml", Stilde=TRUE, niter=50, maxstepprop=.1)
  #coef.cml2[ii,] = coef3$coef; ser.cml2[ii,] = coef3$ser
  
  ### Raking
  data_rak <- data
  data_rak$id = 1:nrow(data_rak)
  des2ph <- twophase(id = list(~1, ~1), strata = list(NULL, ~dy), data = data_rak, subset = R1==1)
  des2ph_cal <- survey::calibrate(des2ph, phase = 2, calfun = 'raking', formula = ~ y*inf_naive_x0 + y*inf_naive_x1 + y*inf_naive_x2)
  rak2ph <- svyglm(regformula, design = des2ph_cal)
  
  des2ph_calb <- survey::calibrate(des2ph, phase = 2, calfun = 'raking', formula = ~ y*inf_naive_x0 + y*inf_naive_x1 + y*inf_naive_x2 + dx2*x1*y)
  rak2phb <- svyglm(regformula, design = des2ph_calb)
  
  coef.rak[ii,] = summary(rak2ph)$coef[,1]; ser.rak[ii,] = summary(rak2ph)$coef[,2]
  coef.rakb[ii,] = summary(rak2phb)$coef[,1]; ser.rakb[ii,] = summary(rak2phb)$coef[,2]
  ###
  
  #  print(ii)
  #  
  #  if (ii > 1){
  #    
  #    cat('\nWGT-1', round(apply(coef.wgt1a[1:ii,], 2, var)*1000, 3))
  #    cat('\nWGT-2', round(apply(coef.wgt1b[1:ii,], 2, var)*1000, 3))
  #    cat('\nWGT-3', round(apply(coef.wgt1c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nCML-1', round(apply(coef.cml1a[1:ii,], 2, var)*1000, 3))
  #    cat('\nCML-2', round(apply(coef.cml1b[1:ii,], 2, var)*1000, 3))
  #    cat('\nCML-3', round(apply(coef.cml1c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nACML-1', round(apply(coef.cml2a[1:ii,], 2, var)*1000, 3))
  #    cat('\nACML-2', round(apply(coef.cml2b[1:ii,], 2, var)*1000, 3))
  #    cat('\nACML-3', round(apply(coef.cml2c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nRak-1', round(apply(coef.rak[1:ii,], 2, var)*1000, 3))
  #    cat('\nRak-2', round(apply(coef.rakb[1:ii,], 2, var)*1000, 3))
  #  }
  #}
  
  
  
  indX2_0 <- which(sx2 == 0)
  nX2_0   <- length(indX2_0)
  indX2_1 <- which(sx2 == 1)
  nX2_1   <- length(indX2_1)
  
  N_SIEVE <- 10
  Bspline_0 <- bs(data$x1[indX2_0], df=N_SIEVE, degree=2, Boundary.knots=range(data$x1[indX2_0]), intercept=TRUE)
  Bspline_1 <- bs(data$x1[indX2_1], df=N_SIEVE, degree=2, Boundary.knots=range(data$x1[indX2_1]), intercept=TRUE)
  Bspline <- matrix(0, N, 2*N_SIEVE)
  Bspline[indX2_0,1:N_SIEVE]               <- Bspline_0
  Bspline[indX2_1,(N_SIEVE+1):(2*N_SIEVE)] <- Bspline_1
  colnames(Bspline) <- paste("bs", 1:(2*N_SIEVE), sep="")
  
  dat = data.frame(Y=data$y, X=data$x2, Z=data$x1, R1=data$R1, Bspline)
  dat[dat$R1==0,"X"] = NA
  res_linear = smle(Y="Y", X="X", Z="Z", Bspline_Z=colnames(Bspline), data=dat)
  coefs_smle[ii,] <- res_linear$coefficients[c(1,3,2),1]
  ses_smle[ii,] <- res_linear$coefficients[c(1,3,2),2]
  
  
  #  if (ii > 1){
  #    
  #    cat('\nWGT-1', round(apply(coef.wgt1a[1:ii,], 2, var)*1000, 3))
  #    cat('\nWGT-2', round(apply(coef.wgt1b[1:ii,], 2, var)*1000, 3))
  #    cat('\nWGT-3', round(apply(coef.wgt1c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nCML-1', round(apply(coef.cml1a[1:ii,], 2, var)*1000, 3))
  #    cat('\nCML-2', round(apply(coef.cml1b[1:ii,], 2, var)*1000, 3))
  #    cat('\nCML-3', round(apply(coef.cml1c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nACML-1', round(apply(coef.cml2a[1:ii,], 2, var)*1000, 3))
  #    cat('\nACML-2', round(apply(coef.cml2b[1:ii,], 2, var)*1000, 3))
  #    cat('\nACML-3', round(apply(coef.cml2c[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nRak-1', round(apply(coef.rak[1:ii,], 2, var)*1000, 3))
  #    cat('\nRak-2', round(apply(coef.rakb[1:ii,], 2, var)*1000, 3))
  #    
  #    cat('\nSPMLE', round(apply(coefs_smle[1:ii,], 2, var)*1000, 3))
  #  }
  #}
  
  
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
  ddvar <- matrix(c(ser.wgt10[ii,], 
                    ser.wgt1a[ii,], 
                    ser.wgt1b[ii,], 
                    ser.wgt1c[ii,],
                    ser.cml10[ii,], 
                    ser.cml1a[ii,], 
                    ser.cml1b[ii,], 
                    ser.cml1c[ii,],
                    ser.cml20[ii,], 
                    ser.cml2a[ii,], 
                    ser.cml2b[ii,], 
                    ser.cml2c[ii,],
                    ser.rak[ii,], 
                    ser.rakb[ii,],
                    ses_smle[ii,]), nrow=1)
  
  if (ii == 1) {
    save2 <- paste('N = ', N, 'n1 = ', n, 'beta = ', paste(beta, collapse=" "), ' njob = ', njob)
    if (!file.exists(file_est)) {
      file.create(file_est)
      write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
      file.create(file_var)
      write.table(save2, file = file_var,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
    }
  }
  
  write.table(ddest, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)
  write.table(ddvar, file = file_var, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)  
}
