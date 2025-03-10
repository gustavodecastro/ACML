#library(Hmisc)
#source("func.bin.R")
get2 = function(n,yu,v1,v2){
  sample(ind[R1==1 & y==yu & x1d==v1 & x2d==v2],n)}

N = 2000; beta = c(-3.3,1); alpha <- c(-3.5,3.5)
rho12 = 0.9

nmeth = 5; meth.names=c("Weighted simple", "Cond. simple",
                        "Stilde simple", "Stilde complete", "error_prone")

Nsims <- 1000
nsimprint=10
thetas.off = thetas.wtd = ses.off = ses.wtd = alpha2s = alpha3s = list(0)
for (i in 1:nmeth) thetas.off[[i]] = thetas.wtd[[i]] = ses.off[[i]] = ses.wtd[[i]] = matrix(0,Nsims,length(beta))
nconvfails.old = nconvfails.new = out.off = 0

for (ii in 1:Nsims){
  x1 = rnorm(N)
  x2 = rho12*x1 + sqrt(1-rho12^2)*rnorm(N)
  ptrue = exp(beta[1] + beta[2]*x2)/(1+exp(beta[1] + beta[2]*x2))
  y = rbinom(N, 1, ptrue)
  
  ## 2 variables
  psampl = exp(alpha[1] + alpha[2]*y)/(1+exp(alpha[1] + alpha[2]*y))
  R1 = rbinom(N, 1, psampl)
  
  x2[R1==0] = NA
  
  counts=rep(1,N)
  always1 = rep(FALSE, N)
  
  datcsv = data.frame(R1=R1,y=y,x1=x1,x2=x2,counts=counts,always1=always1,always2=always1)
  
  regformula = y~x2
  desnformula1 = R1~y
  
  out.off.1 = func(regformula, R1~y*x1, desnformula2=NULL, desnformula3=NULL, data=datcsv, method="weighted")
  thetas.off[[1]][ii,] = out.off.1$coefficients
  ses.off[[1]][ii,] = sqrt(diag(out.off.1$cov))[1:length(beta)]
  
  out.off.2 = func(regformula, R1~y*x1, desnformula2=NULL, desnformula3=NULL, data=datcsv)
  thetas.off[[2]][ii,] = out.off.2$coefficients
  ses.off[[2]][ii,] = sqrt(diag(out.off.2$cov))[1:length(beta)]
  
  out.off.3 = func(regformula, R1~y, desnformula2=NULL, desnformula3=NULL, data=datcsv, addStilde=TRUE, niter=50)
  thetas.off[[3]][ii,] = out.off.3$coefficients
  ses.off[[3]][ii,] = sqrt(diag(out.off.3$cov))[1:length(beta)]
  
  out.off.4 = func(regformula, R1~y+x1, desnformula2=NULL, desnformula3=NULL, data=datcsv, addStilde=TRUE, niter=50)
  thetas.off[[4]][ii,] = out.off.4$coefficients
  ses.off[[4]][ii,] = sqrt(diag(out.off.4$cov))[1:length(beta)]
  
  mod_all <- glm(y ~ x1, family = binomial)
  thetas.off[[5]][ii,] = mod_all$coefficients
  
  if (ii %% 100 == 0) {
    cat("\n\nNiter=",ii)
    
    cat('\n wgt:', colMeans(thetas.off[[1]][1:ii,]))
    cat('\n cml:', colMeans(thetas.off[[2]][1:ii,]))
    cat('\n cmlS:', colMeans(thetas.off[[3]][1:ii,]))
    cat('\n cmlS:', colMeans(thetas.off[[4]][1:ii,]))
    cat('\n ErrorProne:', colMeans(thetas.off[[5]][1:ii,]))
    
    cat('\n Var_wgt:', apply(thetas.off[[1]][1:ii,], 2, sd))
    cat('\n Var_cml:', apply(thetas.off[[2]][1:ii,], 2, sd))
    cat('\n Var_cmlS:', apply(thetas.off[[3]][1:ii,], 2, sd))
    cat('\n Var_cmlS:', apply(thetas.off[[4]][1:ii,], 2, sd))
    cat('\n ErrorProne:', apply(thetas.off[[5]][1:ii,], 2, sd))
  }
}
