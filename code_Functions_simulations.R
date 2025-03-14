f = function(regformula,desformula,data,method,niter=20,maxstepprop=0.1,tol=1e-6,Stilde=FALSE) {
  
  M1 = model.apply(regformula,data,counts,always1)
  y = M1$y; x=M1$X; y.aux=y; x.aux=x
  
  M2 = model.apply(desformula,data,counts,always1)
  R1 = M2$y; z=M2$X; w=z; use=(R1==1)
  fit = glm(R1 ~ w-1, family=binomial)
  alpha1 = as.matrix(as.vector(coefficients(fit))); pi1 = plogis(w%*%alpha1)
  
  #############################################################################
  if (method == "weighted") { ############################### Weighted method
    coef = lm(y ~ x-1, weights=1/pi1, subset=use)$coef
    sds = (summary(lm(y ~ x-1, weights=1/pi1, subset=use))$sigma)
    sds = sd(y - fitted(lm(y ~ x-1, weights=1/pi1, subset=use)))
    #sds = sd(y - fitted(lm(y ~ x-1, subset=use)))
    x = x[use,]; y = y[use]; pi11 = pi1[use]
    
    I00 = t(x/sds^2) %*% diag(as.vector(1/pi11)) %*% x
    I01 = t(x/sds^2) %*% diag(as.vector((y-x%*%coef)*(1-pi11)*(1/pi11))) %*% w[R1==1,]
    w.aux = w*as.vector(sqrt(pi1*(1-pi1)))
    I11 = t(w.aux) %*% w.aux
    I10 = matrix(rep(0,ncol(I00)*nrow(I11)),nrow=nrow(I11),ncol=ncol(I00))
    I22 = sum(as.vector(2/(sds^2*pi11) - 3*(y-x%*%coef)^2/(sds^(4)*pi11)))
    I20 = -t(as.vector(as.vector(2*(y-x%*%coef)/pi11)/sds^(3))) %*% x
    I21 = -t(as.vector(-2*(1-pi11)/(sds*pi11) + ((y-x%*%coef)^2*(1-pi11)/pi11)/sds^(3))) %*% w[R1==1,]
    I = rbind(cbind(I00,I01,t(I20)),
              cbind(I10,I11,t(I21)),
              cbind(I20,I21,I22))
    S0 = x*as.vector((y-x%*%coef)/(sds^2*pi11))
    S1 = w*as.vector(R1-pi1)
    S2 = as.vector(-2/(sds*pi11) + (y-x%*%coef)^2/(sds^(3)*pi11))
    S2 = as.vector(-2/(sds) + (y-x%*%coef)^2/(sds^(3)))
    S = t(cbind(t(colSums(S0)), t(colSums(S1)), sum(S2)))
    
    C00 = empcov(S0,nreps=counts[use])$SS
    C22 = empcov(S2,nreps=counts[use])$SS
    C01 = empcov(x=S0,y=S1[use,],nreps=counts[use])$SS
    C21 = empcov(x=S2,y=S1[use,],nreps=counts[use])$SS
    C02 = empcov(x=S0,y=S2,nreps=counts[use])$SS
    C11 = empcov(S1,nreps=counts)$SS
    C  = rbind(cbind(C00,C01,C02),
               cbind(t(C01),C11,t(C21)),
               cbind(t(C02),C21,C22))
    ser = sqrt(diag(solve(I,t(solve(I,C))))[1:ncol(x)])
    
    list(coef=coef,ser=ser)
    
    #################################################################
  } else if (method == "cml") { ############################# cml
    
    conv <- FALSE
    coef = lm(y ~ x-1, weights=1/pi1, subset=use)$coef; coef.fixed = coef
    
    sig <- c(sqrt(var(y[use] - x[use,]%*%coef)))
    data$sig <- sig
    nxi = length(xi)
    results = tempfunction(data,x,xi,nxi,coef,desformula)
    ztemp = results$ztemp; xtemp = results$xtemp

    pi11 = pi1[use]
    phi = c(coef, sig, alpha1)
    #phi = c(coef, alpha1)

    for (jj in 1:niter) {
      if (jj>1) {
        coef = phi[1:ncol(x)]
        sig  <- phi[ncol(x)+1]
        data$sig <- sig
        #data$sig <- c(sqrt(var(y[use] - x[use,]%*%coef)))
        alpha1 = phi[-(1:(ncol(x)+1))]; pi1 = plogis(z%*%alpha1)
        #alpha1 = phi[-(1:(ncol(x)))]; pi1 = plogis(z%*%alpha1)
        results = tempfunction(data,x,xi,nxi,coef,desformula)
        ztemp = results$ztemp; xtemp = results$xtemp
      }
      
      S0 <- S0func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use)
      S1 <- S2func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use)
      S2 = z*as.vector(R1-pi1)
      S = t(cbind(t(colSums(S0)), t(colSums(S1)), t(colSums(S2))))
      #S = t(cbind(t(colSums(S0)), t(colSums(S2))))
      
      I00 <- -hessian_betabeta(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use)
      I01 <- -t(hessian_betasigma(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use))
      I02 <- -t(hessian_betaalpha(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use))
      I10 <- t(I01)
      I11 <- -hessian_sigmasigma(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use)
      I12 <- -t(hessian_sigmaalpha(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use))

      y.aux <- z*as.vector(sqrt(pi1*(1-pi1)))
      I22 <- t(y.aux) %*% y.aux
      I20 <- matrix(c(rep(0,nrow(I22)*ncol(I00))), ncol=ncol(I00), nrow=nrow(I22))
      I21 <- matrix(c(rep(0,nrow(I22)*ncol(I11))), ncol=ncol(I11), nrow=nrow(I22))
      
      if (Stilde) {
        Itilde20 = t(I02)
        Itilde21 = t(I12)
        Itilde22 <- -hessian_alphaalpha(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use)
        I22 = I22 - Itilde22
        I20 = - Itilde20
        I21 = - Itilde21
        
        S2tilde <- Stilde2func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use)
        
        S = t(cbind(t(colSums(S0)), t(colSums(S1)), t(colSums(S2) - colSums(S2tilde))))
        #S = t(cbind(t(colSums(S0)), t(colSums(S2) - colSums(S2tilde))))
      }
      
      II = rbind(cbind(I00,I01,I02),
                 cbind(I10,I11,I12),
                 cbind(I20,I21,I22))
      
      #II = rbind(cbind(I00,I02),
      #           cbind(I20,I22))
      
      Step = 0
      if (niter >1) {
        Step = solve(II,S)
      }
      bigstep = max(abs(Step))
      maxstep = maxstepprop*max(abs(phi))
      if (bigstep > maxstep) {
        Step = Step * maxstep/bigstep
      }
        
      phi2 = phi + Step
      if (niter>1) if (max(abs(phi-phi2)[1:ncol(x)]) < tol*max(abs(phi)[1:ncol(x)])) {
        conv = TRUE
        break
      }
      phi = phi2
    } ##End loop

    coef = phi[1:ncol(x)]
    if (niter>1 & !conv) {
      error = error + 1
      Error[ii] = ii
      print(paste("Erros = ",error))
    }
    
    list(coef=coef,conv=conv)
    
  } else print(paste("Method not available"))
  
} ##End function


####################
## TEMP FUNCTIONS ##
####################
tempfunction = function(data,x,xi,nxi,coef,desformula) {
  yrep = ytemp = R1temp = dxtemp = dytemp = dx2temp = rep(0,nxi*sum(data$R1))
  x3temp = rep(0,nxi*sum(data$R1))
  inf_naive_temp = inf_naive_x1_temp = inf_naive_x2_temp = rep(0,nxi*sum(data$R1))
  xtemp = matrix(ncol=ncol(x), nrow=nxi*sum(data$R1))
  rows = as.vector(which(data$R1==1))
  for (cont in 1:sum(R1)) {
    index = rows[cont]
    ytemp[(1:nxi)+(cont-1)*nxi]   = rep((x%*%coef)[index],nxi)
    yrep[(1:nxi)+(cont-1)*nxi]    = rep(rep(data$y[index]),nxi)
    xtemp[(1:nxi)+(cont-1)*nxi,]  = matrix(c(rep(x[index,],nxi)), ncol=ncol(x), byrow=TRUE)
    R1temp[(1:nxi)+(cont-1)*nxi]  = rep(data$R1[index],nxi)
    dxtemp[(1:nxi)+(cont-1)*nxi]  = rep(data$dx1[index],nxi)
    dx2temp[(1:nxi)+(cont-1)*nxi] = rep(data$dx2[index],nxi)
    x3temp[(1:nxi)+(cont-1)*nxi] = rep(data$x3[index],nxi)
  }
  ytemp = ytemp + rep(sqrt(2)*data$sig[1]*xi,sum(data$R1))
  #dytemp = 1*(ytemp <= data$ycut[1])
  #dytemp = 2*(ytemp >= data$ycut[2])
  dytemp <- ifelse(ytemp <= data$ycut[1], 1,
                   ifelse(ytemp >= data$ycut[2], 2, 0))
  #dytemp <- rep(0, length(ytemp))
  #dytemp = 1*(ytemp <= data$ycut[1] | ytemp >= data$ycut[2])
  #bbeta <<- coef(lm(ytemp ~ xtemp[,2] + x3temp))
  inffun1 = (ytemp - cbind(xtemp[,c(1:2)], x3temp) %*% bbeta[1:3])*xtemp[,2]
  inffun2 = (ytemp - cbind(xtemp[,c(1:2)], x3temp) %*% bbeta[1:3])*x3temp
  countstemp = rep(1,length(ytemp)); always1temp = rep(FALSE,length(ytemp))
  datatemp = data.frame(R1=R1temp, y=ytemp, x1=xtemp[,2], x2=xtemp[,3], x3=x3temp,
                        inf_naive=inffun1, inf_naive_x1=inffun1, inf_naive_x2=inffun2,
                        dx1=dxtemp, dy=dytemp, dx2=dx2temp, countstemp=countstemp, always1temp=always1temp)
  M4 = model.apply(desformula,datatemp,countstemp,always1temp)
  if (attr(M4$Terms,"response") != 1) stop("regformula must have a response variable")
  R1temp = M4$y; ztemp=M4$X
  list(datatemp=datatemp,xtemp=xtemp,ytemp=ytemp,ztemp=ztemp,
       countstemp=countstemp,always1temp=always1temp,dxtemp=dxtemp,yrep=yrep)
}

###############
## FUNCTIONS ##
###############
pi.func    = function(z,alpha1) {
  exp(z%*%alpha1)/(1+exp(z%*%alpha1))
}

###############################################################################
## Gauss-Hermite quadrature to calculate the integrals (40 points were used) ##
###############################################################################
ghq20 = ghq(40,FALSE)
wi  = ghq20$weights
xi  = ghq20$zeros

###############
## INTEGRAIS ##
###############
Integral = function(alpha1,x,ztemp,xtemp,coef,xi,wi,sig) {
  
  nxi = length(xi); witemp = rep(wi,nrow(x)); xitemp = rep(xi,nrow(x))
  
  Int.1 = Int.2 = Int.3 = matrix(ncol=1,nrow=sum(R1))
  fxi.1 = fxi.2 = fxi.3 = rep(0,length(ztemp))
  fxi.4 = fxi.5 = fxi.6 = fxi.6.ps = t(ztemp)
  
  # Integral Denominator
  fxi.1 = pi.func(ztemp, alpha1)*sqrt(2*sig)*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.1 = matrix(colSums(fxi.1))
  
  # Integral 1 derivative wrt beta
  fxi.2 = pi.func(ztemp, alpha1)*2*xi*wi
  fxi.2 = matrix(c(fxi.2), nrow=nxi, byrow=FALSE)
  Int.2 = matrix(colSums(fxi.2))
  
  # Integral 2 derivative wrt beta
  fxi.3 = pi.func(ztemp, alpha1)*2*sqrt(2)*xi*xi*wi/sig
  fxi.3 = matrix(c(fxi.3), nrow=nxi, byrow=FALSE)
  Int.3 = matrix(colSums(fxi.3))
  
  a = matrix(rep(Int.1,nxi),ncol=nxi)
  b = t(a)
  b = as.vector(b)
  
  for (i in 1:ncol(ztemp)) {
    fxi.4[i,] = t(ztemp)[i,] * (wi*(pi.func(ztemp, alpha1)*(1-pi.func(ztemp, alpha1)))*sqrt(2*sig))
    fxi.5[i,] = t(ztemp)[i,] * (wi*(pi.func(ztemp, alpha1)-3*pi.func(ztemp, alpha1)^2+2*pi.func(ztemp, alpha1)^3)
                                *sqrt(2*sig)/b)
    fxi.6[i,] = t(ztemp)[i,] * ((wi*(pi.func(ztemp, alpha1)*(1-pi.func(ztemp, alpha1)))
                                 *2*xi)/b)
  }
  a = t(fxi.4)
  dim(a) = c(nxi,nrow(a)/nxi*ncol(a))
  b = as.vector(colSums(a))
  c = matrix(c(b),nrow=ncol(ztemp),byrow=TRUE)
  fxi.4 = c
  fxi.5 = fxi.5 %*% ztemp
  fxi.6 = fxi.6 %*% xtemp
  
  fxi.sigma = wi*prob*
    b*(1/sig)^(b+1)*(xi^2)^(b/2)*exp((1/sigma)^b*(-(xi^2)^(b/2)) + xi^2)
  fxi.sigma = matrix(c(fxi.sigma),nrow=nxi,byrow=FALSE)
  Int.S = matrix(colSums(fxi.sigma))
  
  list(Int.1=Int.1, Int.2=Int.2, Int.3=Int.3,
       Int.4=fxi.4, Int.5=fxi.5, Int.6=fxi.6)
}

S.0.opt = function(coef,R1,y,x,alpha1,xi=xi,wi=wi,data=data,sig=sig,desformula=desformula) {
  
  nxi = length(xi); witemp = rep(wi,nrow(x)); xitemp = rep(xi,nrow(x))
  results = tempfunction(data,x,xi,nxi,coef,desformula)
  ztemp = results$ztemp; xtemp = results$xtemp
  
  Int.1 = Int.2 = Int.3 = matrix(ncol=1,nrow=sum(R1))
  fxi.1 = fxi.2 = fxi.3 = rep(0,length(ztemp))
  
  # Integral Denominator
  fxi.1 = pi.func(ztemp, alpha1)*sqrt(2*sig)*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.1 = matrix(colSums(fxi.1))
  
  # Integral 1 derivative wrt beta
  fxi.2 = pi.func(ztemp, alpha1)*2*xi*wi
  fxi.2 = matrix(c(fxi.2), nrow=nxi, byrow=FALSE)
  Int.2 = matrix(colSums(fxi.2))  
  
  S = x*as.vector(y-x%*%coef) - x*as.vector(Int.2/Int.1)
  A = colSums(S)
  sum(A^2)
  #   list(Int.1=Int.1, Int.2=Int.2, factor=factor, S=S)
}

##########################
## EMPIRICAL COVARIANCE ##
##########################
empcov = function(x, y=NULL, nreps = rep(1, nrow(x)),center=TRUE)
{
  #   iid cov estimate that allows for relications
  #   inputs matrices x and y with same number of rows, nreps= replications vector of length=nrow(x)
  x=as.matrix(x); nreps=as.vector(nreps); if (!is.null(y)) y=as.matrix(y)
  if (length(nreps) != nrow(x)) stop("x, nreps length  mismatch")
  n = sum(nreps)
  xmean = numeric(ncol(x)); if (!is.null(y)) ymean = numeric(ncol(y))
  if (center){
    xmean = colSums(x*nreps, na.rm = FALSE)/n
    if (!is.null(y)) ymean = colSums(y*nreps, na.rm = FALSE)/n
  }
  xdiff = x - outer(rep(1, nrow(x)), xmean)
  if (is.null(y)) ydiff=xdiff
  else ydiff =y - outer(rep(1, nrow(y)), ymean)
  #    SS = t(xdiff) %*% diag(nreps) %*% ydiff
  r.nreps = sqrt(nreps)
  SS = t(xdiff*r.nreps) %*% (ydiff*r.nreps)
  cov = SS/(n-1)
  list(cov=cov,SS=SS,n=n)
}

na.keep = function(x) x

model.apply = function(formula,data,counts=NULL,always1=NULL){
  mf <- call <- match.call()
  mf[[1]] <- as.name("model.frame")
  names(mf)[1] <- "model"
  mf$na.action = as.name("na.keep")
  mf <- eval(mf[c("model","formula","data","counts","always1","na.action")], sys.frame(sys.parent()))
  Terms <- attr(mf,"terms")
  X <- model.matrix(Terms,mf)
  if(!is.matrix(X)) X <- matrix(X)
  counts <- model.extract(mf, counts)
  always1 <- model.extract(mf, always1)
  y <- model.extract(mf,response)
  list(y=y,X=X,counts=counts,always1=always1,Terms=Terms)
}

logfc <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use=use) {
  
  nxi = length(xi)
  witemp = rep(wi,nrow(x))
  xitemp = rep(xi,nrow(x))
  
  fxi.1 = pi.func(ztemp, alpha1)*wi/sqrt(pi)
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.0 = matrix(colSums(fxi.1))
  
  fy <- exp(-(y - x%*%coef)^2/(2*sig^2))*(1/sqrt(2*pi*sig^2))
  pi1 <-  plogis(z%*%alpha1)
  
  log(pi1[use]) -(y - x%*%coef)^2/(2*sig^2) - 2*log(sig) - log(Int.0)
}

S0func <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use=use) {
  
  nxi = length(xi)
  witemp = rep(wi,nrow(x))
  xitemp = rep(xi,nrow(x))
  
  fxi.1 = pi.func(ztemp, alpha1)*sqrt(2)*sig*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.0 = matrix(colSums(fxi.1))
  
  fxi.1 = pi.func(ztemp, alpha1)*2*xi*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.1 = matrix(colSums(fxi.1))
  
  x[use,]*as.vector((y[use] - x[use,]%*%coef)/(sig^2)) - x[use,]*as.vector(Int.1/Int.0)
}

S2func <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use=use) {
  
  nxi = length(xi)
  witemp = rep(wi,nrow(x))
  xitemp = rep(xi,nrow(x))
  
  fxi.1 = pi.func(ztemp, alpha1)*sqrt(2)*sig*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.0 = matrix(colSums(fxi.1))
  
  fxi.1 = pi.func(ztemp, alpha1)*xi*xi*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.3 = matrix(colSums(fxi.1))
  
  (y[use] - x[use,]%*%coef)^2/(sig^3) - as.vector(2*sqrt(2)*Int.3/Int.0)
}

Stilde2func <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use=use) {
  
  pi1 <-  plogis(z%*%alpha1)
  
  nxi = length(xi)
  witemp = rep(wi,nrow(x))
  xitemp = rep(xi,nrow(x))
  
  fxi.1 = pi.func(ztemp, alpha1)*sqrt(2)*sig*wi
  fxi.1 = matrix(c(fxi.1), nrow=nxi, byrow=FALSE)
  Int.0 = matrix(colSums(fxi.1))
  
  fxi.4 = t(ztemp)  
  for (i in 1:ncol(ztemp)) {
    fxi.4[i,] = t(ztemp)[i,] * (wi*(pi.func(ztemp, alpha1)*(1-pi.func(ztemp, alpha1)))*sqrt(2*sig^2))
  }
  a = t(fxi.4)
  dim(a) = c(nxi,nrow(a)/nxi*ncol(a))
  b = as.vector(colSums(a))
  c = matrix(c(b),nrow=ncol(ztemp),byrow=TRUE)
  fxi.4 = c
  
  Int.4 = as.matrix(fxi.4); ratioInt = matrix(nrow=sum(R1), ncol=nrow(fxi.4))
  for (f in 1:sum(R1))
    ratioInt[f,] = t(Int.4)[f,]*as.vector((1/Int.0))[f]
  
  (z*as.vector(1-pi1))[use,] - ratioInt
}


hessian_betabeta <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use, delta_step = 1e-4) {
  b <- matrix(NA, ncol = length(coef), nrow = length(coef))
  for (i in 1:length(coef)){
    coef_new1 <- coef_new2 <- coef
    coef_new1[i] <- coef_new1[i] + delta_step
    coef_new2[i] <- coef_new2[i]
    a1 <- S0func(y,x,z,alpha1,ztemp,xtemp,coef_new1,xi,wi, sig,use)
    a0 <- S0func(y,x,z,alpha1,ztemp,xtemp,coef_new2,xi,wi, sig,use)
    b[i,] <- colSums((a1 - a0)/(delta_step))
  }
  b
}
hessian_sigmasigma <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use, delta_step = 1e-4) {
  b <- matrix(NA, ncol = length(sig), nrow = length(sig))
  for (i in 1:length(sig)){
    sig_new1 <- sig_new2 <- sig
    sig_new1[i] <- sig_new1[i] + delta_step
    sig_new2[i] <- sig_new2[i]
    a1 <- S2func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig_new1,use)
    a0 <- S2func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig_new2,use)
    b[i,] <- colSums((a1 - a0)/(delta_step))
  }
  b
}
hessian_alphaalpha <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use, delta_step = 1e-4) {
  b <- matrix(NA, ncol = length(alpha1), nrow = length(alpha1))
  for (i in 1:length(alpha1)){
    alpha1_new1 <- alpha1_new2 <- alpha1
    alpha1_new1[i] <- alpha1_new1[i] + delta_step
    alpha1_new2[i] <- alpha1_new2[i]
    a1 <- Stilde2func(y,x,z,alpha1_new1,ztemp,xtemp,coef,xi,wi, sig,use)
    a0 <- Stilde2func(y,x,z,alpha1_new2,ztemp,xtemp,coef,xi,wi, sig,use)
    b[i,] <- colSums((a1 - a0)/(delta_step))
  }
  b
}
hessian_betasigma <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use, delta_step = 1e-4) {
  b <- matrix(NA, ncol = length(beta), nrow = length(sig))
  for (i in 1:length(sig)){
    sig_new1 <- sig_new2 <- sig
    sig_new1[i] <- sig_new1[i] + delta_step
    sig_new2[i] <- sig_new2[i]
    a1 <- S0func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig_new1,use)
    a0 <- S0func(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig_new2,use)
    b[i,] <- colSums((a1 - a0)/(delta_step))
  }
  b
}
hessian_betaalpha <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use, delta_step = 1e-4) {
  b <- matrix(NA, ncol = length(beta), nrow = length(alpha1))
  for (i in 1:length(alpha1)){
    alpha_new1 <- alpha_new2 <- alpha1
    alpha_new1[i] <- alpha_new1[i] + delta_step
    alpha_new2[i] <- alpha_new2[i]
    a1 <- S0func(y,x,z,alpha_new1,ztemp,xtemp,coef,xi,wi, sig,use)
    a0 <- S0func(y,x,z,alpha_new2,ztemp,xtemp,coef,xi,wi, sig,use)
    b[i,] <- colSums((a1 - a0)/(delta_step))
  }
  b
}
hessian_sigmaalpha <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use, delta_step = 1e-4) {
  b <- matrix(NA, ncol = length(sig), nrow = length(alpha1))
  for (i in 1:length(alpha1)){
    alpha_new1 <- alpha_new2 <- alpha1
    alpha_new1[i] <- alpha_new1[i] + delta_step
    alpha_new2[i] <- alpha_new2[i]
    a1 <- S2func(y,x,z,alpha_new1,ztemp,xtemp,coef,xi,wi, sig,use)
    a0 <- S2func(y,x,z,alpha_new2,ztemp,xtemp,coef,xi,wi, sig,use)
    b[i,] <- colSums((a1 - a0)/(delta_step))
  }
  b
}





score_beta <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi,sig,use,h = 1e-10) {
  b <- matrix(NA, ncol = length(beta), nrow = nrow(x))
  for (i in 1:length(coef)){
    beta_new1 <- beta_new2 <- coef
    beta_new1[i] <- beta_new1[i] + h
    beta_new2[i] <- beta_new2[i] - h
    a1 <- logfc(y,x,z,alpha1,ztemp,xtemp,beta_new1,xi,wi, sig,use)
    a0 <- logfc(y,x,z,alpha1,ztemp,xtemp,beta_new2,xi,wi, sig,use)
    b[,i] <- (a1 - a0)/(2*h)
  }
  b
}

score_alpha <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use,h = 1e-10) {
  b <- matrix(NA, ncol = length(alpha1), nrow = nrow(x))
  for (i in 1:length(alpha1)){
    alpha_new1 <- alpha_new2 <- alpha1
    alpha_new1[i] <- alpha_new1[i] + h
    alpha_new2[i] <- alpha_new2[i] - h
    a1 <- logfc(y,x,z,alpha_new1,ztemp,xtemp,coef,xi,wi, sig,use)
    a0 <- logfc(y,x,z,alpha_new2,ztemp,xtemp,coef,xi,wi, sig,use)
    b[,i] <- (a1 - a0)/(2*h)
  }
  b
}

score_sigma <- function(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig,use,h = 1e-10) {
  b <- matrix(NA, ncol = length(sig), nrow = nrow(x))
  for (i in 1:length(sig)){
    sig_new1 <- sig_new2 <- sig
    sig_new1[i] <- sig_new1[i] + h
    sig_new2[i] <- sig_new2[i] - h
    a1 <- logfc(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig_new1,use)
    a0 <- logfc(y,x,z,alpha1,ztemp,xtemp,coef,xi,wi, sig_new2,use)
    b[,i] <- (a1 - a0)/(2*h)
  }
  b
}
