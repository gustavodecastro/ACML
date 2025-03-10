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
    x = x[use,]; y = y[use]; pi11 = pi1[use]
    
    I00 = t(x) %*% diag(as.vector(1/pi11)) %*% x
    I01 = t(x) %*% diag(as.vector((y-x%*%coef)*(1-pi11)*(1/pi11))) %*% w[R1==1,]
    w.aux = w*as.vector(sqrt(pi1*(1-pi1)))
    I11 = t(w.aux) %*% w.aux
    I10 = matrix(rep(0,ncol(I00)*nrow(I11)),nrow=nrow(I11),ncol=ncol(I00))
    I = rbind(cbind(I00,I01),
              cbind(I10,I11))
    S0 = x*as.vector((y-x%*%coef)/pi11)
    S1 = w*as.vector(R1-pi1)
    S = t(cbind(t(colSums(S0)), t(colSums(S1))))
    
    C00 = empcov(S0,nreps=counts[use])$SS
    C01 = empcov(x=S0,y=S1[use,],nreps=counts[use])$SS
    C11 = empcov(S1,nreps=counts)$SS
    C  = rbind(cbind(C00,C01),cbind(t(C01),C11))
    ser = sqrt(diag(solve(I,t(solve(I,C))))[1:ncol(x)])
    
    list(coef=coef,ser=ser)
    
    #################################################################
  } else if (method == "cml") { ############################# cml
    
    coef = lm(y ~ x-1, weights=1/pi1, subset=use)$coef; coef.fixed = coef
    
    sig <- c(sqrt(var(y[use] - x[use,]%*%coef)))
    data$sig <- sig
    nxi = length(xi)
    results = tempfunction(data,x,xi,nxi,coef,desformula)
    ztemp = results$ztemp; xtemp = results$xtemp
    Integrais = Integral(alpha1,x,ztemp,xtemp,coef,xi,wi, sig)
    
    #    x = x[use,]
    pi11 = pi1[use]
    phi = c(coef, alpha1); phi.fixed = phi; coef.fixed = coef
    phi.old = step = matrix(ncol = length(phi), nrow = niter)
    S0.old = matrix(ncol = length(beta), nrow = niter)
    Cycle = FALSE; Repeated.step = FALSE; OPT = FALSE; Div = 1
    conv = FALSE
    
    
    for (jj in 1:niter) {
      if (jj>1) {
        coef = phi[1:ncol(x)]
        alpha1 = phi[(ncol(x)+1):(ncol(x)+length(alpha1))]; pi1 = plogis(z%*%alpha1)
        results = tempfunction(data,x,xi,nxi,coef,desformula)
        ztemp = results$ztemp; xtemp = results$xtemp
        Integrais = Integral(alpha1,x,ztemp,xtemp,coef,xi,wi, sig)
      }
      sig <- c(sqrt(var(y[use] - x[use,]%*%coef)))
      data$sig <- sig
      
      I00 = t(x[use,]) %*% diag(as.vector(Integrais$Int.3/Integrais$Int.1)) %*% x[use,] -
        t(x[use,]) %*% diag(as.vector(Integrais$Int.2/Integrais$Int.1)*(as.vector(Integrais$Int.2/Integrais$Int.1))) %*% x[use,]
      I01 = t(Integrais$Int.6) -
        t(x[use,]) %*% diag(as.vector(Integrais$Int.2/Integrais$Int.1/Integrais$Int.1)) %*% t(Integrais$Int.4)
      colnames(I01) = colnames(z)
      y.aux = z*as.vector(sqrt(pi1*(1-pi1)))
      I11 = t(y.aux) %*% y.aux
      I10 = matrix(c(rep(0,nrow(I11)*ncol(I00))), ncol=ncol(I00), nrow=nrow(I11))
      
      if (Stilde) {
        Itilde10 = Integrais$Int.6 -
          Integrais$Int.4 %*% diag(as.vector(Integrais$Int.2 * 1/Integrais$Int.1 * 1/Integrais$Int.1)) %*% x[use,]
        rownames(Itilde10) = colnames(z)
        z.aux = z*as.vector(sqrt(pi1*(1-pi1)))
        Itilde11 =  t(z.aux[use,]) %*% z.aux[use,] +
          Integrais$Int.5 -
          Integrais$Int.4 %*% diag(as.vector(1/Integrais$Int.1/Integrais$Int.1)) %*% t(Integrais$Int.4)
        I11 = I11 - Itilde11; I10 = I10 - Itilde10
      }
      
      II = rbind(cbind(I00,I01),
                 cbind(I10,I11))
      S0 = x[use,]*as.vector(y[use]-x[use,]%*%coef) - x[use,]*as.vector(Integrais$Int.2/Integrais$Int.1)
      S1 = z*as.vector(R1-pi1)
      
      if (Stilde) {
        Integrais$Int.4 = as.matrix(Integrais$Int.4); fraction = matrix(nrow=sum(R1), ncol=nrow(Integrais$Int.4))
        for (f in 1:sum(R1)) fraction[f,] = t(Integrais$Int.4)[f,]*as.vector((1/Integrais$Int.1))[f]
        Stilde1 = (z*as.vector(1-pi1))[use] - fraction
        S = t(cbind(t(colSums(S0)), t(colSums(S1) - colSums(Stilde1))))
      } else S = t(cbind(t(colSums(S0)), t(colSums(S1))))
      
      Step = 0
      if (niter >1) {
        Step = solve(II,S)
      }
      bigstep = max(abs(Step))
      maxstep = maxstepprop*max(abs(phi))
      if (bigstep > maxstep) {
        Step = Step * maxstep/bigstep
      }
      step[jj,] = Step; phi.old[jj,] = phi; S0.old[jj,] = colSums(S0); #print(paste("S0 = ",colSums(S0)))
      
      if (jj > 5) { 
        if (!Cycle) {
          for (jjj in 1:(jj-4)) { ##Check if we are stuck in a loop
            if (all(round(step[jj,1:3],5) == round(step[jj-jjj,1:3],5)) & max(abs(Step[1:3]))>1e-4) {
              #              Repeated.step = TRUE
              atJ = jj; Cycle = TRUE
            }
          }
        }
        if (Repeated.step) { ##Inside a loop, going back and forth
          print(paste("Back and forth..."))
          S0.test = -S0.old[jj-1,]; Step.test = Step
          print(paste("Optimizing"))
          S0.opt = optim(coef,S.0.opt,R1=R1,y=y,x=x,alpha1=alpha1,xi=xi,wi=wi,sig=sig,data=data,desformula=desformula)
          S0.conv = S0.opt$conv
          if (S0.conv == 0) conv=TRUE
          coef = S0.opt$par
          phi = c(phi,alpha1)
          break
        } else { ##If there are NO repeated steps, is not going back and forth
          phi2 = phi + Step
          if (niter>1) if (max(abs(phi-phi2)[1:ncol(x)]) < tol*max(abs(phi)[1:ncol(x)])) {
            conv = TRUE
            break
          }
          phi = phi2
        } ##Finishing checking if there are or there aren't repeated steps
      } else { ##If jj <= 5
        phi2 = phi + Step
        if (niter>1) if (max(abs(phi-phi2)[1:ncol(x)]) < tol*max(abs(phi)[1:ncol(x)])) {
          conv = TRUE
          break
        }
        phi = phi2
      }
      #print(paste("jj=",jj," coef=",round(phi,4)," coef.fixed=",round(phi.fixed,4)," Step=",round(Step,7)))
    } ##End loop
    
    coef = phi[1:ncol(x)]
    if (niter>1 & !conv) {
      error = error + 1
      Error[ii] = ii
      print(paste("Erros = ",error))
    }
    
    C00 = empcov(S0,nreps=counts[use])$SS
    
    if (Stilde) {
      C01 = empcov(x=S0,y=S1[use,],nreps=counts[use])$SS - empcov(x=S0,Stilde1,nreps=counts[use])$SS
      C11 = empcov(S1,nreps=counts)$SS - 2*empcov(x=S1[use,],Stilde1,nreps=counts[use])$SS + empcov(Stilde1,nreps=counts[use])$SS
    } else {
      C01 = empcov(x=S0,y=S1[use,],nreps=counts[use])$SS
      C11 = empcov(S1,nreps=counts)$SS
    }
    C   = rbind(cbind(C00,C01),cbind(t(C01),C11))
    ser = sqrt(diag(solve(II,t(solve(II,C))))[1:ncol(x)])
    
    list(coef=coef,ser=ser,conv=conv)
    
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
