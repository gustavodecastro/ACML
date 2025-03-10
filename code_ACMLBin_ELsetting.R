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
# model.apply(y~x1d+x2+x1d*x2,datcsv)


func <- function(regformula, desnformula1, desnformula2=NULL, desnformula3=NULL, data, method="offset", srs=NULL, addStilde=FALSE, niter=1, tol=1.0e-6, maxstepprop=0.1, print.iters=FALSE){
    # always1: TRUE if sampled with probability 1 at design phase, FALSE else
    error = 0; ow <- options("warn"); mzero=function(n,m) matrix(0,n,m) # utility function

 #  Model of interest
    M1 = model.apply(regformula,data,counts,always1)
    if (attr(M1$Terms,"response") != 1) stop("regformula must have a response variable")
    y = M1$y; X=M1$X; counts=M1$counts; always1=M1$always1

 #  Design Model 1
    M2 = model.apply(desnformula1,data,counts,always2)
    if (attr(M2$Terms,"response") != 1) stop("desnformula1 must have a response variable")
    R1 = M2$y; z1=M2$X; always2=M2$always1; R1aux = R1; z1aux = z1

    if (!is.null(srs)) {
      use=(R1==1); use.srs = (R1.srs==1)
    }

 #  Design Model 2
    if (!is.null(desnformula2)) {
      M3 = model.apply(desnformula2,data,counts)
      if (attr(M3$Terms,"response") != 1) stop("desnformula2 must have a response variable")
      R2 = M3$y; z2=as.matrix(M3$X)
    }

 #  Design Model 3
    if (!is.null(desnformula3)) {
      M4 = model.apply(desnformula3,data,counts)
      if (attr(M4$Terms,"response") != 1) stop("desnformula3 must have a response variable")
      R3 = M4$y; z3=M4$X
    }

    use1 = !always1; fit1 = fit2 = fit3 = 0; use1 = rep(TRUE,length(R1))
    if (!is.null(srs)){
      fit1.srs = glm(R1.srs ~ 1, family=binomial, weights=counts, subset=use1, data=data)
      alpha1.srs = as.vector(coefficients(fit1.srs)); pi1.srs = plogis(1 * alpha1.srs)
#      pi1 = pi1*(1-pi1.srs)
      pi1[srs] = pi1.srs
      use1 = !use.srs
    }
    fit1 = glm(R1 ~ z1-1, family=binomial, weights=counts, subset=use1, data=data)
    alpha1 = as.vector(coefficients(fit1)); names(alpha1) = colnames(z1); pi1 = plogis(z1 %*% alpha1)

    if (!is.null(desnformula2)) {
      use2 = (R1==1 & !is.na(R2))
      fit2 = glm(R2 ~ z2-1, family=binomial, weights=counts, subset=use2, data=data)
      alpha2 = as.vector(coefficients(fit2)); names(alpha2) = colnames(z2); pi2 = plogis(z2 %*% alpha2)
    }

    if (!is.null(desnformula3)) {
      use3 = (R1==1 & R2==1)
      fit3 = glm(R3 ~ z3-1, family=binomial, weights=counts, subset=use3, data=data)
      alpha3 = as.vector(coefficients(fit3)); pi3 = plogis(z3 %*% alpha3)
    }

    if (method=="offset") {
      yvarname= names(M1[1])
      yvarname= "y"
      temdata=data
      temdata[,yvarname] = rep(1,nrow(data)) # Make y-val always =1
      z1.1 = model.apply(desnformula1,temdata,counts)$X
      if (!is.null(desnformula2)) z2.1 = model.apply(desnformula2,temdata,counts)$X
      if (!is.null(desnformula3)) z3.1 = model.apply(desnformula3,temdata,counts)$X
      temdata[,yvarname] = rep(0,nrow(data)) # Make y-val always =0
      z1.0 = model.apply(desnformula1,temdata,counts)$X
      if (!is.null(desnformula2)) z2.0 = model.apply(desnformula2,temdata,counts)$X
      if (!is.null(desnformula3)) z3.0 = model.apply(desnformula3,temdata,counts)$X

      lpi1.1 = plogis(z1.1 %*% alpha1, log.p=TRUE)
      lpi1.0 = plogis(z1.0 %*% alpha1, log.p=TRUE)
      oo = lpi1.1 - lpi1.0
      if (!is.null(srs)) {
        lpi1.1[srs] = 0
        lpi1.0[srs] = 0
        oo = lpi1.1 - lpi1.0
      }
      if (!is.null(desnformula2)) {
        lpi2.1 = plogis(z2.1 %*% alpha2, log.p=TRUE)
        lpi2.0 = plogis(z2.0 %*% alpha2, log.p=TRUE)
        oo = lpi1.1 + lpi2.1 - lpi1.0 - lpi2.0
      }
      if (!is.null(desnformula3)) {
        lpi3.1 = plogis(z3.1 %*% alpha3, log.p=TRUE)
        lpi3.0 = plogis(z3.0 %*% alpha3, log.p=TRUE)
        oo = lpi1.1 + lpi2.1 + lpi3.1 - lpi1.0 - lpi2.0 - lpi3.0
      }
      weights = counts
    } else if (method=="weighted"){
      oo = rep(0,nrow(data))
      pi1pi2pi3 = as.vector(pi1)
      if (!is.null(desnformula2)) pi1pi2pi3 = as.vector(pi1*pi2)
      if (!is.null(desnformula3)) pi1pi2pi3 = as.vector(pi1*pi2*pi3)
      weights = R1/pi1pi2pi3 + (1-R1)/(1-pi1pi2pi3)
      options(warn = -1) # supress warnings (complains about non-integer weights)
    } else stop("Method not available")

    use0 = (R1==1)# those who have full data
    if (!is.null(desnformula2)) use0 = (R1==1 & R2==1)
    if (!is.null(desnformula3)) use0 = (R1==1 & R2==1 & R3==1)
    if (!is.null(srs)) {
      use0 = (R1==1) | R1.srs==1
      if (!is.null(desnformula2)) use0 = (R1==1 & R2==1) | R1.srs==1
      if (!is.null(desnformula3)) use0 = (R1==1 & R2==1 & R3==1) | R1.srs==1
    }
    fit = glm(y ~ X-1, family=binomial, weights=weights, offset=oo, subset=use0)
    coefs = coefficients(fit); names(coefs)=colnames(X)
    options(ow) # reset    

    error = 0; conv = FALSE
    n0=ncol(X); n1=ncol(z1)
    phi = c(coefs,alpha1)

    if (!is.null(desnformula2)) {
      n2=ncol(z2)
      phi = c(coefs,alpha1,alpha2)
    }
    if (!is.null(desnformula3)) {
      n3=ncol(z3)
      phi = c(coefs,alpha1,alpha2,alpha3)
    }

    for (iter in 1:niter){ # niter should = 1 except when addStilde=TRUE
 
      if (addStilde) { # Recalculate some basic quantities
        coefs = phi[1:n0]
        alpha1 = phi[n0+(1:n1)]
        pi1 = plogis(z1 %*% alpha1)
        if (!is.null(srs)) {
#          pi1 = pi1*(1-pi1.srs)
          pi1[srs] = pi1.srs
        } 
        if (sum(always1)>0) {
          alpha1.new = alpha1
          z1 = z1.old; z1.0 = z1.0.old; z1.1 = z1.1.old
          z1 = as.matrix(z1); z1.1 = as.matrix(z1.1); z1.0 = as.matrix(z1.0)
          pi1 = plogis(z1 %*% alpha1)
        }

        if (!is.null(desnformula2)) {
          alpha2 = phi[n0+n1+(1:n2)]
          if (sum(always2)>0) {
            alpha2.new = alpha2
            alpha2.old[c(eq)] = alpha2.new; alpha2 = alpha2.old
            z2 = z2.old; z2.0 = z2.0.old; z2.1 = z2.1.old
          }
          pi2 = plogis(z2 %*% alpha2)
        }
        if (!is.null(desnformula3)) {
          alpha3 = phi[n0+n1+n2+(1:n3)]
          pi3 = plogis(z3 %*% alpha3)
        }
        if (method=="offset"){
          lpi1.1 = plogis(z1.1 %*% alpha1, log.p=TRUE)
          lpi1.0 = plogis(z1.0 %*% alpha1, log.p=TRUE)
          oo = lpi1.1 - lpi1.0
          if (!is.null(srs)) {
            lpi1.1 = plogis(z1.1[!use.srs,] %*% alpha1, log.p=TRUE)
            lpi1.0 = plogis(z1.0[!use.srs,] %*% alpha1, log.p=TRUE)

            lpi1.1.aux = rep(log(pi1.srs),nrow(data))
            lpi1.1.aux[!use.srs] = lpi1.1
            lpi1.1 = lpi1.1.aux

            lpi1.0.aux = rep(log(pi1.srs),nrow(data))
            lpi1.0.aux[!use.srs] = lpi1.0
            lpi1.0 = lpi1.0.aux

            oo = lpi1.1 - lpi1.0
          }
          if (!is.null(desnformula2)) {
            lpi2.1 = plogis(z2.1 %*% alpha2, log.p=TRUE)
            lpi2.0 = plogis(z2.0 %*% alpha2, log.p=TRUE)
            oo = lpi1.1 + lpi2.1 - lpi1.0 - lpi2.0
          }
          if (!is.null(desnformula3)) {
            lpi3.1 = plogis(z3.1 %*% alpha3, log.p=TRUE)
            lpi3.0 = plogis(z3.0 %*% alpha3, log.p=TRUE)
            oo = lpi1.1 + lpi2.1 + lpi3.1 - lpi1.0 - lpi2.0 - lpi3.0
          }
        }
        if (sum(always1)>0) {
          z1 = z1.new; z1.0 = z1.0.new; z1.1 = z1.1.new; alpha1 = alpha1.new
          z1 = as.matrix(z1); z1.1 = as.matrix(z1.1); z1.0 = as.matrix(z1.0)
        }
        if (!is.null(desnformula2)&sum(always2)>0) {
          z2 = z2.new; z2.0 = z2.0.new; z2.1 = z2.1.new; alpha2 = alpha2.new
        }
      }
      S1  = z1*(R1-as.vector(pi1)) #colSums(S1[use1,])
      if (!is.null(srs)) S1[srs] = 0
      pi1vec = as.vector(pi1)
      if (!is.null(desnformula2)) {
        S2  = z2*(R2-as.vector(pi2)) #colSums(S2[use2,])
        pi2vec = as.vector(pi2)
      }
      if (!is.null(desnformula3)) {
        S3  = z3*(R3-as.vector(pi3)) #colSums(S3[use3,])
        pi3vec = as.vector(pi3)
      }

      if (method=="offset"){
        pstar = as.vector(plogis(oo+X%*%coefs))
        S0    = X*(y-pstar)
        I00   = t(X[use0,,drop=FALSE]) %*% diag((counts[use0]*pstar[use0]*(1-pstar[use0]))) %*% X[use0,,drop=FALSE]
        pi1.1     = as.vector(exp(lpi1.1)); pi1.0 = as.vector(exp(lpi1.0))
        dodalpha1 = as.matrix((z1.1*(1-pi1.1) - z1.0*(1-pi1.0)))
        if (!is.null(srs)) dodalpha1[srs] = 0
        I01   = t(X[use0,,drop=FALSE]) %*% diag((counts[use0]*pstar[use0]*(1-pstar[use0]))) %*% dodalpha1[use0,,drop=FALSE]
        z1aux = z1*sqrt(counts*pi1vec*(1-pi1vec))
        if (!is.null(srs)) z1aux = z1aux[!use.srs,]
        I11   = t(z1aux) %*% z1aux
        I10 = mzero(n1,n0)
        II  = rbind(cbind(I00,I01),
                    cbind(I10,I11))
        S = t(cbind(t(colSums((counts*S0)[use0,,drop=FALSE])),
            t(colSums((counts*S1)[use1,,drop=FALSE])) ))

        if (!is.null(desnformula2)) {
          pi2.1     = as.vector(exp(lpi2.1)); pi2.0 = as.vector(exp(lpi2.0))
          dodalpha2 = (z2.1*(1-pi2.1) - z2.0*(1-pi2.0))
          I02   = t(X[use0,,drop=FALSE]) %*% diag((counts[use0]*pstar[use0]*(1-pstar[use0]))) %*% dodalpha2[use0,,drop=FALSE]
          z2aux = (z2*sqrt(counts*pi2vec*(1-pi2vec)))[use2,]
          I22   = t(z2[use2,,drop=FALSE]) %*% diag((counts*pi2*(1-pi2))[use2]) %*% z2[use2,,drop=FALSE]
          I20 = mzero(n2,n0); I12 = mzero(n1,n2); I21 = t(I12)
          II  = rbind(cbind(I00,I01,I02),
                      cbind(I10,I11,I12),
                      cbind(I20,I21,I22))
          S = t(cbind(t(colSums((counts*S0)[use0,,drop=FALSE])),
              t(colSums((counts*S1)[use1,,drop=FALSE])),
              t(colSums((counts*S2)[use2,,drop=FALSE])) ))
        }
        if (!is.null(desnformula3)) {
          pi3.1     = as.vector(exp(lpi3.1)); pi3.0 = as.vector(exp(lpi3.0))
          dodalpha3 = (z3.1*(1-pi3.1) - z3.0*(1-pi3.0))
          I03   = t(X[use0,,drop=FALSE]) %*% diag((counts[use0]*pstar[use0]*(1-pstar[use0]))) %*% dodalpha3[use0,,drop=FALSE]
          z3aux = (z3*sqrt(counts*pi3vec*(1-pi3vec)))[use3,]
          I33   = t(z3aux) %*% z3aux
          I30 = mzero(n3,n0); I13 = mzero(n1,n3);  I23 = mzero(n2,n3); I31 = t(I13); I32 = t(I23)
          II  = rbind(cbind(I00,I01,I02,I03),
                      cbind(I10,I11,I12,I13),
                      cbind(I20,I21,I22,I23),
                      cbind(I30,I31,I32,I33))
          S = t(cbind(t(colSums((counts*S0)[use0,,drop=FALSE])),
              t(colSums((counts*S1)[use1,,drop=FALSE])),
              t(colSums((counts*S2)[use2,,drop=FALSE])),
              t(colSums((counts*S3)[use3,,drop=FALSE])) ))
        }

        if (addStilde){
          use0 = (R1==1)
          Stilde1 = dodalpha1*(y-pstar) # Subsetting done later when required
          Itilde10 = t(dodalpha1[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% X[use0,,drop=FALSE]
          Itilde11 = t(dodalpha1[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% dodalpha1[use0,,drop=FALSE] +
                     t(z1.1[use0,,drop=FALSE]) %*% diag((counts*pi1.1*(1-pi1.1)*(y-pstar))[use0]) %*% z1.1[use0,,drop=FALSE] -
                     t(z1.0[use0,,drop=FALSE]) %*% diag((counts*pi1.0*(1-pi1.0)*(y-pstar))[use0]) %*% z1.0[use0,,drop=FALSE]
          I10 = I10 - Itilde10; I11 = I11 - Itilde11; 
          II  = rbind(cbind(I00,I01),
                      cbind(I10,I11))
          S = t(cbind(t(colSums((counts*S0)[use0,,drop=FALSE])),
              t(colSums((counts*S1)[use1,,drop=FALSE]) - colSums((counts*Stilde1)[use0,,drop=FALSE])) ))

          if (!is.null(desnformula2)) {
            Stilde2 = dodalpha2*(y-pstar) # Subsetting done later when required
            Itilde20 = t(dodalpha2[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% X[use0,,drop=FALSE]
            Itilde12 = t(dodalpha1[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% dodalpha2[use0,,drop=FALSE]
            Itilde22 = t(dodalpha2[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% dodalpha2[use0,,drop=FALSE] +
                       t(z2.1[use0,,drop=FALSE]) %*% diag((counts*pi2.1*(1-pi2.1)*(y-pstar))[use0]) %*% z2.1[use0,,drop=FALSE] -
                       t(z2.0[use0,,drop=FALSE]) %*% diag((counts*pi2.0*(1-pi2.0)*(y-pstar))[use0]) %*% z2.0[use0,,drop=FALSE]
            Itilde21 = t(Itilde12)
            I12 = I12 - Itilde12; I20 = I20 - Itilde20; I21 = I21 - Itilde21; I22 = I22 - Itilde22
            II  = rbind(cbind(I00,I01,I02),
                        cbind(I10,I11,I12),
                        cbind(I20,I21,I22))
            S = t(cbind(t(colSums((counts*S0)[use0,,drop=FALSE])),
                t(colSums((counts*S1)[use1,,drop=FALSE]) - colSums((counts*Stilde1)[use0,,drop=FALSE])),
                t(colSums((counts*S2)[use2,,drop=FALSE]) - colSums((counts*Stilde2)[use0,,drop=FALSE])) ))
          }

          if (!is.null(desnformula3)) {
            Stilde3 = dodalpha3*(y-pstar) # Subsetting done later when required
            Itilde30 = t(dodalpha3[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% X[use0,,drop=FALSE]
            Itilde13 = t(dodalpha1[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% dodalpha3[use0,,drop=FALSE]
            Itilde23 = t(dodalpha2[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% dodalpha3[use0,,drop=FALSE]
            Itilde33 = t(dodalpha3[use0,,drop=FALSE]) %*% diag((counts*pstar*(1-pstar))[use0]) %*% dodalpha3[use0,,drop=FALSE] +
                       t(z3.1[use0,,drop=FALSE]) %*% diag((counts*pi3.1*(1-pi3.1)*(y-pstar))[use0]) %*% z3.1[use0,,drop=FALSE] -
                       t(z3.0[use0,,drop=FALSE]) %*% diag((counts*pi3.0*(1-pi3.0)*(y-pstar))[use0]) %*% z3.0[use0,,drop=FALSE]
            Itilde31 = t(Itilde13); Itilde32 = t(Itilde23)
            I30 = I30 - Itilde30; I31 = I31 - Itilde31; I32 = I32 - Itilde32; I33 = I33 - Itilde33; I13 = I13 - Itilde13; I23 = I23 - Itilde23
            II  = rbind(cbind(I00,I01,I02,I03),
                        cbind(I10,I11,I12,I13),
                        cbind(I20,I21,I22,I23),
                        cbind(I30,I31,I32,I33))
            S = t(cbind(t(colSums((counts*S0)[use0,,drop=FALSE])),
                t(colSums((counts*S1)[use1,,drop=FALSE]) - colSums((counts*Stilde1)[use0,,drop=FALSE])),
                t(colSums((counts*S2)[use2,,drop=FALSE]) - colSums((counts*Stilde2)[use0,,drop=FALSE])),
                t(colSums((counts*S3)[use3,,drop=FALSE]) - colSums((counts*Stilde3)[use0,,drop=FALSE])) ))
          }
          if (!is.null(srs)) use0 = (R1==1) | (R1.srs==1)
        } #End addStilde
      } else if (method=="weighted"){
        p   = as.vector(plogis(X%*%coefs))
        S0  = X*((y-p)/pi1pi2pi3) # colSums(S0[use0,,drop=FALSE])
        I10 = mzero(n1,n0)
        I00 = t(X[use0,,drop=FALSE]) %*% diag((counts*p*(1-p)/pi1pi2pi3)[use0]) %*% X[use0,,drop=FALSE]
        I01 = t(S0[use0,,drop=FALSE]) %*% diag((counts*(1-pi1vec))[use0]) %*% z1[use0,,drop=FALSE]
        z1aux = z1*sqrt(counts*pi1vec*(1-pi1vec))
        I11 = t(z1aux) %*% z1aux
        if (!is.null(desnformula2)) {
          I20 = mzero(n2,n0); I12 = mzero(n1,n2); I21 = t(I12)
          I02 = t(S0[use0,,drop=FALSE]) %*% diag((counts*(1-pi2vec))[use0]) %*% z2[use0,,drop=FALSE]
          z2aux = z2*sqrt(counts*pi2vec*(1-pi2vec))
          I22 = t(z2aux[use2,]) %*% z2aux[use2,]
        }
        if (!is.null(desnformula3)) {
          I30 = mzero(n3,n0); I13 = mzero(n1,n3);  I23 = mzero(n2,n3); I31 = t(I13); I32 = t(I23)
          I03 = t(S0[use0,,drop=FALSE]) %*% diag((counts*(1-pi3vec))[use0]) %*% z3[use0,,drop=FALSE]
          z3aux = z3*sqrt(counts*pi3vec*(1-pi3vec))
          I33 = t(z3aux) %*% z3aux
        }
      } else stop("Method not available")

      Step = 0
      if (niter >1) Step = solve(II,S)
      bigstep = max(abs(Step))
      maxstep = maxstepprop*max(abs(phi))
      if (bigstep > maxstep) {
          Step = Step * maxstep/bigstep
          if (print.iters) print(paste("Step shortened by",maxstep/bigstep))
      }

      phi2 = phi + Step
      if (print.iters & niter>1) {
        print(paste("Iter", iter))
        print(paste("phi= ",phi))
        print(paste("Step= ",Step))
        print(paste("S=", S))
      }

      if (niter>1) if (max(abs(phi-phi2)[1:ncol(X)]) < tol*max(abs(phi[1:ncol(X)]))){#& max(abs(S)[1:ncol(X)]) < tol) {
#      if (niter>1) if (max(abs(phi-phi2)) < tol*max(abs(phi)) ){#& max(abs(S)[1:ncol(X)]) < tol) {
          conv = TRUE
          break
      }
      phi = phi2
    } # END OF ITERATION LOOP
    if (niter>1 & !conv) return(invisible(list(msg="Failed to Converge", phi=phi,S=S,II=II,error=1,pi1=pi1)))

    if (addStilde) {
        dimS1 = dim(S1)
        lines0 = (1:nrow(S1))[use0]
        S1tilde = matrix(rep(0, dim(S1)[1]*dim(S1)[2]), nrow = dim(S1)[1])
        if (!is.null(srs)) Stilde1[srs] = 0
        S1tilde[lines0,] = Stilde1[use0,,drop=FALSE]
        names = colnames(S1)
        S1 = S1 - S1tilde #colSums(S1[use1,])
        colnames(S1) = names

        if (!is.null(desnformula2)) {
          lines2 = (1:nrow(S1))[use2]
          S2star = matrix(rep(0, dim(S1)[1]*dim(S2)[2]), nrow = dim(S1)[1])
          S2star[lines2,] = S2[use2,,drop=FALSE]
          linesS2 = (1:nrow(S1))[use0]
          S2tilde = matrix(rep(0, dim(S1)[1]*dim(S2star)[2]), nrow = dim(S1)[1])
          S2tilde[linesS2,] = Stilde2[use0,,drop=FALSE]
          S2 = S2star - S2tilde #colSums(S2[use2,])
          colnames(S2) = colnames(Stilde2)
        }
    }

    C00 = empcov(S0[use0,,drop=FALSE],nreps=counts[use0])$SS
    C01 = empcov(x=S0[use0,,drop=FALSE],y=S1[use0,,drop=FALSE],nreps=counts[use0])$SS
    C11 = empcov(S1[use1,,drop=FALSE],nreps=counts[use1])$SS
    CC = rbind(cbind(C00,C01),cbind(t(C01),C11))
    if (!is.null(desnformula2)) {
      C02 = empcov(x=S0[use0,,drop=FALSE],y=S2[use0,,drop=FALSE],nreps=counts[use0])$SS
      C12 = empcov(x=S1[use2,,drop=FALSE],y=S2[use2,,drop=FALSE],nreps=counts[use2])$SS
      C22 = empcov(S2[use2,,drop=FALSE],nreps=counts[use2])$SS
      CC = rbind(cbind(C00,C01,C02),cbind(t(C01),C11,C12),cbind(t(C02),t(C12),C22))
    }
    if (method=="weighted") {
      Mid = C00 - I01%*%solve(I11,t(I01))
      if (!is.null(desnformula2)) Mid = C00 - I01%*%solve(I11,t(I01)) - I02%*%solve(I22,t(I02))
      Var = solve(I00,t(solve(I00,Mid)))
    } else Var = solve(II,t(solve(II,CC))) 
    ses = sqrt((diag(Var))[1:ncol(X)])
   
    ll <- coefs + qnorm(.025)*ses; ul <- coefs + qnorm(.975)*ses
    results = round(matrix(c(coefs, exp(coefs), ses, ll, ul, exp(ll), exp(ul)),ncol=7),4)
    colnames(results) = c("Estimates", "OR", "ses", "ll", "ul", "OR.ll", "OR.ul"); rownames(results) = colnames(X)

    invisible(list(coefficients=coefs,cov=Var,fit=fit,desfit=fit1,fit2=fit2,fit3=fit3,error=error,results=results))
}

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

cover = function(thetas,ses,theta.true){
    cin=as.numeric(0,ncol(thetas))
   for (j in (1:ncol(thetas))) {
       diffs <- (thetas[,j]-theta.true[j])/ses[,j]
       cin[j] <- length(diffs[abs(diffs) < 1.96])/nrow(thetas)
    }
    cin
}