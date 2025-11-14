# dips.R
# Auxiliary functions and function for calculating DiPS

require(glmnet)
require(glmpath)
#require(np)

# AUXILIARY FUNCTIONS
sd.wgt <- function(v,w) {
  mean.wgt <- sum(v*w)/sum(w)
  var.wgt <- sum(w*(v-mean.wgt)^2)/sum(w)
  return(sqrt(var.wgt))
}

VTM <- function(vc, dm) {
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

estIPW <- function(Yi, Ti, PS, Wi=rep(1,length(Yi))) {
  ps.inv <- 1/PS
  ps.inv[is.infinite(ps.inv) & ps.inv>0] <- .Machine$double.xmax
  
  ps.invc <- 1/(1-PS)
  ps.invc[is.infinite(ps.invc) & ps.invc>0] <- .Machine$double.xmax
  
  mu1 <- sum(Yi*Ti*Wi*ps.inv)/sum(Ti*Wi*ps.inv) 
  mu2 <- sum(Wi*Yi*(1-Ti)*ps.invc)/sum(Wi*(1-Ti)*ps.invc)
  
  return(mu1-mu2)
}

#NINA EDIT weights last var as 0, so last var should be tailoring
# ADAPTIVE LASSO
Est.ALASSO.GLM = function(data,Wi=NULL,rtn="EST",ind.unp=NULL,reg=T,adap=T,BIC.factor=0.5,offset=NULL,fam0="binomial",gam=3,nu=1,tune=c("ERIC"),relax=TRUE){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x); if(is.null(Wi)){Wi=rep(1,nn)}; gc()
  if(reg){
    if(adap){
      bini = glmnet(x,y,weights=Wi,family=fam0,alpha=0,offset=offset,standardize=T)
      lam.xx=svd(t(x)%*%x)$d; tmpdf = apply(lam.xx/(lam.xx+VTM(bini$lambda,pp)),2,sum)
      bini = c(bini$beta[,which.min(deviance(bini)+2*tmpdf)]); sd.x = apply(x,2,sd)
      w.b = 1/abs(bini); w.b[bini==0] = 1/(min(abs(bini*sd.x)[abs(bini)>0])/sd.x[bini==0])
      w.b = w.b^gam
      
      #NINA EDIT
      w.b[pp] <- 0
      
    }else{
      w.b = rep(1,pp)
      
      #NINA EDIT
      w.b[pp] <- 0
    }
    lam.all = 10^seq(-10,5,0.01)
    tmpfit = glmnet(x=x, y=y,family=fam0,penalty.factor=w.b,alpha=1,lambda = lam.all,offset=offset,
                    standardize=T, weights=Wi); 
    if(tune=="BIC") {
      BIC.lam = deviance(tmpfit)+min(nn^BIC.factor,log(nn))*tmpfit$df
    } else if(tune=="ERIC") {
      if(length(unique(y))==2){
        phi = 1
      } else{
        ymat = replicate(length(lam.all),y)
        int = t(replicate(nn,tmpfit$a0))
        phi = apply((ymat - int - (x %*% tmpfit$beta))^2,2,sum)/nn
      }
      #BIC.lam = deviance(tmpfit)+2*nu*tmpfit$df*log(nn*phi/(tmpfit$lambda*nn))
      #susan edit
      BIC.lam = deviance(tmpfit)+2*nu*tmpfit$df*log(nn*phi/(tmpfit$lambda))
    }
    m.opt = which.min(BIC.lam); bhat = c(tmpfit$a0[m.opt],tmpfit$beta[,m.opt]); lamhat = tmpfit$lambda[m.opt]
    
    if(relax){
      ind.out = bhat[-1] == 0;
      tmpb = glm(y~cbind(1,x[,!ind.out])-1,family=fam0,weight=Wi)$coef
      bhat[c(T,!ind.out)] = tmpb; bhat[c(F,ind.out)]=0
    }
  } else{
    bhat = glm(y~x,weights=Wi,family=fam0,offset=offset)$coef
  }
  if(length(bhat)==1) bhat <- c(unlist(bhat),rep(0,pp))
  bhat
}

# DIPS CALCULATION
# DiPS() returns DiPS estimates given outcome (Yi), binary treatment (Ti), covariates (Xi), and weights (Wi)
# bwd: bandwidth either "plug-in" or cross-validation ("cv")
# q: kernel order
# alp: order of plug-in bandwidth h=O(n^-alph)
# C: constant multiplied for plug-in bandwidth
# constEff: whether beta0=beta1 is assumed (TRUE) or not (FALSE); if not returns list with pi_1 and pi_0
DiPS <- function(Yi,Ti,Xi,Wi=rep(1,length(Yi)),bwd="plug",q=4,alp=1/(2+q),C=1,
                 constEff=TRUE,
                 t_ind) {
  
  
  
  #Nina added these
  n <- nrow(Xi)
  p <- ncol(Xi)
  
  #NINA EDIT reorder Xi so that tailoring ind is last
  Xi <- cbind(Xi[,-(t_ind)],Xi[,t_ind])
  
  
  if(constEff) {
    betahat.ps <- Est.ALASSO.GLM(data=cbind(Ti,Xi),fam0="binomial",Wi=Wi)[1:(p+1)]
    betahat.rp <- Est.ALASSO.GLM(data=cbind(Yi,Ti,Xi),fam0="gaussian",Wi=Wi)[1:(p+2)]
    
    betahat.ps[is.na(betahat.ps)] <- 0
    betahat.rp[is.na(betahat.rp)] <- 0
    
    #NINA EDIT - keeping track of which vars selected
    selected <-  (abs(betahat.ps[2:(p+1)]) > 10^(-6)) | (abs(betahat.rp[3:(p+2)]) > 10^(-6))
    
    Si.ps <- cbind(1,Xi) %*% betahat.ps
    Si.rp <- cbind(1,Xi) %*% betahat.rp[-2]
    
    Si.ps <- as.numeric(pnorm(Si.ps, mean=mean(Si.ps), sd=sd(Si.ps)))
    Si.rp <- as.numeric(pnorm(Si.rp, mean=mean(Si.rp), sd=sd(Si.rp)))
    
    if(all(betahat.ps[-1]==0) & all(betahat.rp[-c(1,2)]==0)) {
      rv <- mean(Ti)
    } else if(all(betahat.ps[-1]==0)) {
      if(bwd=="cv") {
        b <- npregbw(Ti ~ Si.rp, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b <- C*sd.wgt(Si.rp,w=Wi)/n^alp
      } else {
        b <- bwd[2]
      }
      rv <- npreg(bws=b, txdat=Si.rp, tydat=Ti, exdat=cbind(Si.rp), ckerorder=q)$mean
    } else if(all(betahat.rp[-c(1,2)]==0)) {
      if(bwd=="cv") {
        b <- npregbw(Ti ~ Si.ps, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b <- C*sd.wgt(Si.ps,w=Wi)/n^alp
      } else {
        b <- bwd[1]
      }
      rv <- npreg(bws=b, txdat=Si.ps, tydat=Ti, exdat=cbind(Si.ps), ckerorder=q)$mean
    } else {
      if(bwd=="cv") {
        b <- npregbw(Ti ~ Si.ps + Si.rp, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b <- C*apply(cbind(Si.ps,Si.rp),2,sd.wgt,w=Wi)/n^alp
      } else {
        b <- bwd
      }
      
      rv <- npreg(bws=b, txdat=cbind(Si.ps, Si.rp), tydat=Ti, exdat=cbind(Si.ps, Si.rp), ckerorder=q)$mean
      #nina edit
      rv <- list(rv,selected)
       }
  } else {
    betahat.ps <- Est.ALASSO.GLM(data=cbind(Ti,Xi),fam0="binomial",Wi=Wi)[1:(p+1)]
    betahat.rp1 <- Est.ALASSO.GLM(data=cbind(Yi,Xi)[Ti==1,],fam0="gaussian",Wi=Wi[Ti==1])[1:(p+1)]
    betahat.rp0 <- Est.ALASSO.GLM(data=cbind(Yi,Xi)[Ti==0,],fam0="gaussian",Wi=Wi[Ti==0])[1:(p+1)]
    
    betahat.ps[is.na(betahat.ps)] <- 0
    betahat.rp1[is.na(betahat.rp1)] <- 0
    betahat.rp0[is.na(betahat.rp0)] <- 0
    
    #NINA EDIT - keeping track of which vars selected
    selected <-  (abs(betahat.ps[2:(p+1)]) > 10^(-5)) | (abs(betahat.rp1[2:(p+1)]) > 10^(-5))
    
    Si.ps <- cbind(1,Xi) %*% betahat.ps
    Si.rp1 <- cbind(1,Xi) %*% betahat.rp1
    Si.rp0 <- cbind(1,Xi) %*% betahat.rp0
    
    Si.ps <- as.numeric(pnorm(Si.ps, mean=mean(Si.ps), sd=sd(Si.ps)))
    Si.rp1 <- as.numeric(pnorm(Si.rp1, mean=mean(Si.rp1), sd=sd(Si.rp1)))
    Si.rp0 <- as.numeric(pnorm(Si.rp0, mean=mean(Si.rp0), sd=sd(Si.rp0)))
    
    spar.ps <- all(betahat.ps[-1]==0)
    spar.rp1 <- all(betahat.rp1[-c(1,2)]==0)
    spar.rp0 <- all(betahat.rp0[-c(1,2)]==0)
    
    if(spar.ps & spar.rp1 & spar.rp0) {
      rv <- rep(mean(Ti),n)
    } else if(spar.ps & spar.rp1 & !spar.rp0) {
      if(bwd=="cv") {
        b0 <- npregbw(Ti ~ Si.rp0, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b0 <- C*sd.wgt(Si.rp0,w=Wi)/n^alp
      } else {
        b0 <- bwd[2]
      }
      rv <- npreg(bws=b0, txdat=Si.rp0, tydat=Ti, exdat=cbind(Si.rp0), ckerorder=q)$mean
      rv[Ti==1] <- mean(Ti)
    } else if(spar.ps & !spar.rp1 & spar.rp0) {
      if(bwd=="cv") {
        b1 <- npregbw(Ti ~ Si.rp1, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b1 <- C*sd.wgt(Si.rp1,w=Wi)/n^alp
      } else {
        b1 <- bwd[2]
      }
      rv <- npreg(bws=b1, txdat=Si.rp1, tydat=Ti, exdat=cbind(Si.rp1), ckerorder=q)$mean
      rv[Ti==0] <- mean(Ti)
    } else if(!spar.ps & spar.rp1 & spar.rp0) {
      if(bwd=="cv") {
        b <- npregbw(Ti ~ Si.ps, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b <- C*sd.wgt(Si.ps,w=Wi)/n^alp
      } else {
        b <- bwd[1]
      }
      rv <- npreg(bws=b, txdat=Si.ps, tydat=Ti, exdat=cbind(Si.ps), ckerorder=q)$mean
    } else if(!spar.ps & !spar.rp1 & spar.rp0) {
      if(bwd=="cv") {
        b1 <- npregbw(Ti ~ Si.ps + Si.rp1, bw.method="cv.aic", ckerorder=q)$bw
        b0 <- npregbw(Ti ~ Si.ps, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b1 <- C*apply(cbind(Si.ps,Si.rp1),2,sd.wgt,w=Wi)/n^alp
        b0 <- C*apply(cbind(Si.ps,Si.rp1),2,sd.wgt,w=Wi)[1]/n^alp
      } else {
        b1 <- bwd[c(1,2)]
        b0 <- bwd[c(1)]
      }
      rv <- npreg(bws=b1, txdat=cbind(Si.ps, Si.rp1), tydat=Ti, exdat=cbind(Si.ps, Si.rp1), ckerorder=q)$mean
      rv[Ti==0] <- npreg(bws=b0, txdat=Si.ps, tydat=Ti, exdat=cbind(Si.ps[Ti==0]), ckerorder=q)$mean
    } else if(!spar.ps & spar.rp1 & !spar.rp0) {
      if(bwd=="cv") {
        b1 <- npregbw(Ti ~ Si.ps, bw.method="cv.aic", ckerorder=q)$bw
        b0 <- npregbw(Ti ~ Si.ps + Si.rp0, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b1 <- C*apply(cbind(Si.ps,Si.rp0),2,sd.wgt,w=Wi)[1]/n^alp
        b0 <- C*apply(cbind(Si.ps,Si.rp0),2,sd.wgt,w=Wi)/n^alp
      } else {
        b1 <- bwd[c(1)]
        b0 <- bwd[c(1,2)]
      }
      rv <- npreg(bws=b0, txdat=cbind(Si.ps, Si.rp0), tydat=Ti, exdat=cbind(Si.ps,Si.rp0), ckerorder=q)$mean
      rv[Ti==1] <- npreg(bws=b1, txdat=Si.ps, tydat=Ti, exdat=cbind(Si.ps[Ti==1]), ckerorder=q)$mean
    } else if(spar.ps & !spar.rp1 & !spar.rp0) {
      if(bwd=="cv") {
        b1 <- npregbw(Ti ~ Si.rp1, bw.method="cv.aic", ckerorder=q)$bw
        b0 <- npregbw(Ti ~ Si.rp0, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b1 <- C*apply(cbind(Si.rp1,Si.rp0),2,sd.wgt,w=Wi)[1]/n^alp
        b0 <- C*apply(cbind(Si.rp1,Si.rp0),2,sd.wgt,w=Wi)[2]/n^alp
      } else {
        b1 <- bwd[c(2)]
        b0 <- bwd[c(2)]
      }
      rv <- npreg(bws=b1, txdat=cbind(Si.rp1), tydat=Ti, exdat=cbind(Si.rp1), ckerorder=q)$mean
      rv[Ti==0] <- npreg(bws=b0, txdat=cbind(Si.rp0), tydat=Ti, exdat=cbind(Si.rp0[Ti==0]), ckerorder=q)$mean
    } else {
      if(bwd=="cv") {
        b1 <- npregbw(Ti ~ Si.ps + Si.rp1, bw.method="cv.aic", ckerorder=q)$bw
        b0 <- npregbw(Ti ~ Si.ps + Si.rp0, bw.method="cv.aic", ckerorder=q)$bw
      } else if(bwd=="plug") {
        b1 <- C*apply(cbind(Si.ps,Si.rp1),2,sd.wgt,w=Wi)/n^alp
        b0 <- C*apply(cbind(Si.ps,Si.rp0),2,sd.wgt,w=Wi)/n^alp
      } else {
        b1 <- c(bwd[1],bwd[2])
        b0 <- c(bwd[1],bwd[2])
      }
      
      rv1 <- npreg(bws=b1, txdat=cbind(Si.ps, Si.rp1), tydat=Ti, exdat=cbind(Si.ps, Si.rp1), ckerorder=q)$mean
      rv0 <- npreg(bws=b0, txdat=cbind(Si.ps, Si.rp0), tydat=Ti, exdat=cbind(Si.ps, Si.rp0), ckerorder=q)$mean
      
      #NINA EDIT, add in selected vars
      rv <- list(rv1,rv0,selected)
    }
  }
  
  return(rv)
}