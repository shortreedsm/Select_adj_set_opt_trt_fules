#boot library is used to perform bootstrap 
library(boot)

#Gradient for the outcome 
outgradient <- function(X.out, y, Beta){
  u <- -t(X.out)%*%((y-X.out%*%Beta))/length(y)	
  
  return(c(t(u)))
}
#Gradient for the outcome (scalar - used to speed up computation)
outgradient_scalar <- function(X.out, y, Beta, j){
  u <- -t(X.out[,j])%*%((y-X.out%*%Beta))/length(y)
  
  return(u)
}

#Gradient for the treatment
trtgradient <- function(X.trt, A, beta.trt){
  v <- -t(1/(1+exp(A*X.trt%*%beta.trt)))%*%(as.vector(A)*X.trt)/length(A)
  
  return(c(t(v)))
}
#Gradient for the treatment (scalar - used to speed up computation)
trtgradient_scalar <- function(X.trt, A, Beta,j){
  v <- -t(1/(1+exp(A*X.trt%*%Beta)))%*%(as.vector(A)*X.trt[,j])/length(A)
  
  return(v)
}

#Calculates the negative sum of outcome and treatment gradients
negsumgrad <- function(X.outcome, X.treatment, Y.outcome, A.treatment, Beta.out, Beta.trt){
  u <- outgradient(X.outcome, Y.outcome, Beta.out)
  v <- trtgradient(X.treatment, A.treatment, Beta.trt)
  
  return(-1*c(u,v))
}
#Calculates the negative sum of outcome and treatment gradients (scalar - used to speed up computation)
negsumgrad_scalar <- function(X.outcome, X.treatment, Y.outcome, A.treatment, Beta.out, Beta.trt, j){
  u <- outgradient_scalar(X.outcome, Y.outcome, Beta.out,j)
  
  v <- trtgradient_scalar(X.treatment, A.treatment, Beta.trt,j)
  
  return(-1*c(u,v))
}

#"Algorithm 1" of the GMD algorithm: this function is called by the function "groupLassoAlgo"
algorithm1 <- function(Xoutcome, Xtreatment, Ytrans, Atrans, Beta.initial, Trt.intercept, Out.intercept, Alpha_A, lambda, w.k, S.true){
  
  #number of covariates
  p = (length(Beta.initial))/2
  
  #sample size
  n = length(Atrans)
  
  #used for convergence
  eps = 0.0001
  
  #used for calculation of H matrix
  X.des <- cbind(Xoutcome, Xtreatment)
  
  #H matrix
  H <- (5/4)*t(X.des)%*%X.des/n
  
  #Sub-matrix of H corresponding to group k (two covariates for each group ==> each matrix is 2x2)
  H.k <- lapply(1:(p+2), function(i) matrix(c(H[i,i], H[i,i+p+2], H[i+p+2,i], H[i+p+2,i+p+2]), 2, 2,byrow=TRUE))
  
  #Vector consisting of largest eigenvalues of H^(k) for k = 1, ... , p
  gamma.k <- lapply(H.k, function(i) max(eigen(i)$values))
  
  #initialize Beta.tilde
  B.tilde.new <- Beta.initial
  
  ##############################################################################
  ############# This loop repeats until convergence ############################
  ##############################################################################
  repeat{
    B.tilde <- B.tilde.new
    Out.intercept.old <- Out.intercept
    Trt.intercept.old <- Trt.intercept
    Alpha_A.old <- Alpha_A
    
    B.out <- B.tilde[seq(1,p,1)]
    B.trt <- B.tilde[seq(p+1,2*p,1)]
    
    #Only update the groups that satisfy the initial survival set S or fail the KKT conditions
    for(k in S.true){
      B.out <- B.tilde.new[seq(1,p,1)]
      B.trt <- B.tilde.new[seq(p+1,2*p,1)]
      
      #U is the negative gradient of the sum of loss functions at Beta_k
      U <- negsumgrad_scalar(Xoutcome, Xtreatment, Ytrans, Atrans, c(B.out, Out.intercept, Alpha_A), c(B.trt,Trt.intercept,0), k)
      
      #Update Beta.k
      Beta.k.new <- (1/gamma.k[[k]])*(U + gamma.k[[k]]*c(B.tilde.new[k],B.tilde.new[p+k]))*max(0,(
        1 - (lambda*w.k[k] / sqrt(sum((U + gamma.k[[k]]*c(B.tilde.new[k],B.tilde.new[p+k]))^2))))) 
      
      B.tilde.new[k] <- Beta.k.new[1]
      B.tilde.new[p+k] <- Beta.k.new[2]
    }
    
    #U is the negative gradient of the sum of loss functions at each intercept
    U <- negsumgrad_scalar(Xoutcome, Xtreatment, Ytrans, Atrans, c(B.out, Out.intercept, Alpha_A), c(B.trt,Trt.intercept,0), p+1)
    
    #Update intercepts
    Beta.int.new <- (1/gamma.k[[p+1]])*(U + gamma.k[[p+1]]*c(Out.intercept,Trt.intercept)) 
    
    Out.intercept <- Beta.int.new[1]
    Trt.intercept <- Beta.int.new[2]
    
    #U is negative gradient of loss function at coefficient for treatment main effect 
    U <- negsumgrad_scalar(Xoutcome, Xtreatment, Ytrans, Atrans, c(B.out, Out.intercept, Alpha_A), c(B.trt,Trt.intercept,0), p+2)[1]
    
    #Update coefficient for treatment main effect
    Alpha_A <- (1/gamma.k[[p+2]])*(U + gamma.k[[p+2]]*Alpha_A)
    
    # Convergence criteria to decide when to stop 
    B.old <- c(B.tilde, Out.intercept.old, Alpha_A.old, Trt.intercept.old)#[c(TRUE, FALSE)]
    B.new <- c(B.tilde.new, Out.intercept, Alpha_A, Trt.intercept)#[c(TRUE, FALSE)]
    
    #Determine if convergence has been reached
    if(max(abs(B.new - B.old)/(1 + abs(B.old))) < eps)#(eps/max(unlist(gamma.k))))
    {
      break
    }
  }	
  ##############################################################################
  ############################ End of loop #####################################
  ##############################################################################
  
  return(B.new)
  
} #End of algorithm 1 function

#This is "Algorithm 2" of the GMD algorithm
#Function can be used to calculate GLiDeR coefficient estimates at all values of the sequence of lambda
#If no lambda sequence is given, a uniform sequence is generated with 0.00005 as the minimum lambda value
groupLassoAlgo <- function(Xorig, Yorig, Aorig, lambda.seq=NULL){
  #number of covariates
  p <- ncol(Xorig)
  
  #sample size
  n <- length(Aorig)	
  
  #Treatment indicator is coded 1/-1 for gradient calculations 
  A <- ifelse(Aorig==1, 1, -1)
  
  # Standardized design matrix #
  Xstand <- t(t(Xorig - matrix(rep(apply(Xorig, 2, mean), nrow(Xorig)), 
                               nrow(Xorig),byrow=TRUE))*(1/apply(Xorig, 2, sd)))
  
  # Divide outcome by its standard deviation
  Yneworig <- Yorig/sd(Yorig)
  
  #New outcome model design matrix 
  Xoutorig <- cbind(Xstand, 1, Aorig)
  
  #Weights for each group: least squares if p < n-1, ridge regression with penalty parameter = 0.0001 otherwise
  if(p < (n-1)){
    w.k = sqrt(2)/abs(lm(Yneworig ~ Xoutorig - 1)$coef[1:p])
    
    #NINA EDIT - last weight is 0
    w.k[p] <- 0
    
  }else{
    tmp.mat <- diag(1, nrow=p+2)
    tmp.mat[p+1,p+1] = 0
    tmp.mat[p+2,p+2] = 0
    ridge.est <- (solve(t(Xoutorig)%*%Xoutorig + 0.0001*tmp.mat)%*%t(Xoutorig)%*%Yneworig)[1:p]	
    w.k = sqrt(2)/abs(ridge.est)
    
    #NINA EDIT - last weight is 0
    w.k[p] <- 0
  }
  
  #Another function may give columns of all zeroes for a matrix ==> sd = 0 ==> Xstand = NA; so change these values to 0, which is what we want here
  Xstand[is.na(Xstand)] <- 0
  
  #New treatment model design matrix
  Xtrtorig <- cbind(Xstand,1, 0)
  
  #Beta.lambda stores covariate coefficient values at each lambda 
  #each element in the list is a 2*p length vector of coefficients at lambda 
  #the first p elements correspond to the outcome model  
  #the second p elements correspond to the treatment model
  Beta.lambda <- list()
  
  #treatment intercept; each element in list corresponds to different lambda value
  Trt.int <- list()
  
  #outcome intercept; each element in list corresponds to different lambda value
  Out.int <- list()
  
  #treatment main effect term; each element in list corresponds to different lambda value
  alpha_A <- list()
  
  #Initial estimate of the treatment intercept
  Trt.int[[1]] <- log((sum(A == 1))/(sum(A == -1)))
  
  null.mod <- coef(lm(Yneworig ~ Aorig))
  
  #initial estimate of outcome intercept
  Out.int[[1]] <- null.mod[1]
  
  #initial estimate of treatment main effect term
  alpha_A[[1]] <- null.mod[2]
  
  #gradient of model with all covariate coefficients equal to 0
  U.null <- negsumgrad(Xoutorig, Xtrtorig, Yneworig, A, c(rep(0, p), Out.int[[1]], alpha_A[[1]]), c(rep(0, p),Trt.int[[1]], 0))
  
  #If lambda.seq is not provided (i.e. not doing cross-validation) ==> lambda.seq is given integer value
  if(length(lambda.seq) < 2){
    #We first define lambda^1, which is the smallest lambda value such that all predictors
    #(except the intercept) have zero coefficients:
    #NINA EDIT and except the tailoring
    lambda.1 <- max(unlist(lapply(1:(p-1), function(i) sqrt(sum(c(U.null[i],U.null[i+p+2])^2))/w.k[i]))) 
    
    #Define a sequence of points for lambda 
    #Currently taking 100 points between .00005 and lambda.1 
    ##(uniform sequence of lambda, which may not be ideal in some situations)
    lambda.seq <- seq(lambda.1, .00005, length.out=100)
  }
  
  #Initialize covariate coefficients
  B.init <- rep(0, 2*p)
  
  #Call algorithm1 to get solution at first (largest) value of lambda
  temp <- algorithm1(Xoutorig, Xtrtorig, Yneworig, A, B.init, Trt.int[[1]], Out.int[[1]], alpha_A[[1]], lambda.seq[1], w.k, 1:p)
  
  #covariate coefficient values
  B.tilde <- temp[1:(2*p)]
  
  #outcome intercept
  Out.int[[1]] <- temp[2*p+1]
  
  #treatment main effect term
  alpha_A[[1]] <- temp[2*p+2]
  
  #treatment intercept
  Trt.int[[1]] <- temp[2*p+3]
  
  Beta.lambda[[1]] <- B.tilde
  
  #Check which groups need to be included in algorithm 1 and find a solution at each lambda in the
  #sequence by calling algorithm 2 
  ###################################################################################################	 
  for(iter in 1:(length(lambda.seq)-1)){
    
    B.tilde.current <- Beta.lambda[[iter]]
    
    #outcome model covariate coefficient estimates at current value of lambda 
    B.out.lambda <- B.tilde.current[seq(1,p,1)]
    
    #treatment model covariate coefficient estimates at current value of lambda 
    B.trt.lambda <- B.tilde.current[seq(p+1,2*p,1)]
    
    #U.lambda is the negative gradient of the sum of loss functions at lambda^[l]
    U.lambda <- negsumgrad(Xoutorig, Xtrtorig, Yneworig, A, 
                           c(B.out.lambda, Out.int[[iter]],alpha_A[[iter]]), c(B.trt.lambda,Trt.int[[iter]], 0)) 
    
    #S denotes the set of groups that satisfy the strong rule
    ##lambda.new = lambda^[l+1], lambda = lambda^[l]
    S <- lapply(1:p, function(i) sqrt(sum(c(U.lambda[i],U.lambda[p+i+2])^2)) >= 
                  w.k[i]*(2*lambda.seq[iter+1] - lambda.seq[iter]))
    
    indicator <- unlist(S)
    
    #if true then solution is found
    if(sum(indicator) == 0){
      Beta.lambda[[iter+1]] = B.tilde.current
      Trt.int[[iter+1]] = Trt.int[[iter]]
      Out.int[[iter+1]] = Out.int[[iter]]
      alpha_A[[iter+1]] = alpha_A[[iter]]
    } else {
      
      ##################REPEAT UNTIL KKT CONDITIONS ARE SATISFIED############################################
      repeat{
        S <- indicator	
        
        #Outcome and treatment covariates that satisfy strong rule
        Xtempout <- cbind(t(t(Xoutorig[,1:p])*S), 1, Aorig)
        Xtemptrt <- cbind(t(t(Xtrtorig[,1:p])*S), 1, 0)
        
        #call algorithm1
        temp <- algorithm1(Xtempout, Xtemptrt, Yneworig, A, Beta.lambda[[iter]], 
                           Trt.int[[iter]], Out.int[[iter]], alpha_A[[iter]], lambda.seq[iter+1], w.k, which(S==TRUE))
        
        Out.int[[iter+1]] <- temp[2*p+1]
        
        alpha_A[[iter+1]] <- temp[2*p+2]
        
        Trt.int[[iter+1]] <- temp[2*p+3]
        
        #Next check whether Beta.S satisfies the KKT conditions
        B.out.S <- temp[seq(1,p,1)]
        
        B.trt.S <- temp[seq(p+1,2*p,1)]
        
        #gradient of the sum of loss functions at Beta.S
        U.lambda.new <- negsumgrad(Xoutorig, Xtrtorig, Yneworig, A, 
                                   c(B.out.S, Out.int[[iter+1]], alpha_A[[iter+1]]), c(B.trt.S, Trt.int[[iter+1]], 0))
        
        #For each group k NOT in S, if the follow inequality holds, 
        #then Beta.S is the desired solution at lambda^[l+1]:
        #else, we add the group that failed the KKT condition to the survival set S
        KKT <- lapply(1:p, function(i) sqrt(sum(c(U.lambda.new[i], U.lambda.new[i+p+2])^2)) > 
                        w.k[[i]]*lambda.seq[iter+1])
        
        if(sum(unlist(S) == 0) == sum((unlist(KKT)+unlist(S)) == 0)){
          Beta.lambda[[iter+1]] = c(B.out.S, B.trt.S)
          break
        } else {
          indicator = (unlist(KKT) + unlist(S)) != 0 
        }
      } #End of repeat-loop	
      
    } #End of if-else that checks if S is all FALSE		
    
    
  } #End of for-loop
  
  return(list(Beta.lambda, Out.int, alpha_A, Trt.int, lambda.seq, w.k[1:p]))
  
}	#End of function

#Function calculates GLiDeR point estimates
#This function requires a matrix of covariates (with each row representing one observation): Xorig 
#a vector of outcomes: Yorig
#and a vector of binary treatment indicators (coded 1/0)
#This function chooses lambda via Generalized cross validation (GCV) and returns (at the selected lambda)
#the estimated average causal treatment effect, the outcome model covariate coefficients, 
#the treatment model covariate coefficients, the outcome model intercept, the treatment main effect
#term in the outcome model, the treatment model intercept, and the selected lambda value
GLiDeR <- function(Xorig,Yorig, Aorig, lambda=NULL){
  
  #number of covariates
  p<-ncol(Xorig)
  
  #call groupLassoAlgo to obtain coefficient estimates at each value lambda 	
  temp <- groupLassoAlgo(Xorig, Yorig, Aorig, lambda)
  
  #Covariate coefficient estimates
  Beta.lambda <- temp[[1]]
  
  #outcome model intercepts 
  Beta.out.int <- temp[[2]]
  
  #outcome model treatment main effect estimates
  Beta.a <- temp[[3]]
  
  #treatment model intercepts
  Beta.trt.int <- temp[[4]]
  
  
  #Standardized covariates used in model
  Xstand <- t(t(Xorig - matrix(rep(apply(Xorig, 2, mean), nrow(Xorig)), 
                               nrow(Xorig),byrow=TRUE))*(1/apply(Xorig, 2, sd)))
  
  #Standardized outcome
  Ytrans <- Yorig/sd(Yorig)
  
  lambda <- temp[[5]]
  
  #number of lambda
  n.lambda <- length(lambda)
  
  #Vector listing all coefficients of the outcome model
  outCoef <- unlist(lapply(1:n.lambda, function(i) Beta.lambda[[i]][seq(1,p,1)]))
  
  #Vector listing all coefficients of treatment model
  trtCoef <- unlist(lapply(1:n.lambda, function(i) Beta.lambda[[i]][seq(p+1,2*p,1)])) 
  
  #Matrix of outcome coefficients: row i represents the vector of coef at lambda i
  outCoefMat <- matrix(outCoef, length(Beta.lambda), p, byrow=TRUE) 
  
  #Matrix of trtment coefficients: row i represents the vector of coef at lambda i
  trtCoefMat <- matrix(trtCoef, length(Beta.lambda), p, byrow=TRUE)
  
  #sample size
  n <- length(Aorig)	
  
  #estimated probabilities of receiving treatment
  Phat <- lapply(1:n.lambda, function(i)
    exp(rep(Beta.trt.int[[i]], n) + Xstand%*%trtCoefMat[i,])/(1 + exp(rep(Beta.trt.int[[i]], n) + Xstand%*%trtCoefMat[i,]))
  )
  
  #estimated outcomes if treatment not received (if A = 0)
  Yhat0 <- lapply(1:n.lambda, function(i)
    (rep(Beta.out.int[[i]], n) + Xstand%*%outCoefMat[i,])*sd(Yorig) 
  )
  
  #estimated outcomes if treatment is received (if A = 1)
  Yhat1 <- lapply(1:n.lambda, function(i)
    (rep(Beta.out.int[[i]], n) + Beta.a[[i]] + Xstand%*%outCoefMat[i,])*sd(Yorig) 
  )
  
  muhat1 <- lapply(1:n.lambda, function(i)
    (Aorig*Yorig/Phat[[i]]) - ((Aorig - Phat[[i]])/Phat[[i]])*Yhat1[[i]]
  )
  
  muhat0 <- lapply(1:n.lambda, function(i)
    ((1-Aorig)*Yorig)/(1-Phat[[i]]) + ((Aorig-Phat[[i]])/(1-Phat[[i]]))*Yhat0[[i]]
  )
  
  #estimated average causal treatment effects
  drest <- lapply(1:n.lambda, function(i)
    (1/n)*(sum(muhat1[[i]]) - sum(muhat0[[i]]))
  )
  
  #predicted outcomes: used for GCV	
  Y.pred <- lapply(1:n.lambda, function(i)
    ifelse(Aorig == 1, Yhat1[[i]], Yhat0[[i]])
  )
  
  #residual sum of squares: used for GCV	
  RSS <- lapply(1:n.lambda, function(i)
    sum((Y.pred[[i]] - Yorig)^2)
  )
  
  #used in denominator of GCV statistic	
  val <- lapply(1:n.lambda, function(i)
    2 + sum(ifelse(outCoefMat[i,] == 0, 0, 1)) + sum(abs(outCoefMat[i,])/(sqrt(2)/temp[[6]]))
  )
  
  #GCV statistic at each value lambda	
  GCV <- lapply(1:n.lambda, function(i)
    RSS[[i]]/(1 - ((val[[i]])/n))^2
  )
  
  #return the following that correspond to the minimum GCV statistic
  return(list(delta=drest[[which.min(unlist(GCV))]], 
              alpha = outCoefMat[which.min(unlist(GCV)),], 
              gamma = trtCoefMat[which.min(unlist(GCV)),], 
              alpha0 = Beta.out.int[[which.min(unlist(GCV))]], 
              gamma0 = Beta.trt.int[[which.min(unlist(GCV))]], 
              alpha_A = Beta.a[[which.min(unlist(GCV))]], 
              lambda_star = lambda[which.min(unlist(GCV))]))
  
} #End function

#This function is used to obtain bootstrap samples of the average causal treatment effect of GLiDeR
#The function uses the library "boot" and is called by "GLiDeR_bootstrap" (below)
groupLassoDRboot <- function(X,d){
  sz <- ncol(X)
  Xorig <- X[d,1:(sz-3)]
  p<-ncol(Xorig)
  Yorig <- X[d,sz-2]
  Aorig <- X[d,sz-1]
  lambda <- X[1,sz]
  
  temp <- groupLassoAlgo(Xorig, Yorig, Aorig, c(lambda+.001, lambda))
  
  Beta.lambda <- temp[[1]][[2]]
  Beta.trt.int <- temp[[4]][[2]]
  Beta.out.int <- temp[[2]][[2]]
  Beta.a <- temp[[3]][[2]]
  
  #Vector listing all coefficients of the outcome model
  outCoef <- Beta.lambda[seq(1,p,1)]
  
  #Vector listing all coefficients of treatment model
  trtCoef <- Beta.lambda[seq(p+1,2*p,1)]
  
  n <- length(Aorig)
  
  #Standardized covariates used in model
  Xstand <- t(t(Xorig - matrix(rep(apply(Xorig, 2, mean), 
                                   nrow(Xorig)), nrow(Xorig),byrow=TRUE))*(1/apply(Xorig, 2, sd)))
  
  Phat <- exp(rep(Beta.trt.int, nrow(Xstand)) + Xstand%*%trtCoef)/
    (1 + exp(rep(Beta.trt.int, nrow(Xstand)) + Xstand%*%trtCoef))
  
  Yhat0 <- (rep(Beta.out.int, n) + Xstand%*%outCoef)*sd(Yorig) 
  
  
  Yhat1 <- (rep(Beta.out.int, n) + Beta.a + Xstand%*%outCoef)*sd(Yorig) 
  
  
  muhat1 <- (Aorig*Yorig/Phat) - ((Aorig - Phat)/Phat)*Yhat1
  
  
  muhat0 <- ((1-Aorig)*Yorig)/(1-Phat) + ((Aorig-Phat)/(1-Phat))*Yhat0
  
  drest <- mean(muhat1 - muhat0)
  
  return(drest)
  
} #End function

#Function calculates GLiDeR bootstrap samples of the average causal treatment effect
#This function requires a matrix of covariates (with each row representing one observation): Xorig 
#a vector of outcomes: Yorig
#a vector of binary treatment indicators (coded 1/0)
#and the lambda value selected by GCV (found using "GLiDeR" function)
GLiDeR_bootstrap <- function(Xorig, Yorig, Aorig, lambda_star, B){
  boot_samp <- boot(cbind(Xorig,Yorig,Aorig,rep(lambda_star,nrow(Xorig))),groupLassoDRboot,R=B)$t
  return(boot_samp)
}