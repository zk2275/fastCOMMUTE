### Rewrite GhostKnockoff
# Check that the input matrix is positive-definite
is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))
  
  if (p<500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    #options(warn = -1)
    lambda_min = suppressWarnings(RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values)
    options(warn = oldw)
    if( length(lambda_min)==0 ) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A)$values)
    }
  }
  return (lambda_min>tol*10)
}

create.solve_sdp_M <- function(Sigma, M=1, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  G =  Sigma
  p = dim(G)[1]
  print(G)
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),]
  Cs = c((M+1)/M*G) ##change from 2 to (M+1)/M
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p #not sure if it should be changed - may be not as it is the dimention of the linear part.
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef((M+1)/M*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) { ##change from 2 to (M+1)/M
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}


create.solve_asdp_M <- function(Sigma, M=1, max.size=500, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  
  if(ncol(Sigma) <= max.size) return(create.solve_sdp_M(Sigma, M=M, gaptol=gaptol, maxit=maxit, verbose=verbose))
  
  # Approximate the covariance matrix as block diagonal
  if(verbose) cat(sprintf("Dividing the problem into subproblems of size <= %s ... ", max.size))
  cluster_sol = divide.sdp(Sigma, max.size=max.size)
  n.blocks = max(cluster_sol$clusters)
  if(verbose) cat("done. \n")
  
  # Solve the smaller SDPs corresponding to each block
  if(verbose) cat(sprintf("Solving %s smaller SDPs ... \n", n.blocks))
  s_asdp_list = list()
  if(verbose) pb <- utils::txtProgressBar(min = 0, max = n.blocks, style = 3)
  for(k in 1:n.blocks) {
    s_asdp_list[[k]] = create.solve_sdp_M(as.matrix(cluster_sol$subSigma[[k]]), M=M, gaptol=gaptol, maxit=maxit)
    if(verbose) utils::setTxtProgressBar(pb, k)
  }
  if(verbose) cat("\n")
  
  # Assemble the solutions into one vector of length p
  p = dim(Sigma)[1]
  idx_count = rep(1, n.blocks)
  s_asdp = rep(0,p)
  for( j in 1:p ){
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j]+1
  }
  
}

GhostKnockoff.prelim <- function(cor.X, M=r, method = 'asdp', max.size = 500){
  temp.index <- 1:nrow(cor.X)
  n.X <- nrow(cor.X)
  permute.index <- rep(0, n.X)
  permute.index[-temp.index] <- 1
  Normal_50Studies <- matrix(rnorm(n.X * M * 50), n.X * M, 
                             50)
  P.each <- matrix(0, n.X, n.X)
  if (length(temp.index) != 0) {
    Sigma <- cor.X[temp.index, temp.index, drop = F]
    SigmaInv <- ginv(Sigma) ## (X'X)^-1 SparseM::Solve	or Matrix::solve
    if (method == "sdp") {
      # temp not consider about M used to be create.solve_sdp_M
      temp.s <- create.solve_sdp_M(Sigma) ## D(diag_s) is a diagonal matrix obtained by solving the convex optimization problem
    }
    if (method == "asdp") {
      temp.s <- create.solve_asdp_M(Sigma, max.size = max.size)
    }
   #  print("print out each s:")
   #  print(temp.s)
    s <- temp.s
    diag_s <- diag(s, length(s))
    if (sum(s) == 0) {
      V.left <- matrix(0, length(temp.index) * M, length(temp.index) * M)
    }
    else {
      Sigma_k <- 2 * diag_s - s * t(s * SigmaInv)
      V.each <- Matrix(forceSymmetric(Sigma_k - diag_s))
      V <- matrix(1, M, M) %x% V.each
      diag(V) <- diag(V) + rep(s, M)
      V.left <- try(t(chol(V)), silent = T)
      if (class(V.left) == "try-error") {
        svd.fit <- svd(V)
        u <- svd.fit$u
        svd.fit$d[is.na(svd.fit$d)] <- 0
        cump <- cumsum(svd.fit$d)/sum(svd.fit$d)
        n.svd <- which(cump >= 0.999)[1]
        if (is.na(n.svd)) {
          n.svd <- nrow(V)
        }
        svd.index <- intersect(1:n.svd, which(svd.fit$d != 
                                                0))
        V.left <- t(sqrt(svd.fit$d[svd.index]) * t(u[, 
                                                     svd.index, drop = F]))
      }
    } 
    #### Here we make some changes
    P.each[temp.index, temp.index] <-s * SigmaInv - diag(1, length(s))   

    V.index <- rep(temp.index, M) + rep(0:(M - 1), each = length(temp.index)) * 
      n.X
    Normal_50Studies[V.index, ] <- as.matrix(V.left %*% matrix(rnorm(ncol(V.left) * 
                                                                       50), ncol(V.left), 50))
    permute.index[temp.index[s == 0]] <- 1
  }
 
  permute.V.index <- rep(permute.index, M)
 
  P.each[permute.index == 1, ] <- 0
  print("Print out all of P.each:")
  print(P.each)
  Normal_50Studies[permute.V.index == 1, ] <- matrix(rnorm(sum(permute.index) * 
                                                             M * 50), sum(permute.index) * M, 50)
  return(list(P.each = as.matrix(P.each), V.left = V.left, 
              Normal_50Studies = as.matrix(Normal_50Studies), permute.index = permute.index, 
              M = M))
}

#### Get a smaller gamma later
GhostKnockoff.fit <- function(Zscore_0,n.study,fit.prelim,gamma=1,weight.study=NULL){
  if (length(weight.study) == 0) {
    weight.study <- sqrt(n.study)/sqrt(sum(n.study))
  }
  W.missing <- apply(!is.na(Zscore_0), 2, as.numeric)
  W <- t(weight.study * t(W.missing))
  Zscore_0[is.na(Zscore_0)] <- 0
  M <- fit.prelim$M
  n.X <- nrow(Zscore_0)
  P.each <- fit.prelim$P.each ####?

  Normal_50Studies <- fit.prelim$Normal_50Studies ####?
  for (i in 1:ncol(Zscore_0)) {
    Normal_k <- matrix(Normal_50Studies[, i], nrow = n.X)
    Zscore_i0 <- Zscore_0[, i, drop = F]
    Zscore_ik <- as.vector(P.each %*% Zscore_0[, i, drop = F]) + 
      gamma * Normal_k ### here we construct knockoff zscores
    if (i == 1) {
      GK.Zscore_0 <- W[, i] * Zscore_i0
      GK.Zscore_k <- W[, i] * Zscore_ik
    }
    else {
      GK.Zscore_0 <- GK.Zscore_0 + W[, i] * Zscore_i0
      GK.Zscore_k <- GK.Zscore_k + W[, i] * Zscore_ik
    }
  }
  T_0 <- (GK.Zscore_0)^2
  T_k <- (GK.Zscore_k)^2
  return(list(GK.Zscore_0 = GK.Zscore_0, GK.Zscore_k = GK.Zscore_k, 
              T_0 = T_0, T_k = T_k))
}





#simulation SNP data K population 
Coef.gen<- function(s, h, K, sig.delta, p, exact=TRUE, intercept){
  # add 1 for all change
  beta <- c(seq(0.4, 0.1, length.out=(s/2)), -seq(0.4, 0.1, length.out=(s/2)), rep(0, p - s))# +50 -50 1900(all0)
  w <- matrix(rep(beta, K), p, K)  #2000 3  ###h^*=0
  if(exact==TRUE){
    samp <- sample(1:p, h, replace=F)
    sign <- sample(c(-1,1), h, replace=T)
    for(k in 1:K){
      w[samp, k] <- w[samp, k] + sign*rep(as.numeric(sig.delta[k]), h)
    }
  }else{
    total_p = 1:p
    # Temporarily change the value of samp0, sign
    samp0 <- sample(total_p, h[1], replace=F) # previous h=min(h), here I set h[1]
    sign <- sample(c(-1,1), h[1], replace=T) # previous h=h, here I set h[1]  # what's the shape of h? 
    samp <- samp0 # why not samp <- ...
    w[samp, 1] <- w[samp, 1] + sign*rep(as.numeric(sig.delta), h[1])
    if(K>1){
      for(k in 2:K){
        if(h[k] > h[k-1]){
          h_diff = h[k] - h[k-1] # if it's 0, there is no 'more'?
          samp_more <- sample(total_p[-samp], h_diff, replace=F)
          sign_more <- sample(c(-1,1), h_diff, replace=T)
          w[samp, k] <- w[samp, k-1] #w <- matrix(rep(beta, K), p, K)  #2000 3
          w[samp_more, k] <- w[samp_more, k-1] + sign_more*rep(as.numeric(sig.delta), h_diff)
        }else{
          w[, k] <- w[, k-1] # but how about h[k] < h[k-1]
        }
      }
    }
  }
  
  ##why delete it??
  # beta = c(intercept[1], beta)
  # w = rbind(intercept[-1], w)
  
  return(list(beta=beta, w=w))
}

logistic<-function(x){
  1/(1+exp(-x))
}


Data.gen<- function(B, n.tar, n.src, K, n.test, SIM){
  X.tar <- y.tar <- X.all <- y.all <- c()
  X.src <- y.src <- list()
  
  #####generate target population
  SIM$reset(); id = c()
  for(i in 1:n.tar) id[i] = SIM$addUnrelatedIndividual()
  genotypes1 = SIM$gt1[id,] + SIM$gt2[id,]
  
  SIM$reset(); id = c()
  for(i in 1:n.tar) id[i] = SIM$addUnrelatedIndividual()
  genotypes2 = SIM$gt1[id,] + SIM$gt2[id,]
  
  SIM$reset(); id = c()
  for(i in 1:n.tar) id[i] = SIM$addUnrelatedIndividual()
  genotypes3 = SIM$gt1[id,] + SIM$gt2[id,]
  
  SIM$reset(); id = c()
  for(i in 1:n.tar) id[i] = SIM$addUnrelatedIndividual()
  genotypes4 = SIM$gt1[id,] + SIM$gt2[id,]
  ### change the value of y.tar
  X.tar = cbind(genotypes1, genotypes2, genotypes3, genotypes4)
  remove(genotypes1,genotypes2,genotypes3,genotypes4)
  y.tar<- rbinom(n.tar, size=1, prob=logistic(B[1,1]+X.tar%*%B[,1]))
  
  ######generate source population
  for(k in 1:K){
    SIM$reset(); id = c()
    for(i in 1:n.src[k]) id[i] = SIM$addUnrelatedIndividual()
    genotypes1 = SIM$gt1[id,] + SIM$gt2[id,]
    
    SIM$reset(); id = c()
    for(i in 1:n.src[k]) id[i] = SIM$addUnrelatedIndividual()
    genotypes2 = SIM$gt1[id,] + SIM$gt2[id,]
    
    SIM$reset(); id = c()
    for(i in 1:n.src[k]) id[i] = SIM$addUnrelatedIndividual()
    genotypes3 = SIM$gt1[id,] + SIM$gt2[id,]
    
    SIM$reset(); id = c()
    for(i in 1:n.src[k]) id[i] = SIM$addUnrelatedIndividual()
    genotypes4 = SIM$gt1[id,] + SIM$gt2[id,]
    # change the value of y.src
    X.src[[k]] = cbind(genotypes1, genotypes2, genotypes3, genotypes4)
    remove(genotypes1,genotypes2,genotypes3,genotypes4)
    y.src[[k]]<- rbinom(n.src[k], size=1, prob=logistic(X.src[[k]]%*%B[,k+1])) 
    
    X.all = rbind(X.all, X.src[[k]])
    y.all = c(y.all, y.src[[k]])
  }
  X.all = rbind(X.tar, X.all)
  y.all = c(y.tar, y.all)
  
  ##generate test data
  SIM$reset(); id = c()
  for(i in 1:n.test) id[i] = SIM$addUnrelatedIndividual()
  genotypes1 = SIM$gt1[id,] + SIM$gt2[id,]
  
  SIM$reset(); id = c()
  for(i in 1:n.test) id[i] = SIM$addUnrelatedIndividual()
  genotypes2 = SIM$gt1[id,] + SIM$gt2[id,]
  
  SIM$reset(); id = c()
  for(i in 1:n.test) id[i] = SIM$addUnrelatedIndividual()
  genotypes3 = SIM$gt1[id,] + SIM$gt2[id,]
  
  SIM$reset(); id = c()
  for(i in 1:n.test) id[i] = SIM$addUnrelatedIndividual()
  genotypes4 = SIM$gt1[id,] + SIM$gt2[id,]
  
  X.test = cbind(genotypes1, genotypes2, genotypes3, genotypes4)
  remove(genotypes1,genotypes2,genotypes3,genotypes4)
  y.test <- rbinom(n.test, size=1, prob=logistic(X.test%*%B[,1]))
  
  return(list(X.tar=X.tar, y.tar=y.tar, X.src=X.src, y.src=y.src, X.all=X.all, y.all=y.all, X.test=X.test, y.test=y.test))
}



Data.gen.one<- function(n.tar, X.pool){
  X.tar = X.pool[sample(1:nrow(X.pool), n.tar, replace = FALSE),]
  return(list(X.tar=X.tar))
}

## Make a difference
create.synthetic <- function(K, X.tar, n.src, r, B){
  n.syn = round(n.src*r)     #total number of synthetic data for each source: r*source sample
  
  X.syn <- y.syn <- list()
  for(k in 1:K){
    print(k)
    
    X.syn.k = X.tar[sample(1:nrow(X.tar), size = n.syn[k], replace=TRUE),]
    y.syn.k = rbinom(n.syn[k], 1, prob = logistic(c(cbind(1,X.syn.k) %*% B[[k]])))###?should I change the dimensions
    
    X.syn[[k]] = X.syn.k
    y.syn[[k]] = y.syn.k
  }
  return(list(X.syn=X.syn, y.syn=y.syn))
}


## Create synthetic by ghostknoff
create.synthetic.ghostknoff <- function(X.tar, y.tar, r){
  
  # compute correlation among variants
  
  cor.X<-matrix(as.numeric(corpcor::cor.shrink(X.tar)), nrow=ncol(X.tar))
  
  # fit null model, M = r
  fit.prelim <- GhostKnockoff.prelim(cor.X, M=r, method = 'sdp', max.size = 500)
  
  # hypothetical Z-scores
  Zscore_0<- 1/sqrt(nrow(X.tar))*t(X.tar)%*%y.tar
  #data size of X when only one is involve
  n.study <- nrow(X.tar)
  
  # knockoff data
  GK.stat<-GhostKnockoff.fit(Zscore_0,n.study,fit.prelim,gamma=1,weight.study=NULL)
  Zscore_syn <- GK.stat$GK.Zscore_k #df: nrow(X)*r
  
  return(Zscore_syn)
}



## Create synthetic by ghostknoff (output is just everything)
create.synthetic.ghostknoff.PLANB <- function(X.tar, y.tar, r){
  
  # compute correlation among variants

  cor.X<-matrix(as.numeric(corpcor::cor.shrink(X.tar)), nrow=ncol(X.tar))
  # fit null model, M = r
  ### temp use sdp
  fit.prelim <- GhostKnockoff.prelim(cor.X, M=r, method = 'sdp', max.size = 500)
  
  # hypothetical Z-scores
  Zscore_0<- 1/sqrt(nrow(X.tar))*t(X.tar)%*%y.tar
  #data size of X when only one is involve
  n.study <- nrow(X.tar)
  
  # knockoff data, gamma = 1
  GK.stat<-GhostKnockoff.fit(Zscore_0,n.study,fit.prelim,gamma=0.1,weight.study=NULL)
  Zscore_syn <- GK.stat$GK.Zscore_k #df: nrow(X)*r
  
  return(list(Zscore_0,Zscore_syn,GK.stat))
}



ST.init<-function(X.tar,y.tar){ 
  #single-task initialization
  p <- ncol(X.tar)
  n0.tar <- length(y.tar)
  fit <- cv.glmnet(x=X.tar, y=y.tar, nfolds = 4, family='binomial')
  lam.const = fit$lambda.min / sqrt(2*log(p)/n0.tar)
  beta0 = c(fit$glmnet.fit$a0[which(fit$lambda == fit$lambda.min)], fit$glmnet.fit$beta[, which(fit$lambda == fit$lambda.min)])
  return(list(beta0=as.numeric(beta0), lam.const=as.numeric(lam.const)))
}

TL.init<-function(X.tar, y.tar, X.src, y.src, w=NULL){
  p <- ncol(X.tar)
  n0.src <- length(y.src)
  n0.tar <- length(y.tar)
  
  if(is.null(w)){
    #source population estimates
    fit.src <- cv.glmnet(x=X.src, y=y.src, family='binomial', nfolds=4, lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar))
    lam.const = fit.src$lambda.min / sqrt(2*log(p)/n0.src)
    # temp change of w0 : delete: fit.src$glmnet.fit$a0[which(fit.src$lambda == fit.src$lambda.min)],
    w0 <- c(fit.src$glmnet.fit$a0[which(fit.src$lambda == fit.src$lambda.min)], fit.src$glmnet.fit$beta[, which(fit.src$lambda == fit.src$lambda.min)])
  }else{
    w0 = w
  }
  
  #target population estimates
  fit.tar <- cv.glmnet(x=X.tar, y=y.tar, nfolds=5, family='binomial', offset=w0[1]+X.tar%*%w0[-1], lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar))
  delta0 = c(fit.tar$glmnet.fit$a0[which(fit.tar$lambda == fit.tar$lambda.min)], fit.tar$glmnet.fit$beta[, which(fit.tar$lambda == fit.tar$lambda.min)])
  
  return(list(w0=w0, delta0=delta0, beta0=w0+delta0))
}

thres<-function(b, k, p){
  b*(rank(abs(b))>=(p-k))
}


Trans.global<-function(X.tar, y.tar, X.src, y.src, delta=NULL){
  p <- ncol(X.tar)
  n0.tar <- length(y.tar)
  K <- length(y.src)  ###(K>1)
  ##global method
  Xdelta = c()
  if(is.null(delta)){
    for(k in 1:K){
      w.k = as.numeric(ST.init(X.src[[k]], y.src[[k]])$beta0)
      fit.tar <- cv.glmnet(x=X.tar, y=y.tar, nfolds=5, family='binomial', offset=w.k[1]+X.tar%*%w.k[-1], lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar))
      delta.k = c(fit.tar$glmnet.fit$a0[which(fit.tar$lambda == fit.tar$lambda.min)], fit.tar$glmnet.fit$beta[, which(fit.tar$lambda == fit.tar$lambda.min)])
      delta.k.thre = thres(delta.k, sqrt(n0.tar), p) ###+threshold
      Xdelta.k = tcrossprod(delta.k.thre[-1], X.src[[k]])+delta.k.thre[1]
      Xdelta = c(Xdelta, Xdelta.k)
    }
  }else{
    for(k in 1:K){
      Xdelta.k = tcrossprod(delta[[k]][-1], X.src[[k]])+delta[[k]][1]
      Xdelta = c(Xdelta, Xdelta.k)
    }
  }
  
  XX.src <- yy.src <- NULL
  for(k in 1:K){
    XX.src <- rbind(XX.src, X.src[[k]])
    yy.src <- c(yy.src, y.src[[k]])
  }
  
  n0.tar.global <- length(c(y.tar,yy.src))
  offset <- c(rep(0, nrow(X.tar)), -Xdelta)
  fit.global <- cv.glmnet(x=rbind(X.tar,XX.src), y=c(y.tar,yy.src), nfolds=5, family='binomial', offset=offset, lambda=seq(0.25, 0.05, length.out=20)*sqrt(2*log(p)/n0.tar.global))
  beta.hat = c(fit.global$glmnet.fit$a0[which(fit.global$lambda == fit.global$lambda.min)], fit.global$glmnet.fit$beta[, which(fit.global$lambda == fit.global$lambda.min)])
  
  return(beta.hat)
}

Agg.fun<-function(B, X.til, y.til, const=2){
  X.til = cbind(1, X.til)
  loss.B <- apply(B, 2, function(b) - sum(y.til*log(logistic(X.til%*%b))+(1-y.til)*log(1-logistic(X.til%*%b))))
  eta.hat <- exp(-const*loss.B)/sum(exp(-const*loss.B))
  return(eta.hat)
}

Agg.fun.new<-function(B, X.til, y.til, const=2){
  X.til = cbind(1, X.til)
  XX = X.til%*%B
  eta.hat = glm(y.til~XX-1, family = binomial(link = "logit"))$coefficients
  #eta.hat = glm(y.til~XX-1, family = gaussian)$coefficients
  return(eta.hat)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

Agg.fun.Q<-function(B, X.til, y.til, total.step=10){
  p <- nrow(B)
  K <- ncol(B)
  X.til = cbind(1, X.til)
  
  XB = X.til%*%B
  colnames(B) <- NULL
  temp_delta = (rep.col(y.til, K) - XB)^2
  temp_delta <- ifelse(is.na(temp_delta), max(temp_delta), temp_delta)
  
  theta.hat <- exp(-colSums(temp_delta)/2)
  theta.hat <- ifelse(is.na(theta.hat), 0, theta.hat)
  theta.hat = theta.hat/sum(theta.hat)
  theta.old = theta.hat
  beta <- as.numeric(B%*%theta.hat)
  beta.ew <- beta
  
  Xbeta = X.til%*%beta
  for(ss in 1:total.step){
    theta.hat <- exp(-colSums(temp_delta)/2 + colSums((rep.col(Xbeta, K) - XB)^2)/8)
    if(is.na(sum(theta.hat))){
      break
    }
    theta.hat <- theta.hat/sum(theta.hat)
    beta <- as.numeric(B%*%theta.hat*1/4 + 3/4*beta)
    if(sum(abs(theta.hat-theta.old)) < 10^(-3)){
      break
    }
    theta.old = theta.hat
  }
  
  return(theta.hat)
}

mse.fun<- function(beta.true, est, X.test=NULL){
  sum((beta.true[-1]-est[-1])^2)
}

get.auc = function(X, y, beta){
  pred.y = logistic(X%*%beta)
  #mean(abs(y.test1-pred.y))
  pROC::auc(pROC::roc(y, as.numeric(pred.y)))
}