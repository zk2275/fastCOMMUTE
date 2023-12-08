beta_fastCOMMUTE_generation = function(K, method = "commute",X.tar,X.src,y.tar,y.src,delta) {
  
  nZ_score <- Z_score <- list()
  XTX_tar <- t(X.tar)%*%X.tar
  XTY_tar <- t(X.tar)%*%y.tar
  for (i in 1:K) {
    nZ_score[[i]] <- t(X.src[[i]])%*%y.src[[i]]
    Z_score[[i]] <- 1/sqrt(nrow(X.src[[i]]))*t(X.src[[i]])%*%y.src[[i]]
    
  }
  
  
  if(method == "commute"){
    scale_src_tar <- XTX_src <- XTY_src <- XTX_src_delta <- 
      scale_XTY_tar<- XTX_tar_delta <- scale_XTX_tar_delta<- list()
    for (i in 1:K) {
      XTX_src[[i]] <- t(X.src[[i]])%*%X.src[[i]]
      XTY_src[[i]] <- t(X.src[[i]])%*%y.src[[i]]
      scale_src_tar[[i]] <- (nrow(X.src[[i]])-1)/(nrow(X.tar)-1)
      XTX_src_delta[[i]] <- t(X.src[[i]])%*%X.src[[i]]%*%delta.TL[[i]]
      XTX_tar_delta[[i]] <- t(X.tar)%*%X.tar%*%delta.TL[[i]] # may be unnecessary
      
      scale_XTX_tar_delta[[i]] <- scale_src_tar[[i]]*t(X.tar)%*%X.tar%*%delta.TL[[i]]
      scale_XTY_tar[[i]]  <- scale_src_tar[[i]]*t(X.tar)%*%y.tar
    }
    
    ##beta_commute0
    sum0 <- list()
    for (i in 1:K) {
      sum0[[i]] <- nZ_score[[i]] + XTX_src_delta[[i]] #1129 scale_src_tar[[i]]*XTX_tar_delta
    }
    
    #sum_of_list <- Reduce(`+`, my_list)
    beta_commute0 <- ginv((1+Reduce(`+`, scale_src_tar)  )*XTX_tar  )%*%
      ( XTY_tar +  Reduce(`+`, sum0) )
    
    ##beta_commute1-K
    beta_commute <- list()
    for (k in 1:K) {
      
      # calibrated
      sum_calibrated <- Reduce(`+`, XTY_src[1:k])+
        Reduce(`+`, XTX_src_delta[1:k])
      
      if(k==K){
        #combine
        beta_commute[[k]] <- ginv(XTX_tar + 
                                    Reduce(`+`, XTX_src[1:k]) )%*%
          ( XTY_tar + sum_calibrated)
      }
      else{ 
        # syn
        sum_syn <- Reduce(`+`, scale_XTY_tar[(k+1):K])+ # 1129 nz_score
          Reduce(`+`, scale_XTX_tar_delta[(k+1):K]) # 1129 1129 scale_src_tar[[i]]*XTX_tar_delta
        #combine
        beta_commute[[k]] <- ginv((1+Reduce(`+`, scale_src_tar[(k+1):K]) )*XTX_tar + 
                                    Reduce(`+`, XTX_src[1:k]) )%*%
          (XTY_tar + sum_calibrated + sum_syn)
      }
      #+Reduce(`+`, scale[(k+1):K]) )*XTX_tar
      
    }
    ### !!!put k=0 to the last one
    beta_commute[[K+1]]<-beta_commute0
    return (beta_commute)
  }
  
  #fastCOMMUTE
}


data_generation <- function(sim, p, s, K, n0, nt, h, sig.delta, intercept, exact, r, n.test){
  
  
  set.seed(sim)
  n.tar <- n0 
  n.src <- nt 
  coef.all <- Coef.gen(s=s, h=h, K=K, sig.delta, p=p, exact, intercept) #exact = False?
  beta.true <- as.numeric(coef.all$beta)
  B=cbind(beta.true, coef.all$w) #for target then for source(w for source)
  
  X.pool = data.matrix(read.table('./Xpool.txt', header = F))
  # Choose 200 variables for simulation temporarily
  X.pool = X.pool[,1:p]
  # There are 20000 rows, just sample n.tar(100) rows for target
  sample.tar = sample(1:nrow(X.pool), n.tar, replace = FALSE)
  # Sample out
  X.tar = X.pool[sample.tar,] 
  # include intercepts
  X.tar = cbind(rep(1,nrow(X.tar)),X.tar)
  # Generate y.tar using MLR, without intercept, don't forget the epsilon. 
  y.tar<- X.tar %*% beta.true + rnorm(nrow(X.tar))
  
  # temp generate values for the targeted response variable
  # y.tar.fast.M3 <- X.tar%*%beta.true
  
  # renew X.pool
  X.pool = X.pool[-sample.tar,] 
  X.src <- y.src <- list()
  for(k in 1:K){
    # from the rest of X.pool, get all of the source data
    sample.k = sample(1:nrow(X.pool), n.src[k], replace = FALSE) 
    X.src[[k]] = X.pool[sample(1:nrow(X.pool), n.src[k], replace = FALSE),] 
    # include intercepts
    X.src[[k]] = cbind(rep(1,nrow(X.src[[k]])),X.src[[k]])
    y.src[[k]]<- X.src[[k]]%*%B[, k+1]+ rnorm(nrow(X.src[[k]]))
    X.pool = X.pool[-sample.k,]
  }
  
  # Test data
  X.test = X.pool[sample(1:nrow(X.pool), n.test, replace = FALSE),]
  # include intercepts
  X.test = cbind(rep(1,nrow(X.test)),X.test)
  y.test <- X.test%*%B[,1]+rnorm(n.test)
  
  
  ###directly use cv.glmnet to fit lasso
  beta.tar <- as.numeric(ST.init(X.tar, y.tar)$beta0)
  
  w = list()
  TL = list()
  for(k in 1:K){
    w[[k]] <- as.numeric(ST.init(X.src[[k]], y.src[[k]])$beta0)
    print(paste0('w',k))
    TL[[k]] <- TL.init(X.tar, y.tar, X.src=X.src[[k]], y.src=y.src[[k]], w=w[[k]])
    print(paste0('beta.TL',k))
  }
  
  ###calculate delta for each source
  delta.TL = list()
  for(k in 1:K){
    delta.TL[[k]] = thres(TL[[k]]$delta0, n.tar, p) ###add threshold
  }
  
  
  
  #source-only LASSO + threshold
  w1 = thres(w[[1]], sqrt(n.src[1]), p)
  w2 = thres(w[[2]], sqrt(n.src[2]), p)
  w3 = thres(w[[3]], sqrt(n.src[3]), p)
  
  return(list(X.tar=X.tar,X.src=X.src,
              y.tar=y.tar,y.src=y.src,
              delta.TL=delta.TL,beta.true=beta.true))
}