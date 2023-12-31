#.libPaths("/n/home02/tiangu/apps/R_4_1")
library(glmnet)
library(lassoshooting)
library(mvtnorm)
library(pROC)
library(rmutil)
library(grplasso)
library(Matrix)
library(sim1000G)

source('fastCOMMUTE_functions.R')

id <- 1# set 1 first to try. id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sim_setting <- expand.grid(sim = 1:100, delta = 0.8, h = c("(10,10;10)","(30,30;30)","(60,60;60)","(10,60;60)"), r = 25)

sim = sim_setting[id, 'sim']
sig.delta = sim_setting[id, 'delta'] 
h.tmp = sim_setting[id, 'h'] 
h = c(as.numeric(gsub(".*\\((.*)\\,.*", "\\1", h.tmp)), as.numeric(gsub(".*\\,(.*)\\;.*", "\\1", h.tmp)), as.numeric(gsub(".*\\;(.*)\\).*", "\\1", h.tmp)))
r = sim_setting[id, 'r']
dir_out = 'path_to_save_result'

exact = FALSE             #TRUE=S1; FALSE=S2
K = 3                     #total number of source population 
p = 200                  #number of variables
s = p/2
n0 = 100                  #target population size
nt = c(1000, 1500, 2000)    #source population size
n.test = 1000
intercept = rep(1, K+1) 

simu = function(sim, p, s, K, n0, nt, h, sig.delta, intercept, exact, r, n.test){
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
  
  
  # calculate the running time
  start.time <- Sys.time()
  ###fastCOMMUTE M=0
  Z_score1 <- 1/sqrt(nrow(X.src[[1]]))*t(X.src[[1]])%*%y.src[[1]]
  Z_score2 <- 1/sqrt(nrow(X.src[[2]]))*t(X.src[[2]])%*%y.src[[2]]
  Z_score3 <- 1/sqrt(nrow(X.src[[3]]))*t(X.src[[3]])%*%y.src[[3]]
  beta.fastCOMMUTE.M0 = ginv(49.81987*t(X.tar)%*%X.tar)%*%( # cannot use solve
    t(X.tar)%*%y.tar + 
      sqrt(nrow(X.src[[1]]))*Z_score1+10.98082*t(X.tar)%*%X.tar%*%delta.TL[[1]]+
      sqrt(nrow(X.src[[2]]))*Z_score2+16.18668*t(X.tar)%*%X.tar%*%delta.TL[[2]]+
      sqrt(nrow(X.src[[3]]))*Z_score3+21.65237*t(X.tar)%*%X.tar%*%delta.TL[[3]]
  )
  
  
  ###fastCOMMUTE M=1
  
  beta.fastCOMMUTE.M1 = ginv(38.83905*t(X.tar)%*%X.tar+t(X.src[[1]])%*%X.src[[1]])%*%(
    t(X.tar)%*%y.tar + 
      t(X.src[[1]])%*%y.src[[1]]+t(X.src[[1]])%*%X.src[[1]]%*%delta.TL[[1]]+
      sqrt(nrow(X.src[[2]]))*Z_score2+16.18668*t(X.tar)%*%X.tar%*%delta.TL[[2]]+
      sqrt(nrow(X.src[[3]]))*Z_score3+21.65237*t(X.tar)%*%X.tar%*%delta.TL[[3]]
  )
  
  ###fastCOMMUTE M=2
  
  beta.fastCOMMUTE.M2 <- 
    ginv(t(X.tar)%*%X.tar+t(X.src[[1]])%*%X.src[[1]]+
           t(X.src[[2]])%*%X.src[[2]]+
           21.65237*t(X.tar)%*%X.tar)%*%(
             t(X.tar)%*%y.tar + 
               t(X.src[[1]])%*%y.src[[1]]+t(X.src[[1]])%*%X.src[[1]]%*%delta.TL[[1]]+
               t(X.src[[2]])%*%y.src[[2]]+t(X.src[[2]])%*%X.src[[2]]%*%delta.TL[[2]]+
               sqrt(nrow(X.src[[3]]))*Z_score3+21.65237*t(X.tar)%*%X.tar%*%delta.TL[[3]]
           )
  
  ###fastCOMMUTE M=3
  
  beta.fastCOMMUTE.M3 <- 
    ginv(t(X.tar)%*%X.tar+t(X.src[[1]])%*%X.src[[1]]+
           t(X.src[[2]])%*%X.src[[2]]+
           t(X.src[[3]])%*%X.src[[3]]  )%*%(
             t(X.tar)%*%y.tar + 
               t(X.src[[1]])%*%y.src[[1]]+t(X.src[[1]])%*%X.src[[1]]%*%delta.TL[[1]]+
               t(X.src[[2]])%*%y.src[[2]]+t(X.src[[2]])%*%X.src[[2]]%*%delta.TL[[2]]+
               t(X.src[[3]])%*%y.src[[3]]+t(X.src[[3]])%*%X.src[[3]]%*%delta.TL[[3]]
           )
  
  
  
  # plot
  par(mfrow = c(2, 2))
  plot(beta.fastCOMMUTE.M0)
  plot(beta.fastCOMMUTE.M1)
  plot(beta.fastCOMMUTE.M2,ylim = c(-1,1))
  plot(beta.fastCOMMUTE.M3)
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print("How long do we spend time on running the fastCOMMUTE code?:")
  print(time.taken)
  
  
  
  methods = c('beta.tar', 'w1', 'w2', 'w3', 
              'beta.fastCOMMUTE.M0','beta.fastCOMMUTE.M1',
              'beta.fastCOMMUTE.M2','beta.fastCOMMUTE.M3')
  mseout = c()
  for(i in 1:length(methods)){
    mseout = c(mseout, mse.fun(beta.true, as.numeric(get(methods[i]))))
    
  }
  
  return(list(cbind(methods, mseout)))
}


out = tryCatch(simu(sim, p, s, K, n0, nt, h, sig.delta, intercept, exact, r, n.test), error=function(err) list(matrix(NA,8,2)))
save(out, file = paste0(dir_out, 'h=', h.tmp, '_delta=', sig.delta, '_exact=', exact, '_ntimes=', r, '_2size_', sim,'.Rdata'))



