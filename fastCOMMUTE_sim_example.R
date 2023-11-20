#.libPaths("/n/home02/tiangu/apps/R_4_1")
library(glmnet)
library(lassoshooting)
library(mvtnorm)
library(pROC)
library(rmutil)
library(grplasso)
library(Matrix)
library(sim1000G)

source('fastCOMMUTE_functions_h_identical.R')

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
p = 200 #temp pres
s = 100
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
  # choose 200 for simulation temporarily
  X.pool = X.pool[,1:200]
  
  sample.tar = sample(1:nrow(X.pool), n.tar, replace = FALSE) #there are 20000 rows, just sample n.tar(100) rows
  X.tar = X.pool[sample.tar,] # sample out
  # Temporily change the value of y.tar, y.src and y.test (no intercept)
  # y.tar cannot work in the former function, because of the dimensions of B(2000*1)
  y.tar<- rbinom(n.tar, size=1, prob=logistic(X.tar%*%as.matrix(B[, 1]))) # 100*2000 for X.tar 
  
  # temp generate values for the targeted response variable
  y.tar.fast.M3 <- X.tar%*%beta.true
  
  X.pool = X.pool[-sample.tar,] # renew X.pool
  X.src <- y.src <- list()
  for(k in 1:K){
    sample.k = sample(1:nrow(X.pool), n.src[k], replace = FALSE) # from the rest of X.pool, get all of the source data
    X.src[[k]] = X.pool[sample(1:nrow(X.pool), n.src[k], replace = FALSE),] 
    y.src[[k]]<- rbinom(n.src[k], size=1, prob=logistic(X.src[[k]]%*%B[, k+1]))
    X.pool = X.pool[-sample.k,]
  }
  
  X.test = X.pool[sample(1:nrow(X.pool), n.test, replace = FALSE),]
  y.test <- rbinom(n.test, size=1, prob=logistic(X.test%*%B[,1]))
  X.test = cbind(1, X.test) # plus 1, why?
  
  ###directly use glmnet to fit lasso
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
  ww = list(w1=w1, w2=w2, w3=w3)
  
  methods = c('beta.tar', 'w1', 'w2', 'w3')
  mseout = c(); aucout = c()
  for(i in 1:length(methods)){
    mseout = c(mseout, mse.fun(beta.true, as.numeric(get(methods[i]))))
    aucout = c(aucout, get.auc(X.test, y.test, get(methods[i])))
  }
  cbind(methods, mseout, aucout)
  
  # ###SURE Screening ?
  auc.tar = c()
  for(k in 1:K){
    auc.tmp = get.auc(cbind(1, X.tar), y.tar, ww[[k]])
    auc.tar = c(auc.tar, auc.tmp)
  }
  K.1st = which(auc.tar==sort(auc.tar, decreasing = T)[1])
  K.2nd = which(auc.tar==sort(auc.tar, decreasing = T)[2])
  K.3rd = which(auc.tar==sort(auc.tar, decreasing = T)[3])
  
  
  ##{Test Code
  ###fastCOMMUTE M=0
  # calculate the running time
  start.time <- Sys.time()
  
  print("let's generate our first beta_fastCOMMUTE!")
  Z.score.M0 <- 1/sqrt(nrow(X.tar))*t(X.tar)%*%y.tar.fast.M0  # dimension =p*1
  Z.score.syn.M0 <- create.synthetic.ghostknoff(X.tar, y.tar.fast.M0, r) 
  # dimension =p*1, suming up r columns with p*1's dimension
  r.Z.score.syn.M0 <-  colSums(Z.score.syn.M0)
  beta.fastCOMMUTE.M0 =sqrt(n.tar) /(1+3*r) * ginv(t(X.tar)%*%X.tar)%*%((Z.score.M0+3*r.Z.score.syn.M0))
  # temp change: give a value 1 to the vector to keep the dimension same
  # beta.fastCOMMUTE.M0 =c(1,beta.fastCOMMUTE.M0)
  
  # plot
  par(mfrow = c(2, 2))
  plot(beta.fastCOMMUTE.M0)
  
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print("How long do we spend time on running the fastCOMMUTE code?:")
  print(time.taken)
  
  ### MLR
  beta.MLR.M0 = ginv(t(X.tar)%*%X.tar)%*%(t(X.tar)%*%y.tar.fast.M0)
  plot(beta.MLR.M0)
  
  plot(beta.true)
  
  ### just simple double
  Z.score.simple <- 1/sqrt(nrow(X.tar))*t(X.tar)%*%y.tar.fast.M0  # dimension =p*1
  Z.score.syn.simple <- create.synthetic.ghostknoff(X.tar, y.tar.fast.M0, r=1) 
  beta.fastCOMMUTE.simple = ginv(t(X.tar)%*%X.tar)%*%((sqrt(nrow(X.tar))*Z.score.syn.simple))
  plot(beta.fastCOMMUTE.simple)
  ##}
  
  ###singleTL
  beta.TL1 = TL[[1]]$beta0
  beta.TL2 = TL[[2]]$beta0
  beta.TL3 = TL[[3]]$beta0
  
  ###COMMUTE M=0
  data.syn <- create.synthetic(K, X.tar, n.src, r, B=list(beta.TL1, beta.TL2, beta.TL3))
  beta.syn12.M0 = ST.init(rbind(X.tar,data.syn$X.syn[[K.1st]],data.syn$X.syn[[K.2nd]]), c(y.tar,data.syn$y.syn[[K.1st]],data.syn$y.syn[[K.2nd]]))$beta0
  beta.syn123.M0 = ST.init(rbind(X.tar,data.syn$X.syn[[1]],data.syn$X.syn[[2]],data.syn$X.syn[[3]]), c(y.tar,data.syn$y.syn[[1]],data.syn$y.syn[[2]],data.syn$y.syn[[3]]))$beta0
  
  ###COMMUTE M=1
  beta.syn12.M1 = Trans.global(rbind(X.tar,data.syn$X.syn[[K.2nd]]), c(y.tar,data.syn$y.syn[[K.2nd]]), X.src=list(X.src[[1]]), y.src=list(y.src[[1]]), delta=list(delta.TL[[1]]))
  beta.syn123.M1 = Trans.global(rbind(X.tar,data.syn$X.syn[[2]],data.syn$X.syn[[3]]), c(y.tar,data.syn$y.syn[[2]],data.syn$y.syn[[3]]), X.src=list(X.src[[1]]), y.src=list(y.src[[1]]), delta=list(delta.TL[[1]]))
  
  ###COMMUTE M=2
  beta.syn12.M2 = Trans.global(X.tar, y.tar, X.src=list(X.src[[1]],X.src[[2]]), y.src=list(y.src[[1]], y.src[[2]]), delta=list(delta.TL[[1]], delta.TL[[2]]))
  beta.syn123.M2 = Trans.global(rbind(X.tar,data.syn$X.syn[[3]]), c(y.tar,data.syn$y.syn[[3]]), X.src=list(X.src[[1]], X.src[[2]]), y.src=list(y.src[[1]],y.src[[2]]), delta=list(delta.TL[[1]],delta.TL[[2]]))
  
  ###COMMUTE M=3 (pooledTL)
  beta.TL12 = Trans.global(X.tar, y.tar, X.src=list(X.src[[K.1st]],X.src[[K.2nd]]), y.src=list(y.src[[K.1st]], y.src[[K.2nd]]), delta=list(delta.TL[[K.1st]], delta.TL[[K.2nd]]))
  beta.TL123 = Trans.global(X.tar, y.tar, X.src=list(X.src[[1]],X.src[[2]],X.src[[3]]), y.src=list(y.src[[1]], y.src[[2]], y.src[[3]]), delta=list(delta.TL[[1]], delta.TL[[2]], delta.TL[[3]]))
    
  
  
  
  #####################
  ### aggregation
  #####################
  # Temporarily change the value of y.til
  X.til = Data.gen.one(n.tar=100, X.pool)$X.tar# But what does til mean? 
  y.til = rbinom(100, size=1, prob=logistic(X.til%*%beta.true))
  
  ### naive aggregation
  B.naive = matrix(0, nrow = p+1, ncol = K)
  for (k in 1:K) {
    B.naive[,k] = ww[[k]]
  }
  B.naive = cbind(beta.tar, B.naive)
  #ETA<- exp(-const*loss.B)/sum(exp(-const*loss.B))
  wt.naive.agg <- Agg.fun(B.naive, X.til, y.til, const=1) 
  wt.naive.agg.new <- Agg.fun.new(B.naive, X.til, y.til, const=1)
  beta.naive.agg <- B.naive%*%wt.naive.agg
  beta.naive.agg.new <- B.naive%*%wt.naive.agg.new
  
  ### single-source TL
  B.singleTL = cbind(beta.tar, beta.TL1, beta.TL2, beta.TL3)
  wt.singleTL.agg <- Agg.fun(B.singleTL, X.til, y.til, const=1)
  wt.singleTL.agg.new <- Agg.fun.new(B.singleTL, X.til, y.til, const=1)
  beta.singleTL.agg <- B.singleTL%*%wt.singleTL.agg
  beta.singleTL.agg.new <- B.singleTL%*%wt.singleTL.agg.new
  
  ### COMMUTE M=0 (federated)
  B.syn.M0 = cbind(beta.tar, beta.TL1, beta.TL2, beta.TL3, beta.syn12.M0, beta.syn123.M0)
  wt.syn.M0.agg <- Agg.fun(B.syn.M0, X.til, y.til, const=1)
  wt.syn.M0.agg.new <- Agg.fun.new(B.syn.M0, X.til, y.til, const=1)
  beta.syn.M0.agg <- B.syn.M0%*%wt.syn.M0.agg
  beta.syn.M0.agg.new <- B.syn.M0%*%wt.syn.M0.agg.new
  
  ### COMMUTE M=1
  B.syn.M1 = cbind(beta.tar, beta.TL1, beta.TL2, beta.TL3, beta.syn12.M1, beta.syn123.M1)
  wt.syn.M1.agg <- Agg.fun(B.syn.M1, X.til, y.til, const=1)
  wt.syn.M1.agg.new <- Agg.fun.new(B.syn.M1, X.til, y.til, const=1)
  beta.syn.M1.agg <- B.syn.M1%*%wt.syn.M1.agg
  beta.syn.M1.agg.new <- B.syn.M1%*%wt.syn.M1.agg.new
  
  ### COMMUTE M=2
  B.syn.M2 = cbind(beta.tar, beta.TL1, beta.TL2, beta.TL3, beta.syn12.M2, beta.syn123.M2)
  wt.syn.M2.agg <- Agg.fun(B.syn.M2, X.til, y.til, const=1)
  wt.syn.M2.agg.new <- Agg.fun.new(B.syn.M2, X.til, y.til, const=1)
  beta.syn.M2.agg <- B.syn.M2%*%wt.syn.M2.agg
  beta.syn.M2.agg.new <- B.syn.M2%*%wt.syn.M2.agg.new
  
  ### COMMUTE M=3 (pooledTL)
  B.TL = cbind(beta.tar, beta.TL1, beta.TL2, beta.TL3, beta.TL12, beta.TL123)
  wt.TL.agg <- Agg.fun(B.TL, X.til, y.til, const=1)
  wt.TL.agg.new <- Agg.fun.new(B.TL, X.til, y.til, const=1)
  beta.syn.M3.agg <- B.TL%*%wt.TL.agg
  beta.syn.M3.agg.new <- B.TL%*%wt.TL.agg.new
  
  ### fastCOMMUTE M=0 MLR
  # Using R square
  print("Let's work on the rating of performance!")
  
  
  #B.fastCOMMUTE.TL = cbind(beta.tar, beta.TL1, beta.TL2, beta.TL3, beta.fastCOMMUTE.M0)
  #wt.fastCOMMUTE.TL.agg <- Agg.fun(B.fastCOMMUTE.TL, X.til, y.til, const=1)
  #wt.fastCOMMUTE.TL.agg.new <- Agg.fun.new(B.fastCOMMUTE.TL, X.til, y.til, const=1)
  #beta.fast.M3.agg <- B.fastCOMMUTE.TL%*%wt.fastCOMMUTE.TL.agg
  #beta.fast.M3.agg.new <- B.fastCOMMUTE.TL%*%wt.fastCOMMUTE.TL.agg.new
  
  
  #evaluation using test data
  methods = c('beta.tar', 'w1', 'w2', 'w3', 
              'beta.naive.agg','beta.singleTL.agg','beta.syn.M0.agg','beta.syn.M1.agg','beta.syn.M2.agg','beta.syn.M3.agg',
              'beta.naive.agg.new','beta.singleTL.agg.new','beta.syn.M0.agg.new','beta.syn.M1.agg.new','beta.syn.M2.agg.new','beta.syn.M3.agg.new')
  mseout = c(); aucout = c()
  for(i in 1:length(methods)){
    mseout = c(mseout, mse.fun(beta.true, as.numeric(get(methods[i]))))
    aucout = c(aucout, get.auc(X.test, y.test, get(methods[i])))
  }
  cbind(methods, mseout, aucout)
  
  print("For the fastCOMMUTE part")
  
  return(list(cbind(methods, mseout, aucout)))
}


out = tryCatch(simu(sim, p, s, K, n0, nt, h, sig.delta, intercept, exact, r, n.test), error=function(err) list(matrix(NA,16,3)))
save(out, file = paste0(dir_out, 'h=', h.tmp, '_delta=', sig.delta, '_exact=', exact, '_ntimes=', r, '_2size_', sim,'.Rdata'))



