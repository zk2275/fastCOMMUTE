beta_fastCOMMUTE_generation = function(K, method = "commute",X.tar,X.src,y.tar,y.src,delta) {
  
  nZ_score <- Z_score <- list()
  XTX_tar <- t(X.tar)%*%X.tar
  XTY_tar <- t(X.tar)%*%y.tar
  for (i in 1:K) {
    nZ_score[[i]] <- t(X.src[[i]])%*%y.src[[i]]
    Z_score[[i]] <- 1/sqrt(nrow(X.src[[i]]))*t(X.src[[i]])%*%y.src[[i]]
    
  }
  
  
  if(method == "commute"){
     coef <- XTX_src <- XTY_src <- XTX_src_delta <- XTX_tar_delta <- coef_XTX_tar_delta<- list()
     for (i in 1:K) {
       XTX_src[[i]] <- t(X.src[[i]])%*%X.src[[i]]
       XTY_src[[i]] <- t(X.src[[i]])%*%y.src[[i]]
       coef[[i]] <- sum(XTX_src[[i]])/sum(XTX_tar)
       XTX_src_delta[[i]] <- t(X.src[[i]])%*%X.src[[i]]%*%delta.TL[[i]]
       XTX_tar_delta[[i]] <- t(X.tar)%*%X.tar%*%delta.TL[[i]] # may be unnecessary
       coef_XTX_tar_delta[[i]] <- coef[[i]]*t(X.tar)%*%X.tar%*%delta.TL[[i]]
     }
     
     #beta_commute0
     sum0 <- list()
     for (i in 1:K) {
       sum0[[i]] <- nZ_score[[i]] + coef[[i]]*XTX_tar_delta[[i]]
     }
     
     #sum_of_list <- Reduce(`+`, my_list)
     beta_commute0 <- ginv((1+Reduce(`+`, coef)  )*XTX_tar  )%*%
                      ( XTY_tar +  Reduce(`+`, sum0) )
     
     #beta_commute1-K
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
           sum_syn <- Reduce(`+`,nZ_score[(k+1):K])+
           Reduce(`+`, coef_XTX_tar_delta[(k+1):K]) 
           #combine
           beta_commute[[k]] <- ginv((1+Reduce(`+`, coef[(k+1):K]) )*XTX_tar + 
                                   Reduce(`+`, XTX_src[1:k]) )%*%
                                (XTY_tar + sum_calibrated + sum_syn)
         }
 
      
     }
     ### !!!put k=0 to the last one
     beta_commute[[K+1]]<-beta_commute0
     return (beta_commute)
  }
 
#fastCOMMUTE
}