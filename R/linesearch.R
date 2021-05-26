
linesearch1 <- function(OLD,X, Xa,CkrB, A, B, C,nfac, iter) {
  SSE <-0
  
  if(iter>2){
    #одномерная оптимизация
    #htfkmy lfyyst
    #графики ошибок    
    r <- iter^(1/OLD$k)
    
    Aex <- OLD$Amat2 + r * (OLD$Amat1 - OLD$Amat2)
    Bex <- OLD$Bmat2 + r * (OLD$Bmat1 - OLD$Bmat2)
    Cex <- OLD$Cmat2 + r * (OLD$Cmat1 - OLD$Cmat2)
    for(u in 1:nfac) CkrB[,u] <- kronecker(Cex[,u], Bex[,u])
    SSE <- sum((Xa - tcrossprod(Aex,CkrB))^2)
  }
   ## error measure
  for(u in 1:nfac) CkrB[,u] <- kronecker(C[,u], B[,u])
  Xhat <- tcrossprod(A, CkrB)
  SSEout <- sum((Xa - Xhat)^2)
   ##
  
  OLD$Amat2 <- OLD$Amat1
  OLD$Bmat2 <- OLD$Bmat1
  OLD$Cmat2 <- OLD$Cmat1 
 
  if(SSE<SSEout&&iter>2){
    OLD$Amat1<-Aex
    OLD$Bmat1<-Bex
    OLD$Cmat1<-Cex
    #OLD$k<-OLD$k*0.9
  } else{
    OLD$Amat1<-A
    OLD$Bmat1<-B
    OLD$Cmat1<-C
    #OLD$k<-OLD$k*1.01
    
  }
  ##
  OLD$dA.log <-c(OLD$dA.log,sumsq(OLD$Amat1 - OLD$Amat2))
  OLD$dB.log <-c(OLD$dB.log,sumsq(OLD$Bmat1 - OLD$Bmat2))
  OLD$dC.log <-c(OLD$dC.log,sumsq(OLD$Cmat1 - OLD$Cmat2))
  ##
  OLD$k.log <- c(OLD$k.log, iter^(1/OLD$k))
  OLD$SSE.log <- c(OLD$SSE.log,SSE)
  OLD$SSEout.log <- c(OLD$SSEout.log,SSEout)
  return(OLD)
}

linesearch2 <- function(OLD,X, Xa,CkrB, A, B, C,nfac, iter) {
  SSE <-0
  
  if(iter>2){
    #одномерная оптимизация
    #htfkmy lfyyst
    #графики ошибок    
    r <- iter^(1/2)
    
    Aex <- OLD$Amat2 + r * (OLD$Amat1 - OLD$Amat2)
    Bex <- OLD$Bmat2 + r * (OLD$Bmat1 - OLD$Bmat2)
    Cex <- OLD$Cmat2 + r * (OLD$Cmat1 - OLD$Cmat2)
    for(u in 1:nfac) CkrB[,u] <- kronecker(Cex[,u], Bex[,u])
    SSE <- sum((Xa - tcrossprod(Aex,CkrB))^2)
  }
  ## error measure
  for(u in 1:nfac) CkrB[,u] <- kronecker(C[,u], B[,u])
  Xhat <- tcrossprod(A, CkrB)
  SSEout <- sum((Xa - Xhat)^2)
  ##
  
  OLD$Amat2 <- OLD$Amat1
  OLD$Bmat2 <- OLD$Bmat1
  OLD$Cmat2 <- OLD$Cmat1 
  
  if(SSE<SSEout&&iter>2){
    OLD$Amat1<-Aex
    OLD$Bmat1<-Bex
    OLD$Cmat1<-Cex
    #OLD$k<-OLD$k*0.9
  } else{
    OLD$Amat1<-A
    OLD$Bmat1<-B
    OLD$Cmat1<-C
    #OLD$k<-OLD$k*1.01
    
  }
  ##
  OLD$dA.log <-c(OLD$dA.log,sumsq(OLD$Amat1 - OLD$Amat2))
  OLD$dB.log <-c(OLD$dB.log,sumsq(OLD$Bmat1 - OLD$Bmat2))
  OLD$dC.log <-c(OLD$dC.log,sumsq(OLD$Cmat1 - OLD$Cmat2))
  ##
  OLD$k.log <- c(OLD$k.log, iter^(1/2))
  OLD$SSE.log <- c(OLD$SSE.log,SSE)
  OLD$SSEout.log <- c(OLD$SSEout.log,SSEout)
  return(OLD)
}

###################################################################
ELS <- function(OLD,X, Xa,CkrB, A, B, C,nfac, iter) {
  SSE <-0
  
  if(iter>2){

    ferr <-function(R){
      Aex <- OLD$Amat2 + R * (OLD$Amat1 - OLD$Amat2)
      Bex <- OLD$Bmat2 + R * (OLD$Bmat1 - OLD$Bmat2)
      Cex <- OLD$Cmat2 + R * (OLD$Cmat1 - OLD$Cmat2)
      for(u in 1:nfac) CkrB[,u] <- kronecker(Cex[,u], Bex[,u])
      return(sum((Xa - tcrossprod(Aex,CkrB))^2))
    }
    r <- optimise(ferr,c(0,20))$minimum
    OLD$k<-r
    
    Aex <- OLD$Amat2 + r * (OLD$Amat1 - OLD$Amat2)
    Bex <- OLD$Bmat2 + r * (OLD$Bmat1 - OLD$Bmat2)
    Cex <- OLD$Cmat2 + r * (OLD$Cmat1 - OLD$Cmat2)
    for(u in 1:nfac) CkrB[,u] <- kronecker(Cex[,u], Bex[,u])
    SSE <- sum((Xa - tcrossprod(Aex,CkrB))^2)
  }
  ## error measure
  for(u in 1:nfac) CkrB[,u] <- kronecker(C[,u], B[,u])
  Xhat <- tcrossprod(A, CkrB)
  SSEout <- sum((Xa - Xhat)^2)
  ##
  
  OLD$Amat2 <- OLD$Amat1
  OLD$Bmat2 <- OLD$Bmat1
  OLD$Cmat2 <- OLD$Cmat1 
  
  if(SSE<SSEout&&iter>2){
    OLD$Amat1<-Aex
    OLD$Bmat1<-Bex
    OLD$Cmat1<-Cex
  } else{
    OLD$Amat1<-A
    OLD$Bmat1<-B
    OLD$Cmat1<-C
  }
  ##
  OLD$dA.log <-c(OLD$dA.log,sumsq(OLD$Amat1 - OLD$Amat2))
  OLD$dB.log <-c(OLD$dB.log,sumsq(OLD$Bmat1 - OLD$Bmat2))
  OLD$dC.log <-c(OLD$dC.log,sumsq(OLD$Cmat1 - OLD$Cmat2))
  ##
  OLD$k.log <- c(OLD$k.log, OLD$k)
  OLD$SSE.log <- c(OLD$SSE.log,SSE)
  OLD$SSEout.log <- c(OLD$SSEout.log,SSEout)
  return(OLD)
}