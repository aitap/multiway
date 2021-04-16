linesearch <- function(OLD, Xa,CkrB, A, B, C,nfac, iter) {
  SSE <-0
  
  if(iter>2){
    Aex <- OLD$Amat2+iter^(1/OLD$k)*(OLD$Amat1-OLD$Amat2)
    Bex <- OLD$Bmat2+iter^(1/OLD$k)*(OLD$Bmat1-OLD$Bmat2)
    Cex <- OLD$Cmat2+iter^(1/OLD$k)*(OLD$Cmat1-OLD$Cmat2)
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
    OLD$k<-OLD$k*0.9
  } else{
    OLD$Amat1<-A
    OLD$Bmat1<-B
    OLD$Cmat1<-C
    OLD$k<-OLD$k*1.01
    
  }
  OLD$k.log <- c(OLD$k.log, OLD$k)
  print(OLD$k)
  return(OLD)
}