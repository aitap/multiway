linesearch <- function(OLD, Xa,CkrB, A, B, C,nfac, iter) {
  SSE <-0
  
  if(iter>2){
    Aex <- OLD$Amat2+iter^(1/3)*(OLD$Amat1-OLD$Amat2)
    Bex <- OLD$Bmat2+iter^(1/3)*(OLD$Bmat1-OLD$Bmat2)
    Cex <- OLD$Cmat2+iter^(1/3)*(OLD$Cmat1-OLD$Cmat2)
    for(u in 1:nfac) CkrB[,u] <- kronecker(Cex[,u], Bex[,u])
    SSE <- sum((Xa - tcrossprod(Aex,CkrB))^2)
  }
  OLD$Amat2 <- OLD$Amat1
  OLD$Bmat2 <- OLD$Bmat1
  OLD$Cmat2 <- OLD$Cmat1 
  if(SSE<OLD$SSEold){
    OLD$Amat1<-Aex
    OLD$Bmat1<-Bex
    OLD$Cmat1<-Cex
    OLD$SSEold <- SSE
  } else{
    OLD$Amat1<-A
    OLD$Bmat1<-B
    OLD$Cmat1<-C
  }
  return(OLD)
}