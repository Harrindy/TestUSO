#' A TestUSO Function
#'
#' This function allows you to compute the empirical ODC.
#' @param X a vector of the sample of size m from the distribution F
#' @param Y a vector of the sample of size n from the distribution G
#' @param graph plot the empirical ODC or not? Default=FALSE
#' @keywords Empirical ODC
#' @export
#' @return The value of empirical ODC at points 0, 1/n, 2/n, ..., (n-1)/n,1
#' @examples X=rnorm(30,0,1);Y=rnorm(40,1,1);
#' EstODC(X,Y,graph=TRUE)

EstODC=function(X,Y,graph=FALSE){
  X=sort(X)
  Y=sort(Y)
  m = length(X)
  n = length(Y)
  Rmn = c(array(0,(n)),1)
  for (i in 2:(n+1)){
    Rmn[i] = sum(X<=Y[i-1])/m
  }
  if(graph)
  {
    u=seq(0,1,length=n+1)
    plot(u,Rmn,type='s',xlim=c(0,1),ylim=c(0,1),col="blue")
    points(u,Rmn,pch=20,cex=0.5)
    lines(u,u,lty=3)
  }
  return(Rmn)
}