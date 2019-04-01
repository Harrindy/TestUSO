#' A TestUSO Function
#'
#' This function allows you to compute the least star-shape majorant of the empirical ODC.
#' @param X a vector of the sample of size m from the distribution F
#' @param Y a vector of the sample of size n from the distribution G
#' @param graph plot the empirical ODC, the LSM of the empirical ODC, and the equal distribution line. Default=FALSE
#' @param ts caclcualte the test statistics or not? default=FALSE. If true, return the test statistics when p=1, p=2, p=infinity
#' @keywords Least star-shaped majorant
#' @export
#' @return The value of the least star-shaped majorant of the empirical ODC at points 0, 1/n, 2/n, ..., (n-1)/n,1
#' @examples X=rnorm(30,0,1);Y=rnorm(40,1,1);
#' LSM(X,Y,graph=TRUE,ts=TRUE)

LSM=function(X,Y,graph=FALSE,ts=FALSE){
  X=sort(X)
  Y=sort(Y)
  m = length(X)
  n = length(Y)
  Rmn_data = c(array(0,(n)),1)
  for (i in 2:(n+1)){
    Rmn_data[i] = sum(X<=Y[i-1])/m
  }
  u = seq(0,1,by=1/n)
  u_slope = array(0,n+1)
  LSMRmn = array(1,n+1)
  alpha = array(0,n+1)
  u_slope[1:n] = (1-Rmn_data[1:n])/(1-u[1:n])
  alpha = cummin(u_slope)
  LSMRmn = 1-(1-u)*alpha
  if(graph)
  {
    plot(u,Rmn_data,type='s',xlim=c(0,1),ylim=c(0,1),col="black",lwd=1,ylab=expression(paste(R[mn](u),", ",MR[mn](u))),cex.lab=1)
    #points(u,Rmn_data,pch=20,cex=1)
    #points(u,LSMRmn,pch=2,cex=0.5)
    lines(u,u,lty=5)
    for(kk in 1:n)
    {
      polygon(c((kk-1)/n,(kk-1)/n,kk/n,kk/n),
              c(LSMRmn[kk],Rmn_data[kk],Rmn_data[kk],LSMRmn[kk]+alpha[kk]*1/n),
              col="yellow",border=NA)
      segments((kk-1)/n,LSMRmn[kk],kk/n,LSMRmn[kk]+alpha[kk]*1/n,col="blue")
    }
    legend("bottomright",c("Equal distribution line","Empirical ODC","Least star-shaped majorant "),lty=c(5,1,1),lwd=c(1,1,1),col=c("black","blue","black"))
  }
  if(ts)
  {
    l1=sum( (LSMRmn[1:n]-2*Rmn_data[1:n]+LSMRmn[1:n]+alpha[1:n]*1/n)/n/2 )
    l2=sqrt(sum(((alpha[1:n]*1/n)^2/3+
                   (LSMRmn[1:n]-Rmn_data[1:n])*(LSMRmn[1:n]+
                                                  alpha[1:n]*1/n-Rmn_data[1:n]))/n))
    lsup=max(LSMRmn[1:n]+alpha[1:n]/n-Rmn_data[1:n])
    return(list(M=LSMRmn,Mmn_l1=l1*sqrt(m*n/(m+n)),Mmn_l2=l2*sqrt(m*n/(m+n)),Mmn_lsup=lsup*sqrt(m*n/(m+n))))
  }
  return(list(M=LSMRmn))
}