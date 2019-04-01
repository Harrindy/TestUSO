#' A TestUSO Function
#'
#' This function allows you to test against USO based on X and Y via the sample-based-critical-value approach
#' @param X a vector of the sample of size m from the distribution F
#' @param Y a vector of the sample of size n from the distribution G
#' @param alpha significance level. Default is 0.05
#' @param B number of bootstrapped samples. Default is 1000
#' @param L number of Monte Carolo samples to compute sample-based critical values. Default is 1000
#' @param graph plot or not. Default is FALSE
#' @keywords Goodness-of-fit test for USO
#' @export
#' @return Test statistics, critical values, and conclusion; i.e, reject (if 1) or do not reject (if 0)
#' @examples X=rnorm(30,0,1);Y=rnorm(40,1,1);
#' GoF4USO(X,Y,alpha=0.05,graph=TRUE,L=1000,B=1000)

GoF4USO=function(X,Y,alpha=0.05,graph=FALSE,L=1000,B=1000){
  require("OrdMonReg")
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
  alpha.a = array(0,n+1)
  u_slope[1:n] = (1-Rmn_data[1:n])/(1-u[1:n])
  alpha.a = cummin(u_slope)
  LSMRmn = 1-(1-u)*alpha.a
  l1=sum( (LSMRmn[1:n]-2*Rmn_data[1:n]+LSMRmn[1:n]+alpha.a[1:n]*1/n)/n/2 )
  l2=sqrt(sum(((alpha.a[1:n]*1/n)^2/3+
                 (LSMRmn[1:n]-Rmn_data[1:n])*(LSMRmn[1:n]+
                                                alpha.a[1:n]*1/n-Rmn_data[1:n]))/n))
  lsup=max(LSMRmn[1:n]+alpha.a[1:n]/n-Rmn_data[1:n])
  Mmn_l1=l1*sqrt(m*n/(m+n))
  Mmn_l2=l2*sqrt(m*n/(m+n))
  Mmn_lsup=lsup*sqrt(m*n/(m+n))
  
  u_slope=c()
  u_slope= (1-Rmn_data[1:n])/(1-u[1:n])
  Iso_r=BoundedAntiMean(u_slope[2:n],w=rep(1/n,n-1),a=rep(0,n-1),b=rep(1,n-1))
  Iso_R=c(1-(1-u[1:n])*c(1,Iso_r),1)
  
  if(graph)
  {
    #plot(u,Rmn_data,type='s',xlim=c(0,1),ylim=c(0,1),col="black",lwd=1,ylab=expression(paste(R[mn](u),", ",MR[mn](u))),cex.lab=1)
    plot(u,Rmn_data,type='s',xlim=c(0,1),ylim=c(0,1),col="black",lwd=1,ylab="",cex.lab=1)
    lines(u,u,lty=3)
    for(kk in 1:n)
    {
      polygon(c((kk-1)/n,(kk-1)/n,kk/n,kk/n),
              c(LSMRmn[kk],Rmn_data[kk],Rmn_data[kk],LSMRmn[kk]+alpha.a[kk]*1/n),
              col="yellow",border=NA)
      segments((kk-1)/n,LSMRmn[kk],kk/n,LSMRmn[kk]+alpha.a[kk]*1/n,col="blue")
    }
    lines(u,Iso_R,lwd=2,col="red",lty=2)
  }
  
  ind=1:head(which(Iso_R==1),1)
  sample.y=u[ind]
  sample.x=Iso_R[ind]
  inv_iso=approxfun(sample.x,sample.y,method="linear")
  
  U=ecdf(X)(X)
  V=ecdf(Y)(Y)
  
  lsmrate.b=c()
  CK=ceiling(B*alpha)
  for(b in 1:CK)
  {
    X.b=sample(U,m,replace=TRUE)
    Y.b=sample(V,n,replace=TRUE)
    LSMRmn_data.b=LSM(X.b,Y.b)$M
    lsmrate.b=rbind(lsmrate.b,(1-LSMRmn_data.b[1:n])/(1-u[1:n]))
  }
  for(b in (CK+1):B)
  {
    lsmrate.b=apply(lsmrate.b,2,sort)
    temp.b=lsmrate.b[1,]
    X.b=sample(U,m,replace=TRUE)
    Y.b=sample(V,n,replace=TRUE)
    LSMRmn_data.b=LSM(X.b,Y.b)$M
    comp.b=(1-LSMRmn_data.b[1:n])/(1-u[1:n])
    ind=comp.b>temp.b
    lsmrate.b[1,ind]=comp.b[ind]
  }
  lsmrate.b=apply(lsmrate.b,2,sort,decreasing=TRUE)
  alpha.K=tail(which(lsmrate.b[,n]==1),1)
  
  lsmrate.b=c()
  CK=ceiling(B*alpha)
  for(b in 1:CK)
  {
    X.b=sample(X,m,replace=TRUE)
    Y.b=sample(Y,n,replace=TRUE)
    LSMRmn_data.b=LSM(X.b,Y.b)$M
    lsmrate.b=rbind(lsmrate.b,(1-LSMRmn_data.b[1:n])/(1-u[1:n]))
  }
  for(b in (CK+1):B)
  {
    lsmrate.b=apply(lsmrate.b,2,sort)
    temp.b=lsmrate.b[1,]
    X.b=sample(X,m,replace=TRUE)
    Y.b=sample(Y,n,replace=TRUE)
    LSMRmn_data.b=LSM(X.b,Y.b)$M
    comp.b=(1-LSMRmn_data.b[1:n])/(1-u[1:n])
    ind=comp.b>temp.b
    lsmrate.b[1,ind]=comp.b[ind]
  }
  
  lsmrate.b=apply(lsmrate.b,2,sort,decreasing=TRUE)
  fix.rate=lsmrate.b[alpha.K,]
  fix.rate=sort(fix.rate,decreasing = TRUE)
  lb=c(1-fix.rate*(1-(0:(n-1))/n),1)
  ind=1:head(which(lb==1),1)
  sample.y=u[ind]
  sample.x=lb[ind]
  sample.r=fix.rate[1:(length(ind)-1)]
  inv_rs=approxfun(sample.x,sample.y,method="linear")
  
  if(graph==TRUE)
  {
    lines(u,lb,lwd=2,col="green",lty=4)
    legend("bottomright",c(expression(R[mn]),
                           expression(MR[mn]),
                           expression(R[0]),
                           expression(hat(R)[iso]),
                           expression(hat(R)[bs])),
           lty=c(1,1,3,2,4),lwd=c(1,1,1,2,2),col=c("black","blue","black","red","green"))
  }
  
  TS.lf=c()
  TS.iso=c()
  TS.rs=c()
  for(b in 1:L)
  {
    U.B=runif(m)
    X.iso=inv_iso(U.B)
    X.rs=inv_rs(U.B)
    Y.B=runif(n)
    LSMRmn_data.lf=LSM(U.B,Y.B,ts=TRUE)
    LSMRmn_data.iso=LSM(X.iso,Y.B,ts=TRUE)
    LSMRmn_data.rs=LSM(X.rs,Y.B,ts=TRUE)
    TS.lf=rbind(TS.lf,c(LSMRmn_data.lf$Mmn_l1,
                        LSMRmn_data.lf$Mmn_l2,
                        LSMRmn_data.lf$Mmn_lsup))
    TS.iso=rbind(TS.iso,c(LSMRmn_data.iso$Mmn_l1,
                          LSMRmn_data.iso$Mmn_l2,
                          LSMRmn_data.iso$Mmn_lsup))
    TS.rs=rbind(TS.rs,c(LSMRmn_data.rs$Mmn_l1,
                        LSMRmn_data.rs$Mmn_l2,
                        LSMRmn_data.rs$Mmn_lsup))
  }
  lf_cv=apply(TS.lf,2,quantile,prob=1-alpha)
  iso_cv=apply(TS.iso,2,quantile,prob=1-alpha)
  bs_cv=apply(TS.rs,2,quantile,prob=1-alpha)
  
  if(alpha==0.01){Fixed_cv=c(0.751, 0.860, 1.623)
  }else if(alpha==0.05){Fixed_cv=c(0.580, 0.676, 1.353)
  }else if(alpha==0.1){Fixed_cv=c(0.496, 0.586, 1.219)
  }else{print("please use significance level 0.01, 0.05, or 0.1 for the Tang et al. (2017) approach")}
  
  Test_statistic=c(Mmn_l1,Mmn_l2,Mmn_lsup)
  Reject_USO_bs=1*c(Test_statistic>bs_cv)
  Reject_USO_iso=1*c(Test_statistic>iso_cv)
  Reject_USO_lf=1*c(Test_statistic>lf_cv)
  Reject_USO_fix=1*c(Test_statistic>Fixed_cv)
  
  res=cbind(Test_statistic, 
            Fixed_cv,Reject_USO_fix,
            lf_cv,Reject_USO_lf,
            iso_cv,Reject_USO_iso,
            bs_cv, Reject_USO_bs)
  res=data.frame(res)
  rownames(res)=c("p=1","p=2","p=infinity")
  print("1: reject USO; 0: do not rejct USO")
  print("Fixed_cv: critical values based on the limiting distribution at the least favorable configuration")
  print("Reject_USO_fix: reject or not based on the limiting distribution at the least favorable configuration")
  print("lf_cv: critical values based on the finite distribution at the least favorable configuration")
  print("Reject_USO_lf: reject or not based on the finite distribution at the least favorable configuration")
  print("iso_cv: critical values based on the isotonic fitting")
  print("Reject_USO_iso: reject or not based on the isotonic fitting")
  print("bs_cv: critical values based on the resample method")
  print("Reject_USO_bs: reject or not based on the resample method")
  #print(paste("the tuning parameter is ",alpha.K,sep=""))
  return(res)
}
