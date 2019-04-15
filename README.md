# TestUSO

This R package can be used to conduct nonparametric goodness-of-fit tests for uniform stochastic ordering proposed by Tang et al. (2017) and Wang et al. (2019).

# Installation

    library(devtools)
    install_github("Harrindy/TestUSO",force=TRUE) 

# Illustration:

    library(TestUSO)
    
    set.seed(100)
    #Generate data:
    X=rnorm(30,0,1);Y=rnorm(40,1,1)  

    #Calculate the sample ODC Rmn:
    EstODC(X,Y,graph=TRUE)    
    
   ![Optional Text](../master/Rmn.png)
   
    #Calculate the least star-shaped majorant of the sample ODC Rmn:
    LSM(X,Y,graph=TRUE) 
   
   ![Optional Text](../master/MRmn.png)
    
    #Conduct the test via the approaches proposed by Tang et al. (2017) and Wang et al. (2019):
    GoF4USO(X,Y,alpha=0.05,graph=TRUE) 
    
   ![Optional Text](../master/NewMethods.png)
    
    # [1] "Fixed_cv: critical values using Tang et al. (2017)"
    # [1] "Reject_USO_fix: reject (1) or not (0) using Tang et al. (2017)"
    # [1] "lf_cv: critical values using the finie-sample version of Tang et al. (2017)"
    # [1] "Reject_USO_lf: reject (1) or not (0) using the finie-sample version of Tang et al. (2017)"
    # [1] "iso_cv: critical values using method one of Wang et al. (2019)"
    # [1] "Reject_USO_iso: reject (1) or not (0) using method one of Wang et al. (2019)"
    # [1] "bs_cv: critical values using method two of Wang et al. (2019)"
    # [1] "Reject_USO_bs: reject (1) or not (0) using method two of Wang et al. (2019)"
    #               Test_statistic Fixed_cv Reject_USO_fix      lf_cv Reject_USO_lf    iso_cv Reject_USO_iso     bs_cv Reject_USO_bs
    #p=1            0.03968089        0.580              0  0.6137245             0 0.1069128              0 0.1879498             0
    #p=2            0.05587717        0.676              0  0.6911214             0 0.1687542              0 0.2891967             0
    #p=infinity     0.20701967        1.353              0  1.3118565             0 0.5672693              0 0.8625819             0

# References:

Tang, C., Wang, D., and Tebbs, J. (2017). Nonparametric goodness-of-fit tests for uniform stochastic ordering. Annals of Statistic 48, 2565-2589.

Wang, D., Tang, C., and Tebbs, J. (2019+). More powerful goodness-of-fit tests for uniform stochastic ordering. Submitted for publication.

# Author:
Dewei Wang (deweiwang@stat.sc.edu)
