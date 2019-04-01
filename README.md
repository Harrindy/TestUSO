# TestUSO

This R package can be used to conduct nonparametric goodness-of-fit tests for uniform stochastic ordering proposed by Tang et al. (2017) and Wang et al. (2018+).

# Installation

    library(devtools)
    install_github("Harrindy/TestUSO") 

# Illustration:

    library(TestUSO)
    
    #Generate data:
    X=rnorm(30,0,1);Y=rnorm(40,1,1)  

    #Calculate the sample ODC Rmn:
    EstODC(X,Y,graph=FALSE)    

    #Calculate the least star-shaped majorant of the sample ODC Rmn:
    LSM(X,Y) 

    #Conduct the test via the approaches proposed by Tang et al. (2017) and Wang et al. (2018+):
    GoF4USO(X,Y,alpha=0.05,graph=TRUE) 

# References:

Tang, C., Wang, D., and Tebbs, J. (2017). Nonparametric goodness-of-fit tests for uniform stochastic ordering. Annals of Statistic 48, 2565-2589.

Wang, D., Tang, C., and Tebbs, J. (2018+). A more powerful goodness-of-fit test for uniform stochastic ordering. Submitted for publication.

# Author:
Dewei Wang (deweiwang@stat.sc.edu)
