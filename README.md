# GFabs
  A Group Forward and Backward Stagewise (GFabs) algorithm for Group penalized problems.
 
  GFabs uses coordinate descent with a fixed step size which consists of both forward and backward steps. At each step, the first-order Taylor's expansion is used to reflect the main part of the increment. 

# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("XiaoZhangryy/GFabs")

# Usage

   - [x] [GFabs-manual](https://github.com/XiaoZhangryy/GFabs/blob/master/inst/GFabs-manual.pdf) ------------ Details of the usage of the package.
   
# Example

    library(GFabs)
    library(mvtnorm)

    sigma = outer(1:20, 1:20, FUN = function(x, y) 0.3^(abs(x - y)))
    x     = rmvnorm(100, mean = rep(0,20), sigma = sigma)
    u     = runif(100)
    b     = cbind(5*sin(2*pi*u), 5*cos(2*pi*u), 5, -5, matrix(0, 100, 16))
    error = c(0.7*rnorm(100)+0.3*rcauchy(100))
    y     = rowSums(x * b) + error
    fit   <- GFabs_vc(x, y, u)
    
# Replicate simulation results in Zhang et al.(2021)

All the simulation results can be reproduced by using the codes at [simulation](https://github.com/XiaoZhangryy/GFabs/blob/master/simulations). 

# References

Model Selection for Transformation Model with High Dimensional Varying Coefficients. Manuscript.

# Development
The R-package is developed by Xiao Zhang (zhangxiao_0422@163.com) and Xingjie Shi.
