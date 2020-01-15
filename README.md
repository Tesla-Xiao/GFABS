# GFABS
  A Group Forward and Backward Stagewise (GFabs) algorithm for smoothed partial rank loss function (SPR) with the Group Lasso penalty.
 
  GFabs uses coordinate descent with a fixed step size which consists of both forward and backward steps. At each step, the first-order Taylor's expansion is used to reflect the main part of the increment. 

# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("Tesla-Xiao/GFABS")

# Usage

   - [x] [GFabs-manual](https://github.com/Tesla-Xiao/GFabs/blob/master/inst/GFabs-manual.pdf) ------------ Details of the usage of the package.
   
# Example

    library(GFabs)

    n = 100
    p = 100
    sigma = outer(1:p, 1:p, FUN = function(x, y) 0.3^(abs(x - y)))
    x = rmvnorm(n, mean = rep(0,p), sigma = sigma)
    u = runif(100)
    b = cbind(5*sin(2*pi*u), 5*cos(2*pi*u), 5*exp(2*u-1)-5, 5, -5, matrix(0, nrow=n, ncol=p-5))
    eta = x %*% b
    e = c(0.7*rnorm(n)+0.3*rcauchy(n))
    g = function(x) x
    y = g(eta + e)
    fit <- hvcspr(y, x, u)
    
# References

Model Selection for Transformation Model with High Dimensional Varying Coefficients. Manuscript.

# Development
The R-package is developed by Xiao Zhang (zhangxiao_0422@163.com) and Xingjie.
