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

   - [x] [GFabs-manual](https://github.com/Tesla-Xiao/GFabs/blob/master/inst/GFABS-manual.pdf) ------------ Details of the usage of the package.
   
# Example

    library(GFABS)

    W = matrix(rnorm(80), 20, 4)
    b = c(5, 5, -5, -5)
    e = c(0.7*rnorm(20)+0.3*rcauchy(20))
    y = W %*% b + e
    group <- c(1, 1, 2, 2)
    fit <- GFabs(W, y, group)
    
# References

Model Selection for Transformation Model with High Dimensional Varying Coefficients. Manuscript.

# Development
The R-package is developed by Xiao Zhang (zhangxiao_0422@163.com) and Xingjie Shi.
