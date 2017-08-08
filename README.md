# GLMaSPU  
**GLMaSPU** is a new R package that makes it incredibly easy to implement some testing methods for high dimensional parameters in generalized linear models with R. See [www.wuchong.org/GLMaSPU.html] (www.wuchong.org/GLMaSPU.html) for details.

## Benefits
* It provides several tests for testing high dimensional parameters in generalized linear models
* Unlike traditional low dimensional situation with p < n, GLMaSPU package can deal with high dimensional situation with p > n
* A new and powerful test, called aSPU, has been implemented in GLMaSPU

## Features

* Implement resample based version of Goeman (2011) and Guo (2016)
* Implement both asymptotic version and resample version of aSPU test, a powerful test for a wide high dimensional scenarios.

## Installation
To install the stable version from CRAN, simply run the following from an R console: **We already upload the package to CRAN and we expect it will be avaiable in the near future.**

```r
install.packages("GLMaSPU")
```

To install the latest development builds directly from GitHub, run this instead:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("ChongWu-Biostat/GLMaSPU")
```

## Using aSPU
Several tests for high dimensional generalized linear models have been proposed recently, however, they are mainly based on a sum of squares of the score vector and only powerful under certain limited alternative hypotheses. In practice, since the signals in true alternative hypothesis may be sparse or dense or between, the existing tests may not be powerful. In this package, we implement an adaptive test that maintains high power across a wide range of scenarios. To calculate its p-value, its asymptotic null distribution is derived. For more details, see the following two papers.

* Pan, W., Kim, J., Zhang, Y., Shen, X., & Wei, P. (2014). A powerful and adaptive association test for rare variants. *Genetics*, 197(4), 1081-1095.
* Wu, C., Xu, G. & Pan, W. An Adaptive test on high dimensional parameters in generalized linear models, submitted. 

```r
library("GLMaSPU")
## generate the data
p = 200n = 100beta = c(1,3,3)s = 0.15signal.r = 0.08
non.zero = floor(p * s)seed = 2alpha = c(rep(signal.r,non.zero),rep(0,p-non.zero))dat = generate_data(seed, n = n, p = p, beta = beta,alpha = alpha)cov = dat$ZX = dat$XY = dat$YaSPU_perm(Y, X, cov = cov, pow = c(1:6, Inf),resample = "perm", model = "gaussian", n.perm = 1000) 
aSPU_apval(Y, X, cov = cov, pow = c(1:6, Inf),resample = "perm", model = "gaussian", n.perm = 1000)
# The p-values are similar to the asymptotic based one
```
