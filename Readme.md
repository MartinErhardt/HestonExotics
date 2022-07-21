
![GCC](https://img.shields.io/static/v1?logo=github&label=GCC&message=build passing&color=Blue)
![License](https://img.shields.io/static/v1?label=License&message=MPL-2.0&color=blue)
![docs](https://img.shields.io/static/v1?label=docs&message=doxygen&color=green)

hexo
===============================================
Simple derivative pricing tool for educational purposes based on the Heston model
### Installation
In order to generate a distributable tarball from head, type
```
autoreconf --install &&
  ./configure &&
  make dist
```
The autoconf script will detect dependencies and if necessary create a makefile, that builds blas, curl, eigen3, fftw3, lapack, libcurl and sqlite3.
Documention is doxygen is available:
```
$ make doc
/** installs documentation into doc/html/index.html*/
```
### Usage examples
Here are some simple usecases:
```
$ ./hexo -c AMZN
/** calibrates the Heston model to Amazon stock as underlying and stores the resulting parameters in ParamsDB.db*/
$ ./hexo -p asian all GOOG
/** prices all available European Call options on Google stock as if they were arithmetic Asian call Options and prints the results to stdout*/
$ ./hexo -t rng|distr|pricing|gradient|levmar|rng
/** run tests*/
```
### Key observations/evaluation
One can observe a high degree of fluctuation in the rate-reversion and volatility of volatility parameters κ and ρ respectively resulting from the calibration process. One day κ might be one, later that day it could be 20. At first I thought this was an implementation error, but I was able to successfully test my implementation of the calibration algorithm [Ortiz-Garcia(2021)](https://www.mdpi.com/2227-7390/9/5/529/pdf) against that one by [Eudald Romo Grau](https://github.com/eudaldrg/SWIFTOptionCalibration). Moreover as it turns out, this issue is inherent to the optimization problem, in which we try to minimize the squared difference between observed market prices and those resulting if a Heston model to the parameters over which we minimize is assumed. If f is the objective function mapping parameters to said prices in the heston model, then the derivative of that function, which is used to determine the first iterations of a gradient decent is numerically instable in the parameter κ. This is consistent with what the academic literature describes:

> The Hessian matrix is ill-conditioned with a condition number of 3.978×10e+6.  
> The elements ∂2f(θ)/∂κ^2 and ∂2f(θ)/∂ρ^2 are of a much smaller order than the others.  
> This suggests that the objective function, when around the optimum, is less sensitive to changes along κ and ρ.  
> In other words, the objective function is more stretched along these two axes (..)  
> The ratio between ∂2f(θ)/∂κ^2 and ∂2f(θ)/∂v^2 is of order 10−6, which indicates a great disparity in sensitivity: 
> changing 1 unit of v_m is comparable to changing 106 units of κ.  
>-- <cite>[Yiran Cui(2017)](https://arxiv.org/pdf/1511.08718.pdf)</cite> 

However this is not as much of an issue as one might assume, as we then proceed to use the obtained parameters to compute prices of exotic options.
To illustrate this point, say we wanted to determine the prices of European call options (even though this is besides the point).  
Then the same numerical stability of f in κ, that makes it so difficult to determine a minimal κ, works to our advantage in this step.  
Indeed the final results for prices of arithmetic call options are much more consistent, though I have not been able to verify them in any way, because of a lack of software to test against.  
Other issues include:
- 3x performance, when using Intel compiler
- standard out is locked, when compiled in MinGW64
