\name{mBCV}
\alias{mBCV}
\alias{print.mBCV_obj}
\title{
Calculate Optimal Bandwidth in Kernel Density Estimation
}
\description{
Calculate the optimal bandwidth estimate for the kernel density estimator with a Markov chain Monte Carlo sample using modified biased cross-validation method.
}
\usage{
mBCV(Y_in)

\method{print}{mBCV_obj}(x, \dots)
}
\arguments{
  \item{Y_in}{ 
    data from which the estimate is to be computed.
  }
  \item{x}{
    object of class \code{mBCV_obj}; result of a call to \code{mBCV()}.
  }
  \item{...}{
    further arguments passed to or from other methods.
  }
}
\value{
  \code{mBCV} returns a list of the following conmponents:
  
  \item{bw}{optimal bandwidth.}
  \item{IACT}{intergrated autocorrelation time.}
  \item{Y_in}{input data.}
}

\examples{
res = mBCV(simMCMC)
den = density(res$Y_in, bw=res$bw)
hist(res$Y_in, xlim=range(den$x), freq=FALSE, 
      main="Histogram and Density Estimates", xlab="")
lines(den$x, den$y, col='blue', lwd=2)
}
