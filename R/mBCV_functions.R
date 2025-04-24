mBCV <- function(Y_in){
  # check
  if(!all(is.vector(Y_in), is.numeric(Y_in))) stop("Input data mustb be numeric vector")

  # Initialize package
  KDEmain <- new(cKDE, Y_in, 100, 100) 	
  lower = KDEmain$LU_Bound_h[1] ; upper = KDEmain$LU_Bound_h[2] ; tol = 0.05*lower
  
  # mBCV
  f_mbcv <- function(h0){ KDEmain$mBCV_obj_fn(h0) }
  h_hat = optimize( f_mbcv, c(lower, upper), tol = tol)$minimum
  
  result <- list(bw=h_hat, IACT = KDEmain$IACT(), Y_in = Y_in)
  
  class(result) <- "mBCV_obj"
  
  return(result)
}

print.mBCV_obj <- function(x,...){
  cat("    Modified Biased Cross-Validation\n")
  cat("Data Size: ", length(x$Y_in)," obs. \n")

  print(summary(x$Y_in)  )
  
  cat("------------------------------\n")
  cat("Bandwidth (bw):", round(x$bw, 4), "\n")
  cat("Intergrated autocorrelation time (IACT): ", round(x$IACT, 4), "\n")
}

plot.mBCV_obj <- function(x, main=NULL, xlab="", ...){
  den = density(x$Y_in, bw=x$bw)
  
  hist(x$Y_in, freq=FALSE, 
       main=main, xlab=xlab)
  
  lines(den$x, den$y, col='blue', lwd=2)
  
  legend("topright", legend = c(paste0(c("bw = ","IACT = "),  round(c(x$bw,x$IACT), 3))), bty="n", text.width = 1.3)
}