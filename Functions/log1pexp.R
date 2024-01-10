log1pexp <- function(x){
  
  L <- length(x)
  
  out <- vector(length = L)
  
  for (l in 1:L) {
    if(x[l] < -37){
      out[l] <- exp(x[l])
    }else if(x[l] > -37 & x[l] <= 18){
      out[l] <- log1p(exp(x[l]))
    }else if(x[l] > 18 & x[l] <= 33.3){
      out[l] <- x[l] + exp(-x[l])
    }else if(x[l] > 33.3){
      out[l] <- x[l]
    }  
  }
  return(out)  
}
