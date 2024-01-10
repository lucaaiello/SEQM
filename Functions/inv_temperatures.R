inv_temperatures <- function(d_log_temp, L){
  temp <- rep(1,L)
  for(k in 1:(L-1)){
    temp[k+1] <- temp[k] + exp(d_log_temp[k])
  }
  return(1/temp)
}
