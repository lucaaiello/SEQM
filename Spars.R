Spars <- function(lat,lon,lat0,lon0,t0,depth){
  N <- length(lat)
  
  par <- matrix(nrow = N, ncol = 2)
  par[,1] <- hyp_dist(lat,lon,lat0,lon0,depth)/4.5 + t0
  par[,2] <- par[,1] + 3.5
  
  return(par)
}