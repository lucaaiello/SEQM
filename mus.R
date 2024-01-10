mus <- function(lat,lon,lat0,lon0,t0,d0,vP){
  N <- length(lat)
  vS <- vP/sqrt(3)
  
  par <- matrix(nrow = N, ncol = 2)
  par[,1] <- hyp_dist(lat,lon,lat0,lon0,d0)/vP + t0 + 1.75
  par[,2] <- hyp_dist(lat,lon,lat0,lon0,d0)/vS + t0 + 1.75
  
  return(par)
}