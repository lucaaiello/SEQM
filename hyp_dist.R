hyp_dist <- function(lat,lon,lat0,lon0,depth){
  
  epi_dist <- coord_dist(lat,lon,lat0,lon0)
  
  dist <- sqrt(depth^2 + 4*6371*(6371 - depth)*sin(epi_dist/(2*6371))^2)
  
  return(dist)
  
}
