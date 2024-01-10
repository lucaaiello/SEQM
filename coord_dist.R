coord_dist <- function(lat,lon,lat0,lon0){
  lat <- deg2rad(lat)
  lon <- deg2rad(lon)
  
  lat0 <- deg2rad(lat0)
  lon0 <- deg2rad(lon0)
  
  dist <- 2 * 6371 * asin(sqrt((sin((lat-lat0)/2))^2 + 
                                 cos(lat)*cos(lat0)*(sin((lon-lon0)/2))^2))
  return(dist)
  
}
