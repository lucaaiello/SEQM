log_target <- function(y,delta,p_orig,p,alpha_orig,alpha,
                       lats,lons,lat0_orig,lat0,lon0_orig,lon0,
                       t0_orig,t0,d0_orig,d0,vP_orig,vP,
                       tau_orig,tau,invlambda_orig,invlambda,
                       latstar,lonstar,sigmalat,sigmalon,
                       r,s,a,b,c,d){
  
  L <- length(p_orig)
  N <- length(y)
  
  out <- vector(length = L)
  
  log_aux <- matrix(0, nrow = N, ncol = L)
  
  for (l in 1:L) {
    
    log_aux[which(delta!=0),l] <- log_h(y[which(delta!=0)],
                                        p_orig[l],alpha_orig[l],
                                        lats[which(delta!=0)],lons[which(delta!=0)],
                                        lat0_orig[l],lon0_orig[l],
                                        t0_orig[l],d0_orig[l],
                                        vP_orig[l],tau_orig[l],
                                        invlambda_orig[l])
    
    out[l] <- sum(log_aux[,l] + log_s(y,
                                      p_orig[l],alpha_orig[l],
                                      lats,lons,
                                      lat0_orig[l],lon0_orig[l],
                                      t0_orig[l],d0_orig[l],
                                      vP_orig[l],tau_orig[l],
                                      invlambda_orig[l])) - 
      (lat0_orig[l] - latstar)^2/(2*sigmalat^2) - 
      (lon0_orig[l] - lonstar)^2/(2*sigmalon^2) + 
      (s-1)*d0[l] - r*exp(d0[l]) +
      (a-1)*(alpha[l] - log1pexp(alpha[l])) + (b-1)*(-log1pexp(alpha[l])) +
      (c-1)*vP[l] - d*exp(vP[l]) +
      lat0[l] - 2*log1pexp(lat0[l]) + 
      lon0[l] - 2*log1pexp(lon0[l]) + 
      d0[l] +
      (t0[l] - 2*log1pexp(t0[l])) +
      (alpha[l] - 2*log1pexp(alpha[l])) + 
      (p[l] - 2*log1pexp(p[l])) +
      vP[l] +
      (tau[l] - 2*log1pexp(tau[l])) +
      (invlambda[l] - 2*log1pexp(invlambda[l]))
  }
  
  return(out)
  
}
