log_h <- function(y,p,alpha,lat,lon,lat0,lon0,t0,d0,vP,tau,invlambda){

  mu <- mus(lat,lon,lat0,lon0,t0,d0,vP)

  mu_P <- mu[,1]
  mu_S <- mu[,2]

  logh <- log((1/invlambda)*(p + (1-p)*(alpha*pnorm(y,mu_P,tau,lower.tail = FALSE) +
                                          (1-alpha)*pnorm(y,mu_S,tau,lower.tail = FALSE))) +
                ((1-p)*(alpha*dnorm(y,mu_P,tau) +
                          (1-alpha)*dnorm(y,mu_S,tau)))) -
    log(p + (1-p)*(alpha*pnorm(y,mu_P,tau,lower.tail = FALSE) +
                     (1-alpha)*pnorm(y,mu_S,tau,lower.tail = FALSE)))

  return(logh)

}