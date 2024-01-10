mh_surv_quake_2 <- function(G,burnin,thin,
                            eta_1,eta_2,eta_3,eta_4,eta_5,eta_6,
                            s_1,s_2,s_3,s_4,s_5,s_6,
                            th0,y,lats,lons,delta,
                            latstar,lonstar,sigmalat,sigmalon,
                            r,s,a,b,c,d,
                            beta,d_log_temp) 
{
  library(MASS)
  
  N <- length(y)
  p <- 9
  L <- length(beta)
  
  iterations <- burnin+thin*G
  g <- 1
  
  theta <- array(dim = c(G,p,L))
  
  lat0_curr_orig <- rep(th0[1],L)
  lon0_curr_orig <- rep(th0[2],L)
  d0_curr_orig <- rep(th0[3],L)
  t0_curr_orig <- rep(th0[4],L)
  alpha_curr_orig <- rep(th0[5],L)
  p_curr_orig <- rep(th0[6],L)
  vP_curr_orig <- rep(th0[7],L)
  tau_curr_orig <- rep(th0[8],L)
  invlambda_curr_orig <- rep(th0[9],L)
  
  lat0_curr <- trans_par(lat0_curr_orig,-90,90)
  lon0_curr <- trans_par(lon0_curr_orig,-180,180)
  d0_curr <- log(d0_curr_orig)
  t0_curr <- trans_par(t0_curr_orig,0,120)
  alpha_curr <- trans_par(alpha_curr_orig,0,1)
  p_curr <- trans_par(p_curr_orig,0,1)
  vP_curr <- log(vP_curr_orig)
  tau_curr <- trans_par(tau_curr_orig,0,3.5)
  invlambda_curr <- trans_par(invlambda_curr_orig,800,80000)
  
  acc_1 <- rep(0,L)
  acc_2 <- rep(0,L)
  acc_3 <- rep(0,L)
  acc_4 <- rep(0,L)
  acc_5 <- rep(0,L)
  acc_6 <- rep(0,L)
  
  acc_sw <- rep(0,(L-1))
  N_sw <- rep(0,(L-1))
  
  m_1 <- rbind(lat0_curr,lon0_curr,d0_curr,t0_curr)
  m_2 <- alpha_curr
  m_3 <- p_curr
  m_4 <- vP_curr
  m_5 <- tau_curr
  m_6 <- invlambda_curr
  
  R_1 <- array(dim = c(4,4,L))
  R_2 <- vector(length = L)
  R_3 <- vector(length = L)
  R_4 <- vector(length = L)
  R_5 <- vector(length = L)
  R_6 <- vector(length = L)
  
  for (l in 1:L) {
    
    R_1[,,l] <- diag(c(0.1,0.1,10,1))
    R_2[l] <- 0.1
    R_3[l] <- 0.1
    R_4[l] <- 0.1
    R_5[l] <- 0.1
    R_6[l] <- 0.1
    
  }
  
  log_post_curr_6 <- log_target(y,delta,
                                p_curr_orig,p_curr,
                                alpha_curr_orig,alpha_curr,
                                lats,lons,
                                lat0_curr_orig,lat0_curr,
                                lon0_curr_orig,lon0_curr,
                                t0_curr_orig,t0_curr,
                                d0_curr_orig,d0_curr,
                                vP_curr_orig,vP_curr,
                                tau_curr_orig,tau_curr,
                                invlambda_curr_orig,invlambda_curr,
                                latstar,lonstar,sigmalat,sigmalon,
                                r,s,a,b,c,d)
  
  for (iter in 1:iterations) {
    # Proposal for lat0, lon0, d0, t0
    
    prop_1 <- matrix(nrow = 4,ncol = L)
    
    for (l in 1:L) {
      prop_1[,l] <- mvrnorm(n = 1, mu = c(lat0_curr[l],lon0_curr[l],d0_curr[l],t0_curr[l]), Sigma = eta_1[,,l])  
    }
    
    lat0_prop <- prop_1[1,]
    lon0_prop <- prop_1[2,]
    d0_prop <- prop_1[3,]
    t0_prop <- prop_1[4,]
    
    lat0_prop_orig <- inv_trans_par(lat0_prop,-90,90)
    lon0_prop_orig <- inv_trans_par(lon0_prop,-180,180)
    d0_prop_orig <- exp(d0_prop)
    t0_prop_orig <- inv_trans_par(t0_prop,0,120)
    
    log_post_curr_1 <- log_post_curr_6
    
    log_post_prop_1 <- log_target(y,delta,
                                  p_curr_orig,p_curr,
                                  alpha_curr_orig,alpha_curr,
                                  lats,lons,
                                  lat0_prop_orig,lat0_prop,
                                  lon0_prop_orig,lon0_prop,
                                  t0_prop_orig,t0_prop,
                                  d0_prop_orig,d0_prop,
                                  vP_curr_orig,vP_curr,
                                  tau_curr_orig,tau_curr,
                                  invlambda_curr_orig,invlambda_curr,
                                  latstar,lonstar,sigmalat,sigmalon,
                                  r,s,a,b,c,d)
    
    log_ratio_1 <- log_post_prop_1 - log_post_curr_1
    
    log_alpha_1 <- pmin(0, beta*log_ratio_1)
    
    u_1 <- runif(L)
    
    for (l in 1:L) {
      if (log(u_1[l]) < log_alpha_1[l]){
        #Accept the move
        acc_1[l] <- acc_1[l]+1
        
        lat0_curr[l] <- lat0_prop[l]
        lon0_curr[l] <- lon0_prop[l]
        d0_curr[l] <- d0_prop[l]
        t0_curr[l] <- t0_prop[l] 
        
        lat0_curr_orig[l] <- lat0_prop_orig[l]
        lon0_curr_orig[l] <- lon0_prop_orig[l]
        d0_curr_orig[l] <- d0_prop_orig[l]
        t0_curr_orig[l] <- t0_prop_orig[l]
        
        log_post_curr_1[l] <- log_post_prop_1[l]
      }
      
      s_1[l] <- s_1[l] + (iter+1)^(-0.6) * (exp(log_alpha_1[l]) - 0.234)
      
      dy_1 <- c(lat0_curr[l],lon0_curr[l],d0_curr[l],t0_curr[l]) - m_1[,l]
      
      # Covariance
      R_1[,,l] <- (1-(iter+1)^(-0.6))*R_1[,,l] + (iter+1)^(-0.6)*(dy_1%*%t(dy_1))
      
      # Mean update
      m_1[,l] <- m_1[,l] + (iter+1)^(-0.6)*dy_1
      
      eta_1[,,l] <- exp(s_1[l]) * R_1[,,l]
    }
    
    # Proposal for alpha
    
    alpha_prop <- mvrnorm(n = 1, mu = alpha_curr, Sigma = diag(eta_2))
    
    alpha_prop_orig <- inv_trans_par(alpha_prop,0,1)
    
    log_post_curr_2 <- log_post_curr_1
    
    log_post_prop_2 <- log_target(y,delta,
                                  p_curr_orig,p_curr,
                                  alpha_prop_orig,alpha_prop,
                                  lats,lons,
                                  lat0_curr_orig,lat0_curr,
                                  lon0_curr_orig,lon0_curr,
                                  t0_curr_orig,t0_curr,
                                  d0_curr_orig,d0_curr,
                                  vP_curr_orig,vP_curr,
                                  tau_curr_orig,tau_curr,
                                  invlambda_curr_orig,invlambda_curr,
                                  latstar,lonstar,sigmalat,sigmalon,
                                  r,s,a,b,c,d)
    
    log_ratio_2 <- log_post_prop_2 - log_post_curr_2
    
    log_alpha_2 <- pmin(0, beta*log_ratio_2)
    
    u_2 <- runif(L)
    
    for (l in 1:L) {
      if (log(u_2[l]) < log_alpha_2[l]){
        #Accept the move
        acc_2[l] <- acc_2[l]+1
        
        alpha_curr[l] <- alpha_prop[l]
        
        alpha_curr_orig[l] <- alpha_prop_orig[l]
        
        log_post_curr_2[l] <- log_post_prop_2[l]
      }
      
      s_2[l] <- s_2[l] + (iter+1)^(-0.6) * (exp(log_alpha_2[l]) - 0.44)
      
      dy_2 <- alpha_curr[l] - m_2[l]
      
      # Covariance
      R_2[l] <- (1-(iter+1)^(-0.6))*R_2[l] + (iter+1)^(-0.6)*dy_2^2
      
      # Mean update
      m_2[l] <- m_2[l] + (iter+1)^(-0.6)*dy_2
      
      eta_2[l] <- exp(s_2[l]) * R_2[l]  
    }
    
    # Proposal for p
    
    p_prop <- mvrnorm(n = 1, mu = p_curr, Sigma = diag(eta_3))
    
    p_prop_orig <- inv_trans_par(p_prop,0,1)
    
    log_post_curr_3 <- log_post_curr_2
    
    log_post_prop_3 <- log_target(y,delta,
                                  p_prop_orig,p_prop,
                                  alpha_curr_orig,alpha_curr,
                                  lats,lons,
                                  lat0_curr_orig,lat0_curr,
                                  lon0_curr_orig,lon0_curr,
                                  t0_curr_orig,t0_curr,
                                  d0_curr_orig,d0_curr,
                                  vP_curr_orig,vP_curr,
                                  tau_curr_orig,tau_curr,
                                  invlambda_curr_orig,invlambda_curr,
                                  latstar,lonstar,sigmalat,sigmalon,
                                  r,s,a,b,c,d)
    
    log_ratio_3 <- log_post_prop_3 - log_post_curr_3
    
    log_alpha_3 <- pmin(0, beta*log_ratio_3)
    
    u_3 <- runif(L)
    
    for (l in 1:L) {
      if (log(u_3[l]) < log_alpha_3[l]){
        #Accept the move
        acc_3[l] <- acc_3[l]+1
        
        p_curr[l] <- p_prop[l]
        
        p_curr_orig[l] <- p_prop_orig[l]
        
        log_post_curr_3[l] <- log_post_prop_3[l]
      }
      
      s_3[l] <- s_3[l] + (iter+1)^(-0.6) * (exp(log_alpha_3[l]) - 0.44)
      
      dy_3 <- p_curr[l] - m_3[l]
      
      # Covariance
      R_3[l] <- (1-(iter+1)^(-0.6))*R_3[l] + (iter+1)^(-0.6)*dy_3^2
      
      # Mean update
      m_3[l] <- m_3[l] + (iter+1)^(-0.6)*dy_3
      
      eta_3[l] <- exp(s_3[l]) * R_3[l]  
    }
    
    # Proposal for vP
    
    vP_prop <- mvrnorm(n = 1, mu = vP_curr, Sigma = diag(eta_4))
    
    vP_prop_orig <- exp(vP_prop)
    
    log_post_curr_4 <- log_post_curr_3
    
    log_post_prop_4 <- log_target(y,delta,
                                  p_curr_orig,p_curr,
                                  alpha_curr_orig,alpha_curr,
                                  lats,lons,
                                  lat0_curr_orig,lat0_curr,
                                  lon0_curr_orig,lon0_curr,
                                  t0_curr_orig,t0_curr,
                                  d0_curr_orig,d0_curr,
                                  vP_prop_orig,vP_prop,
                                  tau_curr_orig,tau_curr,
                                  invlambda_curr_orig,invlambda_curr,
                                  latstar,lonstar,sigmalat,sigmalon,
                                  r,s,a,b,c,d)
    
    log_ratio_4 <- log_post_prop_4 - log_post_curr_4
    
    log_alpha_4 <- pmin(0, beta*log_ratio_4)
    
    u_4 <- runif(L)
    
    for (l in 1:L) {
      if (log(u_4[l]) < log_alpha_4[l]){
        #Accept the move
        acc_4[l] <- acc_4[l]+1
        
        vP_curr[l] <- vP_prop[l]
        
        vP_curr_orig[l] <- vP_prop_orig[l]
        
        log_post_curr_4[l] <- log_post_prop_4[l]
      }
      
      s_4[l] <- s_4[l] + (iter+1)^(-0.6) * (exp(log_alpha_4[l]) - 0.44)
      
      dy_4 <- vP_curr[l] - m_4[l]
      
      # Covariance
      R_4[l] <- (1-(iter+1)^(-0.6))*R_4[l] + (iter+1)^(-0.6)*dy_4^2
      
      # Mean update
      m_4[l] <- m_4[l] + (iter+1)^(-0.6)*dy_4
      
      eta_4[l] <- exp(s_4[l]) * R_4[l]  
    }
    
    # Proposal for tau
    
    tau_prop <- mvrnorm(n = 1, mu = tau_curr, Sigma = diag(eta_5))
    
    tau_prop_orig <- inv_trans_par(tau_prop,0,3.5)
    
    log_post_curr_5 <- log_post_curr_4
    
    log_post_prop_5 <- log_target(y,delta,
                                  p_curr_orig,p_curr,
                                  alpha_curr_orig,alpha_curr,
                                  lats,lons,
                                  lat0_curr_orig,lat0_curr,
                                  lon0_curr_orig,lon0_curr,
                                  t0_curr_orig,t0_curr,
                                  d0_curr_orig,d0_curr,
                                  vP_curr_orig,vP_curr,
                                  tau_prop_orig,tau_prop,
                                  invlambda_curr_orig,invlambda_curr,
                                  latstar,lonstar,sigmalat,sigmalon,
                                  r,s,a,b,c,d)
    
    log_ratio_5 <- log_post_prop_5 - log_post_curr_5
    
    log_alpha_5 <- pmin(0, beta*log_ratio_5)
    
    u_5 <- runif(L)
    
    for (l in 1:L) {
      if (log(u_5[l]) < log_alpha_5[l]){
        #Accept the move
        acc_5[l] <- acc_5[l]+1
        
        tau_curr[l] <- tau_prop[l]
        
        tau_curr_orig[l] <- tau_prop_orig[l]
        
        log_post_curr_5[l] <- log_post_prop_5[l]
      }
      
      s_5[l] <- s_5[l] + (iter+1)^(-0.6) * (exp(log_alpha_5[l]) - 0.44)
      
      dy_5 <- tau_curr[l] - m_5[l]
      
      # Covariance
      R_5[l] <- (1-(iter+1)^(-0.6))*R_5[l] + (iter+1)^(-0.6)*dy_5^2
      
      # Mean update
      m_5[l] <- m_5[l] + (iter+1)^(-0.6)*dy_5
      
      eta_5[l] <- exp(s_5[l]) * R_5[l]  
    }
    
    # Proposal for invlambda
    
    invlambda_prop <- mvrnorm(n = 1, mu = invlambda_curr, Sigma = diag(eta_6))
    
    invlambda_prop_orig <- inv_trans_par(invlambda_prop,800,80000)
    
    log_post_curr_6 <- log_post_curr_5
    
    log_post_prop_6 <- log_target(y,delta,
                                  p_curr_orig,p_curr,
                                  alpha_curr_orig,alpha_curr,
                                  lats,lons,
                                  lat0_curr_orig,lat0_curr,
                                  lon0_curr_orig,lon0_curr,
                                  t0_curr_orig,t0_curr,
                                  d0_curr_orig,d0_curr,
                                  vP_curr_orig,vP_curr,
                                  tau_curr_orig,tau_curr,
                                  invlambda_prop_orig,invlambda_prop,
                                  latstar,lonstar,sigmalat,sigmalon,
                                  r,s,a,b,c,d)
    
    log_ratio_6 <- log_post_prop_6 - log_post_curr_6
    
    log_alpha_6 <- pmin(0, beta*log_ratio_6)
    
    u_6 <- runif(L)
    
    for (l in 1:L) {
      if (log(u_6[l]) < log_alpha_6[l]){
        #Accept the move
        acc_6[l] <- acc_6[l]+1
        
        invlambda_curr[l] <- invlambda_prop[l]
        
        invlambda_curr_orig[l] <- invlambda_prop_orig[l]
        
        log_post_curr_6[l] <- log_post_prop_6[l]
      }
      
      s_6[l] <- s_6[l] + (iter+1)^(-0.6) * (exp(log_alpha_6[l]) - 0.44)
      
      dy_6 <- invlambda_curr[l] - m_6[l]
      
      # Covariance
      R_6[l] <- (1-(iter+1)^(-0.6))*R_6[l] + (iter+1)^(-0.6)*dy_6^2
      
      # Mean update
      m_6[l] <- m_6[l] + (iter+1)^(-0.6)*dy_6
      
      eta_6[l] <- exp(s_6[l]) * R_6[l]  
    }
    
    # coupling the chains
    
    swap_idx <- sample(1:(L-1),1)
    
    swap_idx <- c(swap_idx,(swap_idx+1))
    
    N_sw[swap_idx[1]] <- N_sw[swap_idx[1]] + 1
    
    log_A <- (beta[swap_idx[1]]-beta[swap_idx[2]])*
      (log_post_curr_6[swap_idx[2]] - log_post_curr_6[swap_idx[1]])
    
    u_A <- runif(1)
    
    log_A <- min(0,log_A)
    
    if(log(u_A)<log_A){
      
      acc_sw[swap_idx[1]] <- acc_sw[swap_idx[1]] + 1
      
      aux_lat0_curr <- lat0_curr[swap_idx[1]]
      aux_lon0_curr <- lon0_curr[swap_idx[1]]
      aux_d0_curr <- d0_curr[swap_idx[1]]
      aux_t0_curr <- t0_curr[swap_idx[1]]
      aux_alpha_curr <- alpha_curr[swap_idx[1]]
      aux_p_curr <- p_curr[swap_idx[1]]
      aux_vP_curr <- vP_curr[swap_idx[1]]
      aux_tau_curr <- tau_curr[swap_idx[1]]
      aux_invlambda_curr <- invlambda_curr[swap_idx[1]]
      aux_log_post_curr <- log_post_curr_6[swap_idx[1]]
      
      aux_lat0_curr_orig <- lat0_curr_orig[swap_idx[1]]
      aux_lon0_curr_orig <- lon0_curr_orig[swap_idx[1]]
      aux_d0_curr_orig <- d0_curr_orig[swap_idx[1]]
      aux_t0_curr_orig <- t0_curr_orig[swap_idx[1]]
      aux_alpha_curr_orig <- alpha_curr_orig[swap_idx[1]]
      aux_p_curr_orig <- p_curr_orig[swap_idx[1]]
      aux_vP_curr_orig <- vP_curr_orig[swap_idx[1]]
      aux_tau_curr_orig <- tau_curr_orig[swap_idx[1]]
      aux_invlambda_curr_orig <- invlambda_curr_orig[swap_idx[1]]
      
      lat0_curr[swap_idx[1]] <- lat0_curr[swap_idx[2]]
      lon0_curr[swap_idx[1]] <- lon0_curr[swap_idx[2]]
      d0_curr[swap_idx[1]] <- d0_curr[swap_idx[2]]
      t0_curr[swap_idx[1]] <- t0_curr[swap_idx[2]]
      alpha_curr[swap_idx[1]] <- alpha_curr[swap_idx[2]]
      p_curr[swap_idx[1]] <- p_curr[swap_idx[2]]
      vP_curr[swap_idx[1]] <- vP_curr[swap_idx[2]]
      tau_curr[swap_idx[1]] <- tau_curr[swap_idx[2]]
      invlambda_curr[swap_idx[1]] <- invlambda_curr[swap_idx[2]]
      log_post_curr_6[swap_idx[1]] <- log_post_curr_6[swap_idx[2]]
      
      lat0_curr_orig[swap_idx[1]] <- lat0_curr_orig[swap_idx[2]]
      lon0_curr_orig[swap_idx[1]] <- lon0_curr_orig[swap_idx[2]]
      d0_curr_orig[swap_idx[1]] <- d0_curr_orig[swap_idx[2]]
      t0_curr_orig[swap_idx[1]] <- t0_curr_orig[swap_idx[2]]
      alpha_curr_orig[swap_idx[1]] <- alpha_curr_orig[swap_idx[2]]
      p_curr_orig[swap_idx[1]] <- p_curr_orig[swap_idx[2]]
      vP_curr_orig[swap_idx[1]] <- vP_curr_orig[swap_idx[2]]
      tau_curr_orig[swap_idx[1]] <- tau_curr_orig[swap_idx[2]]
      invlambda_curr_orig[swap_idx[1]] <- invlambda_curr_orig[swap_idx[2]]
      
      lat0_curr[swap_idx[2]] <- aux_lat0_curr
      lon0_curr[swap_idx[2]] <- aux_lon0_curr
      d0_curr[swap_idx[2]] <- aux_d0_curr
      t0_curr[swap_idx[2]] <- aux_t0_curr
      alpha_curr[swap_idx[2]] <- aux_alpha_curr
      p_curr[swap_idx[2]] <- aux_p_curr
      vP_curr[swap_idx[2]] <- aux_vP_curr
      tau_curr[swap_idx[2]] <- aux_tau_curr
      invlambda_curr[swap_idx[2]] <- aux_invlambda_curr
      log_post_curr_6[swap_idx[2]] <- aux_log_post_curr
      
      lat0_curr_orig[swap_idx[2]] <- aux_lat0_curr_orig
      lon0_curr_orig[swap_idx[2]] <- aux_lon0_curr_orig
      d0_curr_orig[swap_idx[2]] <- aux_d0_curr_orig
      t0_curr_orig[swap_idx[2]] <- aux_t0_curr_orig
      alpha_curr_orig[swap_idx[2]] <- aux_alpha_curr_orig
      p_curr_orig[swap_idx[2]] <- aux_p_curr_orig
      vP_curr_orig[swap_idx[2]] <- aux_vP_curr_orig
      tau_curr_orig[swap_idx[2]] <- aux_tau_curr_orig
      invlambda_curr_orig[swap_idx[2]] <- aux_invlambda_curr_orig
      
    }
    
    d_log_temp[swap_idx[1]] <- d_log_temp[swap_idx[1]]+(iter+1)^(-0.6)*(exp(log_A) - 0.44)
    
    beta <- inv_temperatures(d_log_temp, L)
    
    if( (iter>burnin) && (iter%%thin==0) ){ 
      ## if iter is larger than burn
      ## and iter is multiple of thin
      ## Save the current_state
      theta[g,1,] <- lat0_curr_orig
      theta[g,2,] <- lon0_curr_orig
      theta[g,3,] <- d0_curr_orig
      theta[g,4,] <- t0_curr_orig
      theta[g,5,] <- alpha_curr_orig
      theta[g,6,] <- p_curr_orig
      theta[g,7,] <- vP_curr_orig
      theta[g,8,] <- tau_curr_orig
      theta[g,9,] <- invlambda_curr_orig
      
      g <- g+1
    }
    if((iter/iterations * 100) %in% seq(10,100,by = 10)){
      cat("### Progress:", iter/iterations * 100, " % \n")
    }
  }
  
  cat("I accepted the ", acc_1/iterations*100,"% of the proposed transition for lat0, lon0, d0, t0 \n")
  cat("I accepted the ", acc_2/iterations*100,"% of the proposed transition for alpha \n")
  cat("I accepted the ", acc_3/iterations*100,"% of the proposed transition for p \n")
  cat("I accepted the ", acc_4/iterations*100,"% of the proposed transition for vP \n")
  cat("I accepted the ", acc_5/iterations*100,"% of the proposed transition for tau \n")
  cat("I accepted the ", acc_6/iterations*100,"% of the proposed transition for invlambda \n")
  cat("I accepted the ", acc_sw/N_sw*100,"% of the proposed switches \n")
  
  return(theta)
}