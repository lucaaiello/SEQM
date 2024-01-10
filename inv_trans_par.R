inv_trans_par <-  function(x,lb,ub){
  
  L <- length(x)
  
  out <- vector(length = L)
  
  for (l in 1:L) {
    if(lb == 0 & ub == 1){
      if(x[l] < -745){
        out[l] <- (lb + ub*exp(-745))/(1+exp(-745))
      }else if(x[l] > 16.82){
        out[l] <- ub
      }else if(x[l] >= -745 & x[l] <= 16.82){
        out[l] <- (lb + ub*exp(x[l]))/(1+exp(x[l]))
      }
    }else if(lb == -90 & ub == 90){
      if(x[l] <= -17.4){
        out[l] <- lb
      }else if(x[l] >= 17.4){
        out[l] <- ub
      }else if(x[l] > -17.4 & x[l] < 17.4){
        out[l] <- (lb + ub*exp(x[l]))/(1+exp(x[l]))
      }
    }else if(lb == -180 & ub == 180){
      if(x[l] <= -15.79){
        out[l] <- lb
      }else if(x[l] >= 15.79){
        out[l] <- ub
      }else if(x[l] > -15.79 & x[l] < 15.79){
        out[l] <- (lb + ub*exp(x[l]))/(1+exp(x[l]))
      }
    }else if(lb==0 & ub == 3.5){
      if(x[l] <= -745.2){
        out[l] <- lb
      }else if(x[l] >= 15.8){
        out[l] <- ub
      }else if(x[l] > -745.2 & x[l] < 15.8){
        out[l] <- (lb + ub*exp(x[l]))/(1+exp(x[l]))
      }
    }else if(lb==800 & ub == 80000){
      if(x[l] <= -21.2){
        out[l] <- lb
      }else if(x[l] >= 16.6){
        out[l] <- ub
      }else if(x[l] > -21.2 & x[l] < 16.6){
        out[l] <- (lb + ub*exp(x[l]))/(1+exp(x[l]))
      }
    } else if(lb==0 & ub==120){
      if(x[l] <= -745.13){
        out[l] <- lb
      }else if(x[l] >= 14.7){
        out[l] <- ub
      }else if(x[l] > -745.13 & x[l] < 14.7){
        out[l] <- (lb + ub*exp(x[l]))/(1+exp(x[l]))
      }
    }
  }
  return(out)
}
