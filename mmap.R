#' Marginal Maximum a Posteriori for the Joint Model of Responses and RTs
#'
#' The function estimates item parameters of the joint model of the logistic item-response
#' model and log-normal response-time model via marginal maximum a posteriori
#' 
#' @param resp binary item response matrix
#' @param RT response time matrix (assumed to be scaled so that underlying tau is fixed at 0)
#' @param mu_p mean of latent ability and speedness; The default is set to 0 for both
#' @param cov_p covariance matrix for latent trait parameters; if not known, assume identity matrix (i.e., no collateral information)
#' @param mu_is mean of transformed item parameters (log a, b, logit c, log alpha, beta)
#' @param cov_is covariance matrix for transformed item parameters
#' @param ipar_init user-defined initial item parameter values
#' @param ipar_fix  user-defined fixed item parameter matrix; column = (a, b, c, alpha, beta)
#' @param D scaling constant that approximates the logistic model to normal-ogive model
#' @param tol Stopping criteria; max_em=maximum em iteration; 
#'            tol_em=difference in the item parameter values between EM cycles
#'            loglike=difference in the log-likelihood between EM cycles
#'            max_nr=maximum Newton Raphson iteration
#'            tol_nr=difference in the item parameters between NR iterations
#' @keywords response time IRT joint model 
#' @export
#' @examples
#' results_MBME <- mbme_jt(resp=resp, RT=RT)

mmap <- function(resp, RT, 
                 ppar_prior=list(mu_p=NULL, cov_p=NULL), 
                 iparst_prior=list(mu_is=NULL, cov_is=NULL), 
                 ipar_init=NULL, # ipar_fix = NULL, 
                 D=1.702, 
                 tol=list(max_em=1000, tol_em=0.001, loglike=1e-4, max_nr=500, tol_nr=0.001)){
  
  ## MMAP Estimation for the Joint Model or Item Responses & Response Times
  ## Required packages::
  # library(MASS)
  # library(mvtnorm)
  
  ## H.A. Kang (hkang@austin.utexas.edu)
  ## Last Updated: 05/21/2021
  
 
  ### #~#~#~#~#~#~#~#~#~#~# ###
  ###   Set up variables    ###
  ### #~#~#~#~#~#~#~#~#~#~# ###
  
  nexaminee <- nrow(resp); nitem <- ncol(resp)
  ipar_name <- c("a", "b", "c", "alp", "bet")
  nipar <- length(ipar_name)
  
  ## Priors (if specified) ------------------
  for (v in 1:length(ppar_prior)){
    eval(parse(text=paste0(names(ppar_prior)[v], " <- ppar_prior[[", v, "]]") ))
  }
  
  for (v in 1:length(iparst_prior)){
    eval(parse(text=paste0(names(iparst_prior)[v], " <- iparst_prior[[", v, "]]") ))
  }
  
  for (v in 1:length(tol)){
    eval(parse(text=paste0(names(tol)[v], " <- tol[[", v, "]]") ))
  }
  
  if(is.null(mu_p)){mu_p <- c(0,0)}
  if(is.null(cov_p)){cov_p <- diag(2)}
  if(is.null(mu_is)){
    mu_a <- 1; var_a <- 0.3^2; mu_alp <- 1; var_alp <- 0.3^2; 
    mu_b <- 0; var_b <- 1; mu_bet <- 0; var_bet <- 1; 
    mu_cstr <- -1.386; var_cstr <- 0.04; 
    mu_astr <- log(mu_a^2/sqrt(var_a + mu_a^2))
    mu_alpstr <- log(mu_alp^2/sqrt(var_alp + mu_alp^2))
    mu_i <- c(mu_astr, mu_b, mu_cstr, mu_alpstr, mu_bet)
  }
  
  # if(is.null(ipar_fix)){item_par_fix <- matrix(NA, ncol=5, nrow=nitem)}
  
  ipar_trunc <- matrix(c(0.2, 3, -4, 4, 0.001, 0.5, 
                         0.2, 3, -4, 4), nrow=5, ncol=2, byrow=T)
  rownames(ipar_trunc) <- ipar_name 
  colnames(ipar_trunc) <- c("lb", "ub")
  
  
  ## Functions ---------------
  irf <- function(th, xi){
    # one person multiple items
    (p <- xi[,3] + (1 - xi[,3]) / ( 1 + exp(-D * xi[,1] * (th - xi[,2])) ))
  }
  
  rtf <- function(tau, time, xi){
    # one person multiple items
    (p <- xi[,1] / (time * sqrt(2*pi)) * exp(- 1/2 * xi[,1]^2 * (log(time) - (xi[,2] - tau))^2))
  }
  
  ## Incidental Variables Quadrature
  nquad <- length(seq(-3.5,3.5,0.5))
  nnode <- nquad^2
  nodes <- matrix(0, nrow=nnode, ncol=2)
  nodes[,1] <- rep(seq(-3.5,3.5,0.5), each = nquad)
  nodes[,2] <- rep(seq(-3.5,3.5,0.5), times = nquad)
  A_kl <- mvtnorm::dmvnorm(nodes, mu_p, cov_p)
  Amat <- matrix(rep(A_kl, nexaminee), nnode, nexaminee) # View(Amat)
  
  start_time <- Sys.time()
  
  ### #~#~#~#~#~#~#~# ###
  ###   Initialize    ###
  ### #~#~#~#~#~#~#~# ###
  
  if (is.null(ipar_init)){
    
    ipar_init <- matrix(NA, nitem, nipar, dimnames=list(NULL, ipar_name))
    total_score <- rowSums(resp)
    pidx_low <- subset(1:nexaminee, total_score < quantile(total_score, prob=.05, na.rm=TRUE))
    
    for (j in 1:nitem){#
      
      ## i) disc: transform point biserial correlation coefficient (Eq 1.9 in Cohen & Kim)
      mean_y1 <- mean(total_score[resp[,j] == 1])
      px <- sum(resp[,j])/nexaminee
      pbis <- (mean_y1 - mean(total_score)) / sd(total_score) * sqrt(px/(1-px))
      ipar_init[j,"a"] <- D * sqrt(pbis^2/(1 - pbis^2))
      
      ## ii) diff: sample prop correct
      ipar_init[j,"b"] <- qnorm(1- sum(resp[,j])/nexaminee)
      if (ipar_init[j,"b"] < ipar_trunc["b","lb"]){
        ipar_init[j,"b"] <- ipar_trunc["b","lb"]
      } else if (ipar_init[j,"b"] > ipar_trunc["b","ub"]){
        ipar_init[j,"b"] <- ipar_trunc["b","ub"]
      }
      
      ## iii) low-asymptote
      ipar_init[j, "c"] <- min(c(mean(resp[pidx_low, j]), ipar_trunc["c", "ub"]))
      if (ipar_init[j,"c"]==0){
        ipar_init[j,"c"] <- ipar_trunc["c","lb"]
      }
      
      ## iv) time-disc: correlation
      ipar_init[j,"alp"] <- 1/sd(log(RT[,j]))
      ## Note. Bias expected when approximating empirically
      
      ## v) diff: sample mean
      ipar_init[j,"bet"] <- mean(log(RT[,j]))
    }
       
  }
  
  ## initial item parameters when entering EM
  ipar_start <- ipar_init
  
  ## Item Covariance Matrix
  if (is.null(cov_is)) {
    ## Use the initial item parameters to obtain empirical covariance matrix
    tmp <- cbind(log(ipar_start[,"a"]),
                 ipar_start[,"b"],
                 log(ipar_start[,"c"]/(1-ipar_start[,"c"])),
                 log(ipar_start[,"alp"]),
                 ipar_start[,"bet"])
    tmp <- scale(tmp, scale=FALSE)
    cov_is <- (t(tmp) %*% tmp)/nitem
  }
  cov_is_inv <- solve(cov_is)
  
  
  ### #~#~#~#~# ###
  ###    EM     ###
  ### #~#~#~#~# ###
  
  cyc <- 1
  delta_em <- rep(100, nipar)
  max_radius <- 3 # safeguard in NR
  ipar_curr <- ipar_start  
  
  ## | Debug -------------
  track_em <- list()
  track_em$loglike <- rep(NA, max_em+1); track_em$loglike[1] <- 100
  track_em$delta <- matrix(NA, max_em+1, nipar, dimnames=list(NULL, ipar_name))  
  track_em$nr_conv <- matrix(NA, max_em, nitem)
  ##   ----------------- |
  
  while (cyc < max_em){
    
    # print(paste0("EM Cycle: ", cyc))
    
    ## Likelihood of observing data given the current parameter estimates
    pr_curr <- matrix(NA, nnode, nitem)
    for (k in 1:nnode){
      pr_curr[k, ] <- irf(nodes[k,1], ipar_curr[, 1:3])
    }
    
    L_resp <- matrix(NA, nnode, nexaminee) # likelihood of responses except for item j
    L_rt <- matrix(NA, nnode, nexaminee)   # likelihood of RTs except for item j
    for (i in 1:nexaminee){# i<-1
      L_resp[, i] <- apply( pr_curr ^ matrix(rep(resp[i,], each=nnode), nnode, nitem) *
                              (1 - pr_curr) ^  matrix(rep(1 - resp[i,], each=nnode), nnode, nitem), 1, prod )
      for (l in 1:nnode){# l<-1
        L_rt[l, i] <- prod( rtf(nodes[l,2], RT[i,], ipar_curr[, 4:5]) )
      }
    }
    L_curr <- L_resp * L_rt # conditional independence
    
    
    ### #~#~#~#~#~#~# ###
    ###     Update    ###
    ### #~#~#~#~#~#~# ###
    
    ## Estimate item parameters
    ipar_new <- matrix(NA, nitem, nipar)
    conv_nr <- rep(NA, nitem)
    
    denom <- colSums(L_curr * Amat)
    dmat <- matrix(rep(denom, each=nnode), nnode, nexaminee)
    f_kl <- rowSums(L_curr * Amat / dmat)
    track_em$loglike[cyc+1] <- sum(log(colSums(L_curr * Amat)))
    
    for (j in 1:nitem){# j<-1; j<-j+1
      
      ## Artificial data
      r_kl <- rowSums( matrix(rep(resp[,j], each=nnode), nnode, nexaminee) * L_curr * Amat / dmat )
      h_kl <- rowSums( matrix(rep(log(RT[,j]), each=nnode), nnode, nexaminee) * L_curr * Amat / dmat )
      g_kl <- rowSums( matrix(rep(log(RT[,j])^2, each=nnode), nnode, nexaminee) * L_curr * Amat / dmat )
      
      ## Fisher Scoring for each item
      itr <- 0
      delta <- rep(100, nipar)
      
      a_est   <- ipar_curr[j,1]; astr <- log(a_est)
      b_est   <- ipar_curr[j,2]
      c_est   <- ipar_curr[j,3]; cstr <- log(c_est / (1 - c_est))
      alp_est <- ipar_curr[j,4]; alpstr <- log(alp_est)
      bet_est <- ipar_curr[j,5]       
      
      while ((itr<max_nr) && (max(abs(delta))>tol_nr)) {
        
        itr <- itr + 1
        
        ## Compute Gradient & Hessian
        fsb <- mmap_fsb(nodes, f_kl, r_kl, g_kl, h_kl,
                        par_curr=c(a_est, b_est, c_est, alp_est, bet_est), 
                        mu_is, cov_is_inv) 
        
        delta <- (MASS::ginv(fsb$H)) %*% fsb$G
        
        ## "safe_radius": the maximum stepsize allowed in NR;
        if(norm(delta, type="2") > max_radius){
          delta <- delta * max_radius / norm(delta, type="2")
        }
        
        ## Update item parameters
        if (abs(ipar_start[j,"a"] - exp(astr - delta[1])) > 2 ||  
            abs(ipar_start[j,"b"] - (b_est - delta[2])) > 3 ||
            abs(ipar_start[j,"c"] - (1/(1 + exp(-(cstr- delta[3]) )))) > 0.4){
          
          if (abs( ipar_start[j, "a"] - exp(astr-delta[1]) ) > 2) {
            delta[1] <- delta[1] / 2          
          }
          if (abs( ipar_start[j, "b"] - (b_est - delta[2]) ) > 3){
            delta[2] <- delta[2] / 2
          }
          if ( abs( ipar_start[j, "c"] - (1 / (1 + exp( -(cstr-delta[3] ) ))) ) > 0.4){
            delta[3] <- delta[3] / 2
          }
          
        } # end of if
        
        astr    <- astr - delta[1]; a_est <- exp(astr) ;
        b_est   <- b_est - delta[2]
        cstr    <- cstr - delta[3]; c_est <- 1 / ( 1 + exp(-cstr) ) ;
        alpstr  <- alpstr - delta[4]; alp_est <- exp(alpstr) ;
        bet_est <- bet_est - delta[5]
               
      } # end of FS while
           
      ## NR convergence
      if ( is.nan(a_est) || is.nan(b_est) || is.nan(c_est)|| is.nan(alp_est)||is.nan(bet_est)||
           a_est < 0 || a_est > (2*ipar_trunc["a","ub"]) ||
           b_est < (2*ipar_trunc["b","lb"]) ||
           b_est > (2*ipar_trunc["b","ub"]) ||
           c_est < 0 || c_est > 1||
           alp_est < 0 || alp_est > (2*ipar_trunc["alp","ub"]) ||
           bet_est < (2* ipar_trunc["bet","lb"]) || bet_est > (2*ipar_trunc["b","ub"]) ){
        conv_nr[j] <- 1  ## Not converged
      } else {
        conv_nr[j] <- 0  ## Converged
      }
      
      if (itr==max_nr){conv_nr[j] <- 1}
      
      ipar_new[j,] <- c(a_est, b_est, c_est, alp_est, bet_est)
      # Debug) # track_em$nr_delta[j,,cyc] <- delta
      
    } # End of item optimization
    
    ## | Debug EM -------
    delta_em <- ipar_new - ipar_curr
    track_em$delta[cyc+1,] <- apply(abs(delta_em), 2, max)
    track_em$nr_conv[cyc+1,] <- conv_nr
    ##    -------------- |
    
    ## If not converged
    if ((cyc > 2) && ( (max(delta_em) < tol_em) || abs(diff(track_em$loglike[cyc:(cyc+1)])) < tol$loglike) ){
      ipar_est <- ipar_new
      break ;
    }
    
    ## If not converged
    cyc <- cyc + 1
    ipar_curr <- ipar_new
    
  } # end of EM cycle
  
  colnames(ipar_est) <- ipar_name
  
 
  
  ### #~#~#~#~#~#~#~#~#~# ###
  ###    Standard Error   ###
  ### #~#~#~#~#~#~#~#~#~# ###
  
  if (!exists("ipar_est")){ipar_est <- ipar_curr}
  
  ## Likelihood of observing data given the current parameter estimates
  pr <- matrix(NA, nnode, nitem)
  for (k in 1:nnode){
    pr[k, ] <- irf(nodes[k,1], ipar_est[, 1:3])
  }
  
  L_resp <- matrix(NA, nnode, nexaminee) # likelihood of responses except for item j
  L_rt <- matrix(NA, nnode, nexaminee)   # likelihood of RTs except for item j
  for (i in 1:nexaminee){# i<-1
    L_resp[, i] <- apply( pr ^ matrix(rep(resp[i,], each=nnode), nnode, nitem) *
                            (1 - pr) ^  matrix(rep(1 - resp[i,], each=nnode), nnode, nitem), 1, prod )
    for (l in 1:nnode){# l<-1
      L_rt[l, i] <- prod( rtf(nodes[l,2], RT[i,], ipar_est[, 4:5]) )
    }
  }
  Like <- L_resp * L_rt
  
  # Amat <- matrix(rep(A_kl, nexaminee), nnode, nexaminee)
  denom <- colSums(Like * Amat)
  dmat <- matrix(rep(denom, each=nnode), nnode, nexaminee)
  f_kl <- rowSums( Like * Amat / dmat )
  
  ipar_vcov <- ipar_vcov_str <- array(NA, c(nipar, nipar, nitem))
  ipar_se <- matrix(NA, nitem, nipar, dimnames=list(NULL, ipar_name))
  for (j in 1:nitem){# j<-1
    
    ## Artificial data
    r_kl <- rowSums( matrix(rep(resp[,j], each=nnode), nnode, nexaminee) * Like * Amat / dmat )
    h_kl <- rowSums( matrix(rep(log(RT[,j]), each=nnode), nnode, nexaminee) * Like * Amat / dmat )
    g_kl <- rowSums( matrix(rep(log(RT[,j])^2, each=nnode), nnode, nexaminee) * Like * Amat / dmat )
    
    ## Compute Gradient & Hessian
    fsb <- mmap_fsb(nodes, f_kl, r_kl, g_kl, h_kl,
                    par_curr=ipar_est[j,], 
                    mu_is, cov_is_inv) 
    
    ipar_vcov_str[,,j] <- solve(-(fsb$H + cov_is_inv))
    coeff <- matrix(c(ipar_est[j,"a"], 1, ipar_est[j,"c"] * (1 - ipar_est[j,"c"]), ipar_est[j, "alp"], 1), 5, 1)
    ipar_vcov[,,j] <- coeff %*% t(coeff) * ipar_vcov_str[,,j]
    ipar_se[j,] <- sqrt(diag(ipar_vcov[,,j]))
    
  }
  
  end_time <- Sys.time()
  est_time <- end_time - start_time
  
  ### #~#~#~#~#~#~#~#~#~# ###
  ###    Export results   ###
  ### #~#~#~#~#~#~#~#~#~# ###
  
  out <- list(iest=ipar_est, ise=ipar_se, icov=ipar_vcov, icov_str=ipar_vcov_str,
              track=list(cyc=cyc,
                         em=track_em,
                         est_time=est_time),
              supp=list(quad=nodes, 
                        ppar_prior=ppar_prior,
                        iparst_prior=iparst_prior,
                        tol=tol))  #, input=list(resp=resp, RT=RT)
  
  
  return(out)
  
}


