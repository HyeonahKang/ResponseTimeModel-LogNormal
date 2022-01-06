#' Maximum a Posteriori of person parameters for the Joint Model of IRM and RTM
#'
#' The function estimates the person parameters of the joint model of the logistic
#' item-response model and log-normal response-time model using known item parameters
#' 
#' @param resp item response matrix
#' @param RT response time matrix
#' @param ipar item parameter matrix (nitem * 5)
#' @param D scaling constant that approximates the logistic model to normal-ogive model 
#' @param mu_p mean of latent ability and speedness; The default is set to 0 for both
#' @param cov_p covariance matrix for latent trait parameters; if not known, assume identity matrix (i.e., no collateral information)
#' @keywords MAP joint reponse time model IRT
#' @export
#' @examples
#' out <- map(resp, RT)

map <- function(resp, RT, ipar, 
                ppar_prior=list(mu_p=NULL, cov_p=NULL),
                SE=FALSE, D=1.702, 
                maxitr=200, tol=list(nr=1.0e-4, loglike=1.0e-6),
                trunc=c(-3.5, 3.5) ){
  
 
  ## #~#~#~#~#~#~#~# ##
  ##   Setting up    ##
  ## #~#~#~#~#~#~#~# ##
  
  nexaminee <- nrow(resp)
  nitem <- ncol(resp)
  
  ## Priors (if specified)
  for (v in 1:length(ppar_prior)){
    eval(parse(text=paste0(names(ppar_prior)[v], " <- ppar_prior[[", v, "]]") ))
  }
  
  ## if not specified
  if(is.null(mu_p)){mu_p <- c(0,0)}
  if(is.null(cov_p)){cov_p <- diag(2)}
  
  cov_p_inv <- solve(cov_p)
  
  ## Functions
  irf <- function(th, xi){
    # one person multiple items
    (p <- xi[,3] + (1 - xi[,3]) / ( 1 + exp(-D * xi[,1] * (th - xi[,2])) ))
  }
  
  rtf <- function(tau, time, xi){
    # one person multiple items
    (p <- xi[,1] / (time * sqrt(2*pi)) * exp(- 1/2 * xi[,1]^2 * (log(time) - (xi[,2] - tau))^2))
  }
  
  maxradius <- 2
  
  pest_init <- cbind(qnorm(rowMeans(resp)), -as.numeric(scale(rowMeans(log(RT)))))
  colnames(pest_init) <- c("th", "tau")
  pest_init[pest_init[,"th"] < trunc[1] - 0.5,1] <- trunc[1] - 0.5
  pest_init[pest_init[,"th"] > trunc[2] + 0.5,1] <- trunc[2] + 0.5
  pest_init[pest_init[,"tau"] < trunc[1] - 0.5,2] <- trunc[1] - 0.5
  pest_init[pest_init[,"tau"] > trunc[2] + 0.5,2] <- trunc[2] + 0.5
  
  
  pest <- matrix(NA, nexaminee, 2, dimnames=list(NULL, c("th", "tau")))
  nonconv <- rep(NA, nexaminee)
  
  
  for (i in 1:nexaminee) {# i<-3
    
    itr <- 0
    eta_curr <- pest_init[i,]
    
    track_nr <- list()    
    track_nr$loglike <- rep(NA, maxitr+1); track_nr$loglike[1] <- 100
    
    while (itr < maxitr){
      
      itr <- itr + 1
      
      eta_try <- eta_curr
      for (k in 1:2){
        eta_try[k] <- ifelse(eta_curr[k] <= trunc[1], trunc[1], eta_curr[k])
        eta_try[k] <- ifelse(eta_curr[k] >= trunc[2], trunc[2], eta_curr[k])
      }
      
      ## first and second partial derivatives
      pr <- irf(eta_try[1], ipar)
      
      l1 <- D * sum(ipar[,1] * (resp[i,] - pr) * (pr - ipar[,3]) / (pr * (1 - ipar[,3]))) - 
        sum(cov_p_inv[,1] * t(eta_try - mu_p) )
      
      l2 <- - sum( ipar[,4]^2 * (log(RT[i,]) - (ipar[,5] - eta_try[2]))) - 
        sum(cov_p_inv[,2] * t(eta_try - mu_p))
      
      ## Fisher's scoring
      lam11 <- D^2 * sum(ipar[,1]^2 * (1 - pr) * (pr - ipar[,3]) * (ipar[,3] - pr) /
                           (pr * (1 - ipar[,3])^2 )) - cov_p_inv[1,1]
      lam22 <- - sum(ipar[,4]^2) - cov_p_inv[2,2]
      lam12 <- - cov_p_inv[1,2]
      
      ## hessian
      hess <- matrix(c(lam11, lam12, lam12, lam22), 2, 2, byrow=T)      
      delta <- c(l1, l2) %*% solve(hess)
      
      ## Normalize delta if it exceeds maximum radius.
      if (sqrt(sum(delta^2)) > maxradius){
        delta <- delta * maxradius / abs(sqrt(sum(delta^2)))
      }
            
      track_nr$loglike[itr] <- log( prod ( pr^resp[i,] * (1-pr)^(1-resp[i,]) * rtf(eta_try[2], RT[i,], ipar[,4:5]) ) )
      
      if ( ( max(abs(delta)) < tol$nr ) && 
           (itr > 1 && abs( diff(track_nr$loglike[(itr-1):(itr)])) < tol$loglike ) ){
        eta_hat <- eta_try
        break;
      }
      
      eta_curr <- eta_try - delta
      
    } # end of while   


    if (itr <= maxitr){
      ## If converged within the tolerance criterion
      
      if (any(abs(eta_hat) > (trunc[2] + 0.5))){
        nonconv[i] <- 2 # out of bounds
        eta_hat[eta_hat < (trunc[1] - 0.5)] <- trunc[1] - 0.5
        eta_hat[eta_hat > (trunc[2] + 0.5)] <- trunc[2] + 0.5
      }
      
      pest[i,] <- eta_hat
      
    } else {
      
      ## If not converged
      nonconv[i] <- 1
      eta_curr[eta_curr < (trunc[1] - 0.5)] <- trunc[1] - 0.5
      eta_curr[eta_curr > (trunc[2] + 0.5)] <- trunc[2] + 0.5
      pest[i,] <- eta_curr
    }
    
  } # end of i (examinee)
  
  
  if (isTRUE(SE)){
    
    pest_se <- matrix(NA, nexaminee, 2, dimnames=list(NULL, c("th", "tau")))
    for (i in 1:nexaminee){# i <- 1
      
      pr <- irf(pest[i,"th"], ipar)
      lam11 <- D^2 * sum(ipar[,1]^2 * (1 - pr) * (pr - ipar[,3]) * (ipar[,3] - pr) /
                           (pr * (1 - ipar[,3])^2 )) 
      lam22 <- - sum(ipar[,4]^2)
      lam12 <- 0 
      
      hess <- matrix(c(lam11, lam12, lam12, lam22), 2, 2, byrow=T)
      pest_se[i,] <- diag(sqrt(-solve(hess)))
      
    }
    
  } else {
    
    pest_se <- NULL
    
  }
  
  ### #~#~#~#~#~#~#~#~#~# ###
  ###    Export results   ###
  ### #~#~#~#~#~#~#~#~#~# ###
  
  output <- list(est=pest, 
                 se=pest_se,
                 nonconv=nonconv,
                 init=pest_init,
                 track=track_nr)
  
  return(output)
}
