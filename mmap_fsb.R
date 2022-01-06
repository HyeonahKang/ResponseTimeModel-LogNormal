#' Fisher's Scoring for Joint Model MMAP
#'
#' The function computes gradient and Hessian
#' @param par_curr The estimates of item parameters for current the iteration
#' @keywords response time
#' @export


mmap_fsb <- function(nodes, f_kl, r_kl, g_kl, h_kl, 
                     par_curr, mu_is, cov_is_inv, D=1.702){  
  
  a <- par_curr[1]; b <- par_curr[2]; c <- par_curr[3];
  alp <- par_curr[4]; bet <- par_curr[5]
  
  P <- c + (1 - c) / (1 + exp(-D*a*(nodes[,1] - b)))
  P_star <- 1 / (1 + exp(-D*a*(nodes[,1] - b)))
  W <- P_star * (1 - P_star) / (P * (1 - P))  
  W[is.infinite(W)] <- 1  
  
  par_curr_str <- c(log(par_curr[1]), par_curr[2], log(par_curr[3]/(1 - par_curr[3])),
                    log(par_curr[4]), par_curr[5] )
                    
  
  ## Bayesian Marginal Estimation (Kim & Baker, ch.7.4.1)
  L1 <- D * a * (1-c) * sum((nodes[,1] - b) * W * (r_kl - f_kl * P) ) - sum(cov_is_inv[,1] * (par_curr_str - mu_is))
  
  L2 <- - D * a * (1-c) * sum(W * (r_kl - f_kl * P)) - sum(cov_is_inv[,2] * (par_curr_str - mu_is))
  
  L3 <- c * sum((r_kl - f_kl * P) / P) - sum(cov_is_inv[,3] * (par_curr_str - mu_is))
  
  L4 <- sum(f_kl - alp^2 * (g_kl - 2*(bet - nodes[,2]) * h_kl + (bet - nodes[,2])^2 * f_kl)) - sum(cov_is_inv[,4] * (par_curr_str - mu_is))
  
  L5 <- alp^2 * sum(h_kl - (bet - nodes[,2]) * f_kl) - sum(cov_is_inv[,5] * (par_curr_str - mu_is))
  
  
  ## Expected values of the Hessian matrix => Fisher's scoring
  L11 <- - D^2 * a^2 * sum(f_kl * (nodes[,1] - b)^2 * ((P - c) / (1 - c))^2 * (1 - P) / P)
  
  L22 <- - D^2 * a^2 * sum(f_kl * ((P - c) / (1 - c))^2 * (1 - P) / P)
  
  L33 <- - c^2 * sum(f_kl * (1 - P) / P)
  
  L44 <- - 2 * sum(f_kl)
  
  L55 <- -alp^2 * sum(f_kl)
  
  L12 <- D^2 * a^2 * sum(f_kl * (nodes[,1] - b) * ((P - c) / (1 - c))^2 * ((1 - P)/P))
  
  L13 <- - D * a * (c/(1 - c)) * sum(f_kl * (nodes[,1] - b) * (P-c) * (1-P)/P)
  
  L23 <- D * a * (c/(1 - c)) * sum(f_kl * (P - c) * (1 - P)/P)
  
  H <- matrix(c(L11, L12, L13, 0, 0, 
                L12, L22, L23, 0, 0, 
                L13, L23, L33, 0, 0, 
                0,    0,  0, L44, 0, 
                0,    0,  0,  0, L55), nrow=5, ncol=5, byrow=T) - cov_is_inv
  
  G <- c(L1, L2, L3, L4, L5); 
  names(G) <- colnames(H) <- rownames(H) <- c("a", "b", "c", "alp", "bet")
  
  return(list(G=G, H=H))
  
}


