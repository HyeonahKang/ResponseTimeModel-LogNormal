##########################################
##   Example Code for Joint Estimation  ##
##########################################
## H.-A. Kang (hkang@austin.utexas.edu)
## Last Updated: 05/21/2021
# rm(list=ls())

## Packages needed
# library(MASS)
# library(mvtnorm)
# library(LNIRT) # for comparison


#########################################
## Set hyper-parameters

hpar <- list(
  mu_a = 1,   var_a = 0.3^2,
  mu_b = 0,   var_b = 1,
  mu_cstr = -1.386, var_cstr = 0.04,
  mu_alp = 1, var_alp = 0.3^2,
  mu_bet = 0, var_bet = 1,
  rho_i = 0.3,
  mu_p = c(0, 0),
  rho_p = 0.3,
  trunc_i = matrix(c(0.3, 3, -4, 4, 0.0, 0.3), byrow=T, ncol=2, nrow=3))

## Item hyperparameters on the transformed scales
## astr = log a; alpstr = log alpha; cstr = logit c
## transformation applied to achieve normal prior & to ensure appropriate parameter ranges
hpar$mu_astr <- log(hpar$mu_a^2 / sqrt(hpar$var_a + hpar$mu_a^2))
hpar$var_astr <- log(hpar$var_a / hpar$mu_a^2 + 1)
hpar$mu_alpstr <- log(hpar$mu_alp^2 / sqrt(hpar$var_alp + hpar$mu_alp^2))
hpar$var_alpstr <- log(hpar$var_alp / hpar$mu_alp^2 + 1)

## prior mean of the transformed item parameters
hpar$mu_is <- unlist(hpar[c('mu_astr', 'mu_b', 'mu_cstr', 'mu_alpstr', 'mu_bet')])
nipar <- length(hpar$mu_is)

## prior covariance matrix of the transformed item parameters
var_is <- unlist(hpar[c('var_astr', 'var_b', 'var_cstr', 'var_alpstr', 'var_bet')])
cor_is <- matrix(hpar$rho_i, nipar, nipar); diag(cor_is) <- 1
hpar$cov_is <- sqrt(as.matrix(var_is) %*% t(as.matrix(var_is))) * cor_is
rm(var_is, cor_is)


#########################################
## Generate item parameters

seednum <- 12345
set.seed(seednum)

nitem <- 20
while (1 > 0){
  istr <- mvtnorm::rmvnorm(nitem, mean=hpar$mu_is, sigma=hpar$cov_is)
  ## transform back to the original scale
  ipar <- cbind(exp(istr[,1]), istr[,2], 1/(1 + exp(-istr[,3])),
                exp(istr[,4]), istr[,5])
  
  ## Check if the parameters are within the specified domains
  if (min(ipar[,1]) >= hpar$trunc_i[1,1] && max(ipar[,1]) <= hpar$trunc_i[1,2] &&
      max(abs(ipar[,2])) <= hpar$trunc_i[2,2] && max(abs(ipar[,5])) <= hpar$trunc_i[2,2] &&
      min(ipar[,3]) >= hpar$trunc_i[3,1] && max(ipar[,3]) <= hpar$trunc_i[3,2] &&
      min(ipar[,4]) >= hpar$trunc_i[1,1] && max(ipar[,4]) <= hpar$trunc_i[1,2] ) 
    break;
}
colnames(ipar) <- c("a", "b", "c", "alp", "bet")
rm(istr)

## Fix c at 0 ----
cfix <- 0
hpar$mu_is <- hpar$mu_is[-3]
hpar$cov_is <- hpar$cov_is[-3,]
hpar$cov_is <- hpar$cov_is[,-3]

while (1 > 0){
  istr <- mvtnorm::rmvnorm(nitem, mean=hpar$mu_is, sigma=hpar$cov_is)
  ## transform back to the original scale
  ipar <- cbind(exp(istr[,1]), istr[,2], cfix,
                exp(istr[,3]), istr[,4])
  
  ## Check if the parameters are within the specified domains
  if (min(ipar[,1]) >= hpar$trunc_i[1,1] && max(ipar[,1]) <= hpar$trunc_i[1,2] &&
      max(abs(ipar[,2])) <= hpar$trunc_i[2,2] && max(abs(ipar[,5])) <= hpar$trunc_i[2,2] &&
      # min(ipar[,3]) >= hpar$trunc_i[3,1] && max(ipar[,3]) <= hpar$trunc_i[3,2] &&
      min(ipar[,4]) >= hpar$trunc_i[1,1] && max(ipar[,4]) <= hpar$trunc_i[1,2] ) 
    break;
}
colnames(ipar) <- c("a", "b", "c", "alp", "bet")
rm(istr)


#########################################
## Generate person parameters

nexaminee <- 1000
hpar$cov_p <- matrix(c(1, hpar$rho_p, hpar$rho_p, 1), 2, 2)
while (1 > 0) {
  ppar <- mvtnorm::rmvnorm(nexaminee, mean=hpar$mu_p, sigma=hpar$cov_p)
  if (min(ppar)>=-4 && max(ppar)<=4){break}
}
colnames(ppar) <- c("th", "tau")


#########################################
## Generate data

gendata <- function(ppar, ipar, D=1.702, seednum){
  nperson <- nrow(ppar)
  nitem <- nrow(ipar)
  
  irf <- function(th,xi){
    p <- xi["c"] + (1-xi["c"]) / (1 + exp(-D * xi["a"] * (th - xi["b"]) ))
    return(p)
  }
  
  resp <- rt <- matrix(NA, nperson, nitem)
  for (i in 1:nperson){
    for (j in 1:nitem){
      resp[i,j] <- ( runif(1) < irf(ppar[i,"th"], ipar[j,1:3]) ) * 1
      rt[i,j] <- rlnorm(1, ipar[j,"bet"]-ppar[i,"tau"], 1/ipar[j,"alp"])
    }
  }
  
  return(list(resp=resp, rt=rt))
}

dat <- gendata(ppar, ipar, seednum=12345)
# hist(rowSums(dat$resp))
# hist(rowMeans(log(dat$rt)))


#########################################
## Estimate item parameters

source("mmap.R")
source("mmap_fsb.R")

mod <- mmap(dat$resp, dat$rt, 
            ppar_prior=list(mu_p=hpar$mu_p, cov_p=hpar$cov_p), 
            iparst_prior=list(mu_is=hpar$mu_is, cov_is=hpar$cov_is))

## Evaluate accuracy
(rmse <- sqrt(colMeans((mod$iest - ipar)^2)))


#########################################
## Estimate person parameters

source("map.R")

pout <- map(dat$resp, dat$rt, mod$iest, 
            ppar_prior=list(mu_p=hpar$mu_p, cov_p=hpar$cov_p),SE=T)

sqrt(colMeans((pout$est - ppar)^2))

compare <- par(mfrow=c(1, 2))
plot(ppar[,"th"], map_est$est[,"th"]); abline(0, 1)
plot(ppar[,"tau"], map_est$est[,"tau"]); abline(0, 1)
par(compare)

