
# THESIS PROJECT FINAL VERSION
# Ankith Sirnapally(25253247)
# Poisson Mixture EM vs Anscombe + mclust
# Extension to Zero-Inflated Poisson (ZIP)

library(mclust)
library(matrixStats)
library(e1071)

library(clue)
library(microbenchmark)
library(flexmix)

# ZIP DENSITY FUNCTION
dzip <- function(x, lambda, omega = 0, log = FALSE) {
   lp <- x
   lp[] <- -Inf
   id0 <- x == 0
   id1 <- !id0
  
  if(any(id0)) {lp[id0] <- matrixStats::logSumExp(c(log(omega), log1p(-omega) - lambda)) }
  if(any(id1)) { lp[id1] <- log1p(-omega) + dpois(x[id1], lambda, log = TRUE)}
  if(log) { return(lp)} 
   else {return(exp(lp))}
}
# RANDOM GENERATION FROM ZIP
rzip <- function(n, lambda, omega) {is_zero <- rbinom(n, 1, prob = omega)
   x <- numeric(n)
   x[!is_zero] <- rpois(sum(!is_zero), lambda)
  return(x)
}
# SIMULATE ZIP MIXTURE DATA
simulate_zip_mix <- function(n, pi, lambda, omega) {
  
  G <- length(pi)
  p <- nrow(lambda)
  z <- sample(1:G, size = n, replace = TRUE, prob = pi)
  X <- matrix(0, nrow = n, ncol = p)
  for(i in 1:n){g <- z[i]
  for(j in 1:p){
       X[i, j] <- rzip( n = 1, lambda = lambda[j, g], omega = omega[j, g])
    }
  }
  
  return(list(X = X, z = z))
}
# POISSON MIXTURE EM

pois_mix_EM <- function(X, G, maxit = 500, tol = 1e-6){
    n <- nrow(X)
    p <- ncol(X)
   Z <- matrix(runif(n * G), nrow = n, ncol = G)
   Z <- Z / rowSums(Z)
   pi_hat <- colMeans(Z)
  
  lambda_hat <- sweep(crossprod(Z, X), 1, n * pi_hat, FUN = "/")
  
  loglik_old <- -Inf
  
  for(iter in 1:maxit){ log_Z <- matrix(NA, nrow = n, ncol = G)
  for(g in 1:G){
   log_Z[, g] <- log(pi_hat[g]) + rowSums( dpois(X, lambda = matrix( lambda_hat[g, ],nrow = n, ncol = p, byrow = TRUE), log = TRUE )
        )
    }
    
    lse <- rowLogSumExps(log_Z)
    
    z <- exp(log_Z - lse)
    
    ll <- sum(lse)
    
    pi_hat <- colMeans(z)
    
    lambda_hat <- sweep(crossprod(z, X), 1, n * pi_hat, FUN = "/")
    
    if(abs(ll - loglik_old) < tol){
      break
    }
    
    loglik_old <- ll
  }
  
  return(list(pi = pi_hat,lambda = lambda_hat,z = z,classification = max.col(z),loglik = ll,iter = iter))
}
# UPDATED ZIP MIXTURE EM
zip_mix_EM <- function(X, G, omega_init = 0.2, maxit = 300, tol = 1e-6){
  
  n <- nrow(X)
  p <- ncol(X)
  
  Z <- matrix(runif(n * G), nrow = n, ncol = G)
  
  Z <- Z / rowSums(Z)
  
  pi_hat <- colMeans(Z)
  
  lambda_hat <- sweep(crossprod(Z, X), 1, n * pi_hat, FUN = "/")
  
  omega_hat <- matrix(omega_init, nrow = G, ncol = p)
  
  loglik_old <- -Inf
  
  for(iter in 1:maxit){ log_Z <- matrix(NA, nrow = n, ncol = G)
   for(g in 1:G){ temp <- matrix(0, nrow = n, ncol = p)
   for(j in 1:p){  temp[, j] <- dzip(   X[, j],lambda = lambda_hat[g, j], omega = omega_hat[g, j],log = TRUE )
      }
       log_Z[, g] <- log(pi_hat[g]) + rowSums(temp)
    }
    
    lse <- rowLogSumExps(log_Z)
    
    z <- exp(log_Z - lse)
    
    ll <- sum(lse)
  # EXTRA E-STEP FOR U-HAT
    
    u_hat <- array(0, dim = c(n, p, G))
    
    for(g in 1:G){
      
      for(j in 1:p){
        
        zero_id <- X[,j] == 0
        
        u_hat[zero_id, j, g] <-
          omega_hat[g,j] /
          (
            omega_hat[g,j] +
              (1 - omega_hat[g,j]) * exp(-lambda_hat[g,j])
          )
      }
    }
    
    # M-STEP
    
    pi_hat <- colMeans(z)
    
    for(g in 1:G){
      
      for(j in 1:p){
        
        numerator <- sum(z[,g] *  (1 - u_hat[,j,g]) * X[,j] )
       denominator <- sum(z[,g] * (1 - u_hat[,j,g]) )
         lambda_hat[g,j] <- numerator / denominator
        
      }
    }
    
    for(g in 1:G){
      for(j in 1:p){
         omega_hat[g,j] <- sum(z[,g] * u_hat[,j,g]) / sum(z[,g])
      }
    }
    
    if(abs(ll - loglik_old) < tol){
      break
    }
    
    loglik_old <- ll
  }
  
  return(list( pi = pi_hat,lambda = lambda_hat, omega = omega_hat, z = z, classification = max.col(z),loglik = ll,iter = iter))
}
# TRUE PARAMETERS
set.seed(123)
n <- 300
G <- 2
p <- 5
pi_true <- c(0.4, 0.6)
lambda_true <- matrix(c(2,7,3,8,4,9,2,6,5,10),nrow = p,ncol = G,byrow = TRUE)
omega_true <- matrix(
  c(0.4,0.1, 0.3,0.2, 0.2,0.1, 0.3,0.2,0.4,0.1 ), nrow = p, ncol = G, byrow = TRUE)

# GENERATE DATA
sim <- simulate_zip_mix(n = n,pi = pi_true,lambda = lambda_true,omega = omega_true)
X <- sim$X
true_z <- sim$z

# VISUALIZE DATA
hist(X[,1], breaks = 20,main = "Zero-Inflated Poisson Mixture Data",xlab = "Counts")

# ZERO INFLATION INDEX
ZI <- apply(X, 2, function(x){lambda_hat <- mean(x); mean(x == 0) / exp(-lambda_hat)})
print(ZI)
# ANSCOMBE + MCLUST
start_ans <- Sys.time()
Y <- 2 * sqrt(X + 3/8)
fit_ans <- Mclust( Y, G = G, verbose = FALSE)
end_ans <- Sys.time()
time_ans <- end_ans - start_ans
ARI_ans <- adjustedRandIndex( true_z, fit_ans$classification)

# POISSON EM

start_pois <- Sys.time()
fit_pois <- pois_mix_EM(X, G = G)
end_pois <- Sys.time()
time_pois <- end_pois - start_pois
ARI_pois <- adjustedRandIndex(true_z,fit_pois$classification)

# UPDATED ZIP EM

start_zip <- Sys.time()
fit_zip <- zip_mix_EM(X, G = G)
cat("\nTrue Lambda:\n")
print(lambda_true)

cat("\nEstimated Lambda:\n")
print(round(fit_zip$lambda, 3))

cat("\nTrue Omega:\n")
print(omega_true)

cat("\nEstimated Omega:\n")
print(round(fit_zip$omega, 3))
end_zip <- Sys.time()

time_zip <- end_zip - start_zip

ARI_zip <- adjustedRandIndex(true_z,fit_zip$classification)

# FLEXMIX

df <- data.frame(X1 = X[,1])
fit_flex <- flexmix(X1 ~ 1,data = df,k = G,model = FLXMRglm(family = "poisson"))
ARI_flex <- adjustedRandIndex(true_z,clusters(fit_flex))

# RESULTS TABLE

results <- data.frame( Method = c( "Anscombe + mclust", "Poisson EM", "ZIP EM", "flexmix"),
 ARI = c(ARI_ans,ARI_pois,ARI_zip,ARI_flex),
  Time = c( as.numeric(time_ans), as.numeric(time_pois), as.numeric(time_zip), NA
  )
)

print(results)

# SIMPLE SIMULATION STUDY

settings <- expand.grid(n = c(100, 300),G = c(2, 3),p = c(5, 10))
simulation_results <- data.frame()
for(i in 1:nrow(settings)){
  n_sim <- settings$n[i]
  G_sim <- settings$G[i]
  p_sim <- settings$p[i]
  pi_true <- rep(1/G_sim, G_sim)
  lambda_true <- matrix(sample(2:10, p_sim * G_sim, replace = TRUE),nrow = p_sim, ncol = G_sim)
  omega_true <- matrix(
    runif(p_sim * G_sim, 0.1, 0.4),nrow = p_sim,ncol = G_sim)
  sim <- simulate_zip_mix(n = n_sim,pi = pi_true,lambda = lambda_true,omega = omega_true)
  X_sim <- sim$X
  z_true <- sim$z
  Y_sim <- 2 * sqrt(X_sim + 3/8)
  fit_ans <- Mclust(Y_sim, G = G_sim, verbose = FALSE)
  ARI_ans <- adjustedRandIndex(z_true,fit_ans$classification)
  fit_zip <- zip_mix_EM(X_sim,G = G_sim)
  ARI_zip <- adjustedRandIndex(z_true, fit_zip$classification)
  simulation_results <- rbind( simulation_results,
    data.frame( n = n_sim, G = G_sim, p = p_sim, ARI_Anscombe = ARI_ans, ARI_ZIP = ARI_zip
    )
  )
}

print(simulation_results)

# MICROBENCHMARK TIMING

mb <- microbenchmark(Anscombe_Mclust = Mclust( Y, G = G, verbose = FALSE),
  Poisson_EM = pois_mix_EM( X, G = G), ZIP_EM = zip_mix_EM( X, G = G ), times = 10)
print(mb)

# SUMMARY OUTPUT

cat("RESULTS SUMMARY\n")

cat("\nAdjusted Rand Index:\n")
print(results)

cat("\nZero Inflation Index:\n")
print(ZI)

cat("\nPoisson EM iterations:\n")
print(fit_pois$iter)

cat("\nZIP EM iterations:\n")
print(fit_zip$iter)

cat("\nEstimated ZIP omega values:\n")
print(round(fit_zip$omega, 3))

cat("\nSimulation Study Results:\n")
print(simulation_results)
