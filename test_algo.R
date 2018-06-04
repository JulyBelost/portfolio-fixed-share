price_ratio_5 <- read.table("Ru5_data.csv", sep = ',', header = TRUE)
herst_exp_5 <- read.table("Ru5_hexp.csv", sep = ',', header = TRUE)
# price_ratio_11 <- read.table("Ru11_data.csv", sep = ',', header = TRUE)
# herst_exp_11 <- read.table("Ru11_hexp.csv", sep = ',', header = TRUE)
a = c(0.5, 0.7)
b = c(0.1, 0.25)

xi = function(a){ c(0, a-0.1, a, a+0.1, 1)}
yi = function(b){ c(0, b, 0.5, 1-b, 1)}

herst_to_trust_level = pchipfun(xi(a[1]), yi(b[1]))

trust_level = sapply(herst_exp_5, herst_to_trust_level)

run_pfs = function(price_ratio, trust_level){
  # Initial parameters
  T <- nrow(price_ratio)      # number of steps
  N <- ncol(price_ratio)      # number of instruments
  
  w_ <- c()                      # w^*
  w <- rep(1/N, N)               # w
  w_m <- rep(1/N, N)             # w^m
  
  X_t <- c()
  
  # Algorithm
  for (t in 1:T){
    p_t = trust_level[t,]
    
    w_ = (p_t*w)/sum(p_t*w)
    
    x_t = price_ratio[t,]
    
    X_t[t] <- sum(x_t*w_)
    
    
    # LOSS UPDATE
    w_m = (w * (p_t*x_t + (1 - p_t)*X_t[t]))/X_t[t]
    
    # MIXING(Fixed-Share) UPDATE
    alpha = 1/t
    w = alpha/N + (1 - alpha)*w_m
    
  }
  
  K = cumprod(c(1,X_t)) # portfolio wealth
  
  return(K)
}

K = run_pfs(herst_exp_5, trust_level)
K_n = rowMeans(sapply(price_ratio_5, cumprod))
