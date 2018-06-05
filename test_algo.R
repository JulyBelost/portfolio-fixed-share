price_ratio_5 <- read.table("Ru5_data.csv", sep = ',', header = TRUE)
herst_exp_5 <- read.table("Ru5_hexp.csv", sep = ',', header = TRUE)
price_ratio_11 <- read.table("Ru11_data.csv", sep = ',', header = TRUE)
herst_exp_11 <- read.table("Ru11_hexp.csv", sep = ',', header = TRUE)
a = c(0.5, 0.7)
b = c(0.1, 0.25)

ht_to_pt = function(a, b, h_df){
  xi = function(a){ c(0, a-0.1, a, a+0.1, 1)}
  yi = function(b){ c(0, b, 0.5, 1-b, 1)}
  
  trust_level = sapply(h_df, pchipfun(xi(a), yi(b)))
  return(trust_level)
}

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
  
  K = cumprod(X_t) # portfolio wealth
  
  return(K)
}

run_main = function(price_ratio, herst_exp){ 
  K = run_pfs(price_ratio, ht_to_pt(a[1],b[1], herst_exp))
  K_z = run_pfs(price_ratio, (herst_exp*0)+1)
  K_n = rowMeans(sapply(price_ratio, cumprod))
  
  #par(mfrow = c(1, 2))
  plot(K, pch = "*")
  points(K_z, pch = "*", col = "blue")
  points(K_n, pch = "*", col = "red")
  
  print(sum(K>K_n)/length(K))
  print(sum(K>K_z)/length(K))
  print(sum(K_z>K_n)/length(K))
  print(head(K))
  print(tail(K))
  print(tail(K_n))
}

run_main(price_ratio_11, herst_exp_11)
run_main(price_ratio_5, herst_exp_5)
