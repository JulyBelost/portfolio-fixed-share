# prepared price ratio and herst exponent tables for two sets of instruments
price_ratio_5 <- read.table("Ru5_data.csv", sep = ',', header = TRUE)
herst_exp_5 <- read.table("Ru5_hexp.csv", sep = ',', header = TRUE)
price_ratio_11 <- read.table("Ru11_data.csv", sep = ',', header = TRUE)
herst_exp_11 <- read.table("Ru11_hexp.csv", sep = ',', header = TRUE)

# algorithm parameters
a = c(0.5, 0.7)
b = c(0.1, 0.25)
alpha = c(0.001,0.01,0.1,0.25,1)

# herst exponent transformation into trust levels
ht_to_pt = function(a, b, herst_df){
  xi = function(a){ c(0, a-0.1, a, a+0.1, 1)}
  yi = function(b){ c(0, b, 0.5, 1-b, 1)}
  # plot(xi(a),yi(b))
  
  trust_level_df = sapply(herst_df, pchipfun(xi(a), yi(b)))
  return(trust_level_df)
}

# portfolio fixed share algorithm for unreliable instruments
run_portfolio_fs = function(price_ratio_df, trust_level_df, alpha){
  # Initial parameters
  T <- nrow(price_ratio_df)      # number of steps
  N <- ncol(price_ratio_df)      # number of instruments
  
  w_ <- c()                      # w^*
  w <- rep(1/N, N)               # w
  w_m <- rep(1/N, N)             # w^m
  
  X_t <- c()
  
  # Algorithm
  for (t in 1:T){
    p_t = trust_level_df[t,]
    
    w_ = (p_t*w)/sum(p_t*w)
    
    x_t = price_ratio_df[t,]
    
    X_t[t] <- sum(x_t*w_)
    
    
    # LOSS UPDATE
    w_m = (w * (p_t*x_t + (1 - p_t)*X_t[t]))/X_t[t]
    
    # MIXING(Fixed-Share) UPDATE
    #alpha = 1/t
    
    w = alpha/N + (1 - alpha)*w_m
    
  }
  
  K = cumprod(X_t) # portfolio wealth
  
  return(K)
}

run_main = function(price_ratio_df, herst_exp_df){
  # portfolio wealth vector for Buy and Hold algorithm
  K_n = rowMeans(sapply(price_ratio_df, cumprod))
  
  # algorithms evaluation output
  res = data.frame(a=0, b=0, alpha=0, profit=0, BandH=0, SR=0)
  
  for(l in 1:length(alpha)){
    for(i in 1:length(a)){
      for(j in 1:length(b)){
        # portfolio wealth vector for Portfolio Fixed-Share for unreliable instruments algorithm
        K = run_portfolio_fs(price_ratio_df, ht_to_pt(a[i],b[j], herst_exp_df), alpha[l])
        
        res[nrow(res) + 1,] = list(a[i],b[j],alpha[l],tail(K, n=1), tail(K_n, n=1),sum(K>K_n)/length(K))
      }
    }
    # portfolio wealth vector for Singer Portfolio algorithm
    K_z = run_portfolio_fs(price_ratio_df, (herst_exp_df*0)+1, alpha[l])
    
    res[nrow(res) + 1,] = list(1,1, alpha[l], tail(K_z, n=1), tail(K_n, n=1), sum(K_z>K_n)/length(K_z)) 
  }

  return (res)
  # par(mfrow = c(1, 2))
  # plot(K, pch = "*")
  # points(K_z, pch = "*", col = "blue")
  # points(K_n, pch = "*", col = "red")
}

#res11 = run_main(price_ratio_11, herst_exp_11)
res5 = run_main(price_ratio_5, herst_exp_5)
write.table(res5, file="res.txt")