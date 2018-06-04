price_ratio <- read.table("Ru5_data.csv", sep = ',', header = TRUE)
trust_level <- read.table("Ru5_hexp.csv", sep = ',', header = TRUE)

# Initial parameters
T <- nrow(price_ratio)      # number of steps

X_t <- c()

N <- ncol(price_ratio)    # number of instruments

#alpha = 1/t

w_ <- c()                      # w^*
w <- rep(1/N, N)               # w
w_m <- rep(1/N, N)             # w^m



# Algorithm
for (t in 1:T){
  p_t = trust_level[t,]
  
  w_ = (p_t*w)/sum(p_t*w)
  
  x_t = price_ratio[t,]
  
  X_t[t] <- sum(x_t*w_)
  
  
  # LOSS UPDATE
  w_m = (w * (p_t*x_t + (1 - p_t)*x_t))/X_t[t]
  
  # MIXING(Fixed-Share) UPDATE
  alpha = 1/t
  w = alpha/N + (1 - alpha)*w_m
  
}

K = cumprod(X_t) # portfolio wealth
