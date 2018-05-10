T <- 10      # number of steps

K_0 <- 1      # initial wealth
K <- c()      

N <- 9         # number of experts

alpha <- 0.5

w <- rep(1/N, N)
w_m <- rep(1/N, N)

companies <- levels(stocks$X.TICKER.)

# for (i in 1:nrow(stocks)){
#   date_str <- as.character(stocks$X.DATE.[i])
#   good_str <- ifelse(nchar(date_str)==6, date_str, paste("0", date_str, sep=""))
#   stocks$date[i] <- as.Date(good_str, "%d%m%y")
# }

#stocks$date <- factor(stocks$date)


for (t in 1:T){
  day <- stocks$date[order(levels(stocks$date))[t]]
  #print (day)
  x <- subset(stocks, subset = date == day)[c('X.TICKER.','X.OPEN.')]
  #print (x)
  X_t <- 0
  for (i in 1:N){
    X_t <- X_t + x[i,2]*w[i]
  }
  
  # LOSS UPDATE
  for (i in 1:N){
    w_m[i] <- (w[i] * x[i,2])/X_t
  }
  # MIXING(Fixed-Share) UPDATE
  for (i in 1:N){
    w[i] <- alpha/N + (1- alpha)*w_m[i]
  }
  print(paste("xt", X_t))
  print("wm")
  print(w_m)
  print("w")
  print(w)
}