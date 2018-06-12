#########################################################################################
# setClass('finamDate')
# setAs("character","finamDate", 
#       function(from) {date_str <- as.character(from)
#                       good_str <- ifelse(nchar(date_str)==6, date_str, paste("0", date_str, sep=""))
#                       factor(as.Date(good_str, "%d%m%y")) })
# 
# stocks = read.delim("finam_output.txt", 
#                   sep = ",", 
#                   col.names = c("ticker", "per", "date", "time",  
#                                 "open", "high", "low", "close", "vol"),
#                   colClasses = c("factor", "NULL", "finamDate", "NULL",
#                                  "NULL", "NULL", "NULL", "numeric", "NULL"))
# 
# stocks = transform(stocks, 
#             price_ratio=ave(close, ticker, FUN=function(y) c(1, tail(y, -1) / head(y, -1))))
# 
# #добавить матрицу с р-шками, где они считаются. добавить столбец с ними в stocks
# stocks$trust_level <- rep(1)
################################### DATA LOADING ########################################


companies <- levels(stocks$ticker)

# Initial parameters
T <- nlevels(stocks$date)      # number of steps

X_t <- c()

N <- nlevels(stocks$ticker)    # number of instruments

#alpha = 1/t

w_ <- c()                      # w^*
w <- rep(1/N, N)               # w
w_m <- rep(1/N, N)             # w^m



# Algorithm
for (t in 1:T){
  day <- stocks$date[order(levels(stocks$date))[t]]

  inputs = subset(stocks, subset = date == day)
  p_t = inputs$trust_level
  
  w_ = (p_t*w)/sum(p_t*w)
  
  x_t = inputs$price_ratio
  
  X_t[t] <- sum(x_t*w_)

  
  # LOSS UPDATE
  w_m = (w * (p_t*x_t + (1 - p_t)*X_t[t]))/X_t[t]

  # MIXING(Fixed-Share) UPDATE
  alpha = 1/t
  w = alpha/N + (1 - alpha)*w_m
  
}

K = cumprod(X_t) # portfolio wealth
