library(data.table)
library(ggplot2)
library(pracma)
library(useful)
library(rowr)


################################### DATA LOADING #########################################
setClass('finamDate')
setAs("character","finamDate",
      function(from) {date_str <- as.character(from)
                      good_str <- ifelse(nchar(date_str)==6, date_str, paste("0", date_str, sep=""))
                      as.Date(good_str, "%d%m%y") })
                      #factor(as.Date(good_str, "%d%m%y")) })


Hs = function(ts){
  return(hurstexp(ts, display=FALSE)[1])
}


load_data = function(){
  finam_data = read.delim("stocks_10012012_10012018_short.txt",
                    sep = ",",
                    col.names = c("ticker", "per", "date", "time",
                                  "open", "high", "low", "close", "vol"),
                    colClasses = c("factor", "NULL", "finamDate", "character",
                                   "numeric", "NULL", "NULL", "numeric", "NULL"))
  
# обработать ошибку когда нет данных за какие то дни
# Error in data.frame(head(x, 1)[, 1:2], price_ratio = tail(x$close, 1)/head(x$open,  : 
#                   arguments imply differing number of rows: 0, 1
  # [1] ticker date   time   open   close 
  # <0 rows> (or 0-length row.names)
  #finam_data = finam_data[!(finam_data$ticker=="BANE" | finam_data$ticker=="GAZP"),]
  
  # stocks_raw = rbindlist(lapply(split(finam_data, list(finam_data$ticker, finam_data$date), lex.order = TRUE), 
  #                               function(x) {
  #                                 print (head(x, 1)[,1:2])
  #                                 data.frame(head(x,1)[,1:2], 
  #                                            price_ratio = tail(x$close, 1)/head(x$open, 1), 
  #                                            ts =I(list(rbind(x$close, x$open))))
  #                               }))
  
  stocks_raw = rbindlist(lapply(split(finam_data, finam_data$ticker, lex.order = TRUE), 
                            function(df) { 
                              rbindlist(lapply(split(df, df$date, lex.order = TRUE),
                              function(x) {
                                print (head(x, 1)[,1:2])
                                data.frame(head(x,1)[,1:2], 
                                           price_ratio = tail(x$close, 1)/head(x$open, 1), 
                                           ts =I(list(rbind(x$close, x$open))))
                              }))
                            }))
  
  # tail's argument n = -(window-1)
  stocks_raw = rbindlist(lapply(split(stocks_raw, stocks_raw$ticker), 
                            FUN=function(df) {
                              new_df = tail(subset(df, select = c(ticker, date, price_ratio)), -1)
                              new_df$hurst = as.numeric(rollApply(df$ts, function(x) Hs(c(unlist(x))), 
                                                        window=2,minimum=2, align='right'))
                              new_df
                              }))
  
  return(stocks_raw)
}
##########################################################################################


# hurst exponent transformation into trust levels
ht_to_pt = function(a, b, hurst){
  xi = function(a){ c(0, a-0.1, a, a+0.1, 1)}
  yi = function(b){ c(0, b, 0.5, 1-b, 1)}
  plot(pchipfun(xi(a),yi(b)))
  
  trust_levels = pchip(xi(a), yi(b), hurst)
  return(trust_levels)
}


# portfolio fixed share algorithm for unreliable instruments
run_portfolio_fs = function(stocks, alpha){
  # Initial parameters
  T <- nlevels(stocks$date)      # number of steps
  N <- nlevels(stocks$ticker)    # number of instruments
  
  w_ <- c()                      # w^*
  w <- rep(1/N, N)               # w
  w_m <- rep(1/N, N)             # w^m
  
  X_t <- c()
  
  day_nums = order(levels(stocks$date))
  
  # Algorithm   ###TODO:check for missing values and maybe order of instruments of p_t, x_t
  for (t in 1:T){
    day <- stocks$date[day_nums[t]]
    inputs = subset(stocks, subset = date == day)
    
    p_t = inputs$trust_level
    
    w_ = (p_t*w)/sum(p_t*w)
    
    x_t = inputs$price_ratio
    
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

#TODO: load only if there is no finam_data yet
stocks = load_data()
stocks$date = as.factor(stocks$date)

for (d in levels(stocks$date)){
  if (nrow(stocks[stocks$date==d,]) != 8){
    stocks = stocks[!(stocks$date==d),]
  }
}

# algorithm parameters

a = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
b = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.3, 0.4, 0.5)
alpha = c(0.001,0.01,0.1,0.25,1)

# portfolio wealth vector for Buy and Hold algorithm
K_n = rowMeans(data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod)))
K_n_crp = cumprod(lapply(split(stocks$price_ratio, stocks$date), mean))
                              
# algorithms evaluation output
res <- data.frame(matrix(ncol = 6, nrow = 0))
names <- c("a","b","alpha", "profit", "BandH", "SR")
colnames(res) <- names


for(l in 1:length(alpha)){
  for(i in 1:length(a)){
    for(j in 1:length(b)){
      stocks$trust_level = ht_to_pt(a[i],b[j], stocks$hurst)
      
      # portfolio wealth vector for Portfolio Fixed-Share for unreliable instruments algorithm
      K = run_portfolio_fs(stocks, alpha[l])
      
      res[nrow(res) + 1,] = list(a[i], b[j], alpha[l], tail(K, n=1), tail(K_n, n=1), sum(K>K_n)/length(K))
    }
  }
  stocks$trust_level = rep(1,nrow(stocks))
  
  # portfolio wealth vector for Singer Portfolio algorithm
  K_z = run_portfolio_fs(stocks, alpha[l])
  
  res[nrow(res) + 1,] = list(1, 1, alpha[l], tail(K_z, n=1), tail(K_n, n=1), sum(K_z>K_n)/length(K_z)) 
}

# par(mfrow = c(1, 2))
# plot(K, pch = "*")
# points(K_z, pch = "*", col = "blue")
# points(K_n, pch = "*", col = "red")

write.table(res, file="new_stocks_moex_top_short.txt")
