library(data.table)
library(reshape2)
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


load_data = function(exp_len){
  finam_data = read.delim("stocks_10012012_10012018_short_new.txt",
                    sep = ",",
                    col.names = c("ticker", "per", "date", "time",
                                  "open", "high", "low", "close", "vol"),
                    colClasses = c("factor", "NULL", "finamDate", "character",
                                   "numeric", "NULL", "NULL", "numeric", "NULL"))
  

  
  stocks_raw = rbindlist(lapply(split(finam_data, finam_data$ticker, lex.order = TRUE), 
                            function(df) {
                              print (df[1,1])
                              rbindlist(lapply(split(df, df$date, lex.order = TRUE),
                              function(x) {
                                data.frame(head(x,1)[,1:2], 
                                           price_ratio = tail(x$close, 1)/head(x$open, 1), 
                                           ts =I(list(rbind(x$close, x$open))))
                              }))
                            }))
  
  # tail's argument n = -(window-1)
  stocks_raw = rbindlist(lapply(split(stocks_raw, stocks_raw$ticker), 
                            FUN=function(df) {
                              new_df = tail(subset(df, select = c(ticker, date, price_ratio)), -(exp_len-1))
                              new_df$hurst = as.numeric(rollApply(df$ts, function(x) Hs(c(unlist(x))), 
                                                        window=exp_len,minimum=exp_len, align='right'))
                              new_df
                              }))
  
  return(stocks_raw)
}
##########################################################################################


# hurst exponent transformation into trust levels
ht_to_pt = function(a, b, hurst){
  xi = function(a){ c(0, 0.49, a, a+0.1, 1)}
  yi = function(b){ c(0, 0, 0, 1-b, 1)}
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
    
    w_ = if(sum(p_t)) (p_t*w)/sum(p_t*w) else w
    
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
stocks = data.frame(load_data(20))
stocks$date = as.factor(stocks$date)

for (d in levels(stocks$date)){
    avail_tickers = subset(stocks, subset = date == d)$ticker
    lvls = levels(stocks$ticker)
    
    for (t in lvls[!lvls %in% avail_tickers]  ){
      stocks[nrow(stocks) + 1,] = list(t, d, 1, 0)
    }
}

# delete dates with even one missing instrument
#for (d in levels(stocks$date)){
  #if (nrow(stocks[stocks$date==d,]) != nlevels(stocks$ticker)){
  #stocks = stocks[!(stocks$date==d),]
  # }
#}

stocks$date = factor(stocks$date)


# algorithm parameters

a = c(0.5, 0.6, 0.7, 0.8, 0.9)
b = c(0.5, 0.7, 0.8, 0.9, 0.95, 0.99)
alpha = c(0.001, 0.01, 0.1, 0.25, 1)

# portfolio wealth vector for Buy and Hold algorithm
K_n = rowMeans(data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod)))

# constant rebalanced portfolio with 1/N
K_n_crp = cumprod(lapply(split(stocks$price_ratio, stocks$date), mean))

                              
# algorithms evaluation output
res <- data.frame(matrix(ncol = 11, nrow = 0))
names <- c("a", "b", "alpha", "profit", "B&H profit", ">B&H", 
           "Singer profit", ">Singer", "CRP profit", ">CRP", "best stock")
colnames(res) <- names

b_s = max(data.frame(lapply(split(stocks$price_ratio, stocks$ticker), prod)))

for(l in 1:length(alpha)){
  print(paste("alpha", alpha[l]))
  
  # portfolio wealth vector for Singer Portfolio algorithm
  stocks$trust_level = rep(1,nrow(stocks))
  K_z = run_portfolio_fs(stocks, alpha[l])
  
  res[nrow(res) + 1,] = list(1, 1, alpha[l], tail(K_z, 1), 
                             tail(K_n, 1), sum(K_z>K_n)/length(K_z), 
                             tail(K_z, 1), 0, 
                             tail(K_n_crp, 1), sum(K_z>K_n_crp)/length(K_z), b_s) 
  
  for(i in 1:length(a)){
    for(j in 1:length(b)){

      # portfolio wealth vector for Portfolio Fixed-Share for unreliable instruments algorithm
      stocks$trust_level = ht_to_pt(a[i],b[j], stocks$hurst)
      #stocks$trust_level = stocks$hurst
      K = run_portfolio_fs(stocks, alpha[l])
      
      plot_data_raw = data.frame(x = as.numeric(1:length(K)), 
                                 "с уровнями доверия" = K, "автоматический(1/N)" = K_n, 
                                 "постоянный(1/N)" = K_n_crp, "Зингера" = K_z)
      plot_data = melt(plot_data_raw, id="x")
      
      ggplot(data=plot_data,
             aes(x=x, y=value, colour=variable)) +
        geom_point(size=0.4) +scale_colour_manual(values=c("orange", "blue", "darkgreen", "red")) + 
        labs(x="Торговый день", y="Относительный капитал", 
             title = "Ежедневный доход", subtitle = "AFLT, GMKN, YNDX", color='Портфель') + 
        scale_x_continuous(breaks = seq(0, length(K), by = 150)) +
        theme(panel.background = element_rect(fill = '#ecf7ff'),
              legend.key = element_rect(fill = "white"),
              legend.position = c(.25, .95),
              legend.justification = c("right", "top"),
              legend.box.just = "right") 

      
      res[nrow(res) + 1,] = list(a[i], b[j], alpha[l], tail(K, 1),
                                 tail(K_n, 1), sum(K>K_n)/length(K),
                                 tail(K_z, 1), sum(K>K_z)/length(K),
                                 tail(K_n_crp, 1), sum(K>K_n_crp)/length(K), b_s)
    }
  }
}

  #par(mfrow = c(1, 2))

#write.table(res, file="AFLT_GMKN_YAND_2014_2018_10days.txt")
