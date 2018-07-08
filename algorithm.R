library(data.table)
library(reshape2)
library(ggplot2)
library(pracma)
library(useful)
library(rowr)
library(rlist)

input_path = file.path("input", "finam_raw", "stocks_10012012_10032013_AFLT_BANE.txt")
exp_len = 30

################################################# DATA LOADING #################################################
setClass('finamDate')
setAs("character","finamDate",
      function(from) {date_str <- as.character(from)
                      good_str <- ifelse(nchar(date_str)==6, date_str, paste("0", date_str, sep=""))
                      as.Date(good_str, "%d%m%y") })


Hs = function(ts){
  return(hurstexp(ts, display=FALSE)[1])
}


load_data = function(exp_len){
  finam_data = read.delim(input_path,
                    sep = ",",
                    col.names = c("ticker", "per", "date", "time",
                                  "open", "high", "low", "close", "vol"),
                    colClasses = c("factor", "NULL", "finamDate", "character",
                                   "numeric", "NULL", "NULL", "numeric", "NULL"))
  
  stocks_raw = rbindlist(lapply(split(finam_data, finam_data$ticker), 
                                function(df) {
                                  print (df[1,1])
                                  rbindlist(lapply(split(df, df$date),
                                                   function(x) {
                                                     data.frame(x[1, 1:2], 
                                                                price_ratio = tail(x$close, 1)/head(x$open, 1), 
                                                                ts = I(list(rbind(x$close, x$open))))
                                                   }))
                                }))

  # tail's argument n = -(window-1) to strip days without computed hurst exponent value
  stocks = rbindlist(lapply(split(stocks_raw, stocks_raw$ticker), 
                            function(df) {
                              print (df[1,1])
                              new_df = tail(subset(df, select = c(ticker, date, price_ratio)), -(exp_len-1))
                              new_df$hurst = as.numeric(rollApply(df$ts, function(x) Hs(c(unlist(x))), 
                                                        window=exp_len, minimum=exp_len, align='right'))
                              new_df
                            }))
  
  return(stocks)
}
#################################################################################################################


# hurst exponent transformation into trust levels
ht_to_pt_matlab = function(a, b, hurst){
  xi = function(a){ c(0, a-0.1, a, a+0.1, 1)}
  yi = function(b){ c(0, b, 0.5, 1-b, 1)}
  plot(pchipfun(xi(a),yi(b)))
  
  trust_levels = pchip(xi(a), yi(b), hurst)
  return(trust_levels)
}

ht_to_pt = function(a, b, hurst){
  xi = function(a){ c(0, 0.49, a, a+0.1, 1)}
  yi = function(b){ c(0, 0, 0, b, 1)}
  plot(pchipfun(xi(a),yi(b)))
  
  trust_levels = pchip(xi(a), yi(b), hurst)
  return(trust_levels)
}


# portfolio fixed share algorithm for unreliable instruments
run_portfolio_fs = function(stocks, alpha){
  # Initial parameters
  T = nlevels(stocks$date)      # number of steps
  N = nlevels(stocks$ticker)    # number of instruments
  
  w_ = c()                      # w^*
  w = rep(1/N, N)               # w
  w_m = rep(1/N, N)             # w^m
  W = data.frame(matrix(ncol = N, nrow = 0)) # w^* matrix
  colnames(W) = sort(levels(stocks$ticker))
  
  X_t = c()
  
  days = sort(levels(stocks$date))
  
  # Algorithm   ## !mind dates and tickers order!
  for (t in 1:T){
    inputs = subset(stocks, subset = date == days[t])
    
    p_t = inputs$trust_level
    
    # TRUST UPDATE
    # if trust levels for all instruments are zeros then do not make trust update
    w_ = if(sum(p_t)) (p_t*w)/sum(p_t*w) else w
    W[nrow(W) + 1,] = w_
    
    x_t = inputs$price_ratio
    X_t[t] = sum(x_t*w_)
    
    # LOSS UPDATE
    w_m = (w * (p_t*x_t + (1 - p_t)*X_t[t]))/X_t[t]
    
    # MIXING(Fixed-Share) UPDATE
    w = alpha(t)/N + (1 - alpha(t))*w_m
  }
  
  K = cumprod(X_t) # portfolio wealth
  
  return(K)
}


if (!exists("stocks")){
  stocks = data.frame(load_data(exp_len))
}
stocks$date = as.factor(stocks$date)

# filling missing dates with price_ratio = 1, hurst = 0
for (d in levels(stocks$date)){
    avail_tickers = subset(stocks, subset = date == d)$ticker
    lvls = levels(stocks$ticker)
    
    for (t in lvls[!lvls %in% avail_tickers]  ){
      stocks[nrow(stocks) + 1,] = list(t, d, 1, 0)
    }
}
# sorting for right dates order
stocks[order(stocks$ticker, stocks$date)]

############### write prepared stocks dataframe to file ###############
dump_filename = sprintf("%s_%sd_hurst.txt", 
                        gsub(".txt$", "", basename(input_path)), 
                        exp_len)
dump_path = file.path("input", dump_filename)
write.table(stocks, file=dump_path)
################################################################

# algorithm parameters
a = c(0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
b = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9)

const_alphas = c(0.001, 0.01, 0.1, 0.25, 1)
const_alpha_fun = function(x) { function(t) {x} }
alpha = list.append(sapply(const_alphas, FUN=const_alpha_fun), function(t) {1 / t})
alpha_label = c(const_alphas, "1/t")

# portfolio wealth vector for Buy and Hold algorithm
K_n = rowMeans(data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod)))

# constant rebalanced portfolio with 1/N
K_n_crp = cumprod(lapply(split(stocks$price_ratio, stocks$date), mean))

# best portfolio stock 
b_s = max(data.frame(lapply(split(stocks$price_ratio, stocks$ticker), prod)))

# algorithms evaluation output
res = data.frame(matrix(ncol = 11, nrow = 0))
names = c("a", "b", "alpha", "profit", "B&H profit", ">B&H", 
          "Singer profit", ">Singer", "CRP profit", ">CRP", "best stock")
colnames(res) = names

## TODO make rhis cycle -> function
for(l in 1:length(alpha)){
  print(paste("alpha", alpha_label[l]))
  
  # portfolio wealth vector for Singer Portfolio algorithm
  stocks$trust_level = rep(1,nrow(stocks))
  K_z = run_portfolio_fs(stocks, alpha[[l]])
  
  res[nrow(res) + 1,] = list(1, 1, alpha_label[l], tail(K_z, 1), 
                             tail(K_n, 1), sum(K_z>K_n)/length(K_z), 
                             tail(K_z, 1), 0, 
                             tail(K_n_crp, 1), sum(K_z>K_n_crp)/length(K_z), b_s) 
  
  for(i in 1:length(a)){
    for(j in 1:length(b)){

      # portfolio wealth vector for Portfolio Fixed-Share for unreliable instruments algorithm
      stocks$trust_level = ht_to_pt(a[i],b[j], stocks$hurst)
      #stocks$trust_level = stocks$hurst
      K = run_portfolio_fs(stocks, alpha[[l]])
      
      plot_data_raw = data.frame(x = as.numeric(1:length(K)), 
                                 "С уровнями доверия" = K, "Buy and Hold" = K_n,
                                 "CRP" = K_n_crp, "Зингера" = K_z)
      plot_data = melt(plot_data_raw, id="x")
      
      
      # TODO добавить параметры a,b,alpha на график and generate subtitle
      print(ggplot(data=plot_data,
                   aes(x=x, y=value, colour=variable)) +
              geom_point(size=0.4) +scale_colour_manual(values=c("orange", "blue", "darkgreen", "red")) +
              labs(x="Торговый день", y="Относительный капитал",
                   title = "Ежедневная доходность алгоритмов", subtitle = "AFLT, GMKN, YNDX", color='Портфель') +
              scale_x_continuous(breaks = seq(0, length(K), by = 150)) +
              theme(panel.background = element_rect(fill = '#ecf7ff'),
                    legend.key = element_rect(fill = "white"),
                    legend.position = c(.25, .95),
                    legend.justification = c("right", "top"),
                    legend.box.just = "right"))


      res[nrow(res) + 1,] = list(a[i], b[j], alpha_label[l], tail(K, 1),
                                 tail(K_n, 1), sum(K>K_n)/length(K),
                                 tail(K_z, 1), sum(K>K_z)/length(K),
                                 tail(K_n_crp, 1), sum(K>K_n_crp)/length(K), b_s)
    }
  }
}


# TODO find best algo parameters and best singer alpha and make table with x vectors for them

############################TODO make this part automatic########################################  
#K_bs = data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod))$GMKN
#sum(K>K_bs)/length(K)

# stocks_list = data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod))
# plot_data_raw1 = data.frame(x = as.numeric(1:length(K)),
#                             "Портфель" = K, "PIKK" = stocks_list$PIKK,
#                             "LKOH" = stocks_list$LKOH, "SIBN" = stocks_list$SIBN)
# plot_data1 = melt(plot_data_raw1, id="x")
# 
# ggplot(data=plot_data1, aes(x=x, y=value, colour=variable)) +
#   geom_line() +scale_colour_manual(values=c("red", "green", "lightblue", "yellow" )) +
#   labs(x="Торговый день", y="Относительный капитал",
#        title = "Ежедневная доходность акций составляющих портфель", subtitle = "PIKK, LKOH, SIBN", color='Инструмент') +
#   scale_x_continuous(breaks = seq(0, length(K), by = 150)) +
#   theme(panel.background = element_rect(fill = '#ecf7ff'),
#         legend.key = element_rect(fill = "white"),
#         legend.position = c(.25, .95),
#         legend.justification = c("right", "top"),
#         legend.box.just = "right")

# ggplot(data = diamonds) + 
#   geom_bar(mapping = aes(x = cut))


res_filename = sprintf("%s_%sd_hurst.txt", 
                       gsub(".txt$", "", gsub("^stocks_", "", basename(input_path))), 
                       exp_len)
output_path = file.path("results", res_filename)
write.table(res, file=output_path)
