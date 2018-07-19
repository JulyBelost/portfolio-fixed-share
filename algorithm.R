pacman::p_load(reshape2, ggplot2, pracma, useful, rowr, rlist, data.table)

setClass('finamDate')
setAs("character","finamDate",
      function(from) {date_str <- as.character(from)
                      good_str <- ifelse(nchar(date_str)==6, date_str, paste("0", date_str, sep=""))
                      as.Date(good_str, "%d%m%y") })

# returns hurst exponent value instead of vector
hurstexp_ = function(ts){
  return(hurstexp(ts, display=FALSE)[1])
}


# loading data from file
load_data = function(input_path, exp_len){
  dump_filename = sprintf("%s_%sd_hurst.txt",
                          gsub(".txt$", "", basename(input_path)),
                          exp_len)
  dump_path = file.path(dirname(input_path), "..", dump_filename)
  
  if(file.exists(dump_path)){
    stocks = read.delim(dump_path, sep = " ")
    # sorting for right dates order
    stocks = stocks[order(stocks$ticker, stocks$date),]
    
    return(stocks)
  }
    
  finam_data = read.delim(input_path,
                    sep = ",",
                    col.names = c("ticker", "per", "date", "time",
                                  "open", "high", "low", "close", "vol"),
                    colClasses = c("factor", "NULL", "finamDate", "character",
                                   "numeric", "NULL", "NULL", "numeric", "NULL"))
  
  stocks_raw = rbindlist(lapply(split(finam_data, finam_data$ticker), 
                                function(df) {
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
                              new_df$hurst = as.numeric(rollApply(df$ts, function(x) hurstexp_(c(unlist(x))), 
                                                        window=exp_len, minimum=exp_len, align='right'))
                              new_df
                            }))
  
  # manipulation with stocks to get rid of "data.table" type
  stocks = data.frame(stocks)
  stocks$date = as.factor(stocks$date)
  
  # filling missing dates with price_ratio = 1, hurst = 0
  for(d in levels(stocks$date)){
    avail_tickers = subset(stocks, subset = date == d)$ticker
    lvls = levels(stocks$ticker)
    
    for(t in lvls[!lvls %in% avail_tickers]){
      stocks[nrow(stocks) + 1,] = list(t, d, 1, 0)
    }
  }
  
  # sorting for right dates order
  stocks = stocks[order(stocks$ticker, stocks$date),]
  
  # write prepared stocks dataframe to file
  write.table(stocks, file=dump_path)

  return(stocks)
}

# TODO improve function(add more points)
# hurst exponent transformation into trust levels
ht_to_pt = function(a, b, hurst){
  xi = function(a){ c(0, 0.49, a, a+0.1, 1)}
  yi = function(b){ c(0, 0, 0, b, 1)}
  #plot(pchipfun(xi(a),yi(b)))
  
  trust_levels = pchip(xi(a), yi(b), hurst)
  return(trust_levels)
}

ht_to_pt_matlab = function(a, b, hurst){
  xi = function(a){ c(0, a-0.1, a, a+0.1, 1)}
  yi = function(b){ c(0, b, 0.5, 1-b, 1)}
  plot(pchipfun(xi(a),yi(b)))
  
  trust_levels = pchip(xi(a), yi(b), hurst)
  return(trust_levels)
}


# portfolio fixed share algorithm for unreliable instruments
run_portfolio_fs = function(stocks, alpha, verbose = FALSE){
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
    # if trust levels for all instruments are zeros then do not make trust update TODO!!! check
    w_ = if(sum(p_t)) (p_t*w)/sum(p_t*w) else w
    W[nrow(W) + 1,] = w_
    
    x_t = inputs$price_ratio
    X_t[t] = sum(x_t*w_)
    
    # LOSS UPDATE
    w_m = (w * (p_t*x_t + (1 - p_t)*X_t[t]))/X_t[t]
    
    # MIXING(Fixed-Share) UPDATE
    w = alpha(t)/N + (1 - alpha(t))*w_m
  }
  
  # TODO add p levels as well
  if (verbose){
    plot_W = melt(data.frame(x = as.numeric(1:T), W), id="x")
    print(ggplot(data=plot_W, aes(x=x, y=value, fill=variable)) + geom_area() +
      scale_fill_brewer(palette="Blues"))
  }
  
  K = cumprod(X_t) # portfolio wealth
  return(K)
}


process_portfolio = function(input_path, exp_len, dump_only = FALSE){
  stocks = load_data(input_path, exp_len)
  
  # exit function if only prepared data dump needed
  if(dump_only) { return() }
  
  # portfolio wealth vector for Buy and Hold algorithm
  K_n = rowMeans(data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod)))
  # constant rebalanced portfolio with 1/N
  K_n_crp = cumprod(lapply(split(stocks$price_ratio, stocks$date), mean))
  # best portfolio stock portfolio
  K_stocks = data.frame(lapply(split(stocks$price_ratio, stocks$ticker), cumprod))
  best_stock = names(which.max(tail(K_stocks, 1)))[1]
  K_bs = K_stocks[,best_stock]
  
  # algorithms evaluation output
  res = data.frame(matrix(ncol = 12, nrow = 0))
  names = c("a", "b", "alpha", "profit", "B&H profit", ">B&H", 
            "Singer profit", ">Singer", "CRP profit", ">CRP", "best stock", ">bs")
  colnames(res) = names
  
  best_K = c(-1)
  best_K_z = c(-1)
  
  for(l in 1:length(alpha)){
    print(paste("alpha", alpha_label[l]))
    
    # portfolio wealth vector for Singer Portfolio algorithm
    stocks$trust_level = rep(1,nrow(stocks))
    K_z = run_portfolio_fs(stocks, alpha[[l]])
    
    res[nrow(res) + 1,] = list(1, 1, alpha_label[l], tail(K_z, 1),
                               tail(K_n, 1), sum(K_z>K_n)/length(K_z),
                               tail(K_z, 1), 0,
                               tail(K_n_crp, 1), sum(K_z>K_n_crp)/length(K_z), 
                               tail(K_bs, 1), sum(K_z>K_bs)/length(K_z)) 
    
    if(tail(K_z, 1) > tail(best_K_z, 1)){
      best_K_z = K_z
      bzalpha = alpha[l]
    }
    
    for(i in 1:length(a)){
      for(j in 1:length(b)){
  
        # portfolio wealth vector for Portfolio Fixed-Share for unreliable instruments algorithm
        # consider to try stocks$trust_level = stocks$hurst
        stocks$trust_level = ht_to_pt(a[i],b[j], stocks$hurst)
        K = run_portfolio_fs(stocks, alpha[[l]])
  
        res[nrow(res) + 1,] = list(a[i], b[j], alpha_label[l], tail(K, 1),
                                   tail(K_n, 1), sum(K>K_n)/length(K),
                                   tail(K_z, 1), sum(K>K_z)/length(K),
                                   tail(K_n_crp, 1), sum(K>K_n_crp)/length(K), 
                                   tail(K_bs, 1), sum(K>K_bs)/length(K))
        
        if(tail(K, 1) > tail(best_K, 1)){
          best_K = K
          ba = a[i]
          bb = b[j]
          balpha = alpha[l]
        }
      }
    }
  }
  
  best_res = res[order(res$profit, decreasing=TRUE),][1:5,]
  best_res[,'name'] = paste(colnames(K_stocks), collapse = ", ")
  
  # plot chart with all portfolio performance
  pplot_data_raw = data.frame(x = as.numeric(1:length(K)), 
                             "With trust levels" = K, "Buy and Hold" = K_n, "CRP" = K_n_crp, "Singer" = K_z)
  pplot_data = melt(pplot_data_raw, id="x")
  portf_plot = ggplot(data=pplot_data, aes(x=x, y=value, colour=variable)) +
          geom_point(size=0.4) + scale_colour_manual(values=c("orange", "blue", "darkgreen", "red")) +
          labs(x="Trading day", y="Wealth",
               subtitle = sprintf("%s (a=%s, b=%s, alpha=%s)", 
                                  paste(colnames(K_stocks), collapse = ", "), a[i], b[j], alpha_label[l]),
               color='Portfolio') +
          scale_x_continuous(breaks = seq(0, length(K), by = 150)) +
          theme(panel.background = element_rect(fill = '#ecf7ff'),
                legend.key = element_rect(fill = "white"),
                legend.position = c(.3, .95),
                legend.justification = c("right", "top"),
                legend.box.just = "right")
  print(portf_plot)
  
  # plot chart with individual stocks performance
  splot_data_raw = data.frame(x = as.numeric(1:length(K)), "Portfolio" = K, K_stocks)
  splot_data = melt(splot_data_raw, id="x")
  stocks_plot = ggplot(data=splot_data, aes(x=x, y=value, colour=variable)) +
    geom_line() +
    labs(x="Trading day", y="Wealth", color='Stock') +
    scale_x_continuous(breaks = seq(0, length(K), by = 150)) +
    theme(panel.background = element_rect(fill = '#ecf7ff'),
          legend.key = element_rect(fill = "white"),
          legend.position = c(.25, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right")
  print(stocks_plot)
  
  # result files saving
  output_path = file.path(dirname(input_path), "..", "..", "results")
  res_filename = sprintf("%s_%sd_hurst.txt", 
                         gsub(".txt$", "", gsub("^stocks_", "", basename(input_path))), 
                         exp_len)
  res_path = file.path(output_path, res_filename)
  print(res_path)
  write.table(res, file=res_path)
  
  write.table(pplot_data, file=file.path(output_path, sprintf("pplot_data_%s", res_filename)))
  write.table(splot_data, file=file.path(output_path, sprintf("splot_data_%s", res_filename)))
  
  pplot_filename = sprintf("%s_portfolios.pdf", gsub(".txt$", "", res_filename))
  splot_filename = sprintf("%s_stocks.pdf", gsub(".txt$", "", res_filename))
  pdf(file.path(output_path, pplot_filename))
  portf_plot
  dev.off()
  pdf(file.path(output_path, splot_filename))
  stocks_plot
  dev.off()
  
  return(best_res)
}


############################ algorithm hyperparameters ############################
a = c(0.5, 0.6, 0.7, 0.8, 0.9)
b = c(0.9,0.8,0.7,0.6,0.5, 0.3, 0.2, 0.1, 0.05, 0.01)

const_alphas = c(0.0001) #, 0.001, 0.01, 0.1, 0.25, 1)
const_alpha_fun = function(x) { function(t) {x} }
alpha = list.append(sapply(const_alphas, FUN=const_alpha_fun), function(t) {1 / t})
alpha_label = c(const_alphas, "1/t")
###################################################################################

# if called from command line use $ Rscript algorithm.R exp_len portf_fold
# where exp_len and portf_fold are single values
# otherwise if executed in Rstudio uses exp_len and port_folders defined below in else section
args = commandArgs(trailingOnly = TRUE)
if(len(args) == 2){
  exp_len = c(strtoi(args[1]))
  port_folders = c(args[2])
} else {
  exp_len = c(10, 20, 30)
  port_folders = c("portf_size0") #, "portf_size3", "portf_size4", "portf_size5", "portf_size6")
}

# iterate through files in choosen folders and calls process_portfolio function for them with different exp_len values
process_folders = function(port_folders, exp_len, dump_only = TRUE){
  for (fold in port_folders){
    input_dir = file.path("exp0_market200", fold, "input", "finam_raw")
    files = list.files(path=input_dir, pattern="*.txt", full.names = TRUE)
    
    for (e in exp_len){
      for (file in files){
        print(file)
        print(paste("exp_len =", e))
          best_params_df = process_portfolio(file, e, dump_only)
      }
    }
  }
}

process_folders(port_folders, exp_len, dump_only = FALSE)

# 
# TODO solve case when one stock is way better
# TODO add instrument with zero profit as a way of take out all the money
# TODO check plots format for journal
# TODO make table with x vectors for best K
