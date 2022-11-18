## ------  case with N(0,1) errors 
library(stlplus)
set.seed(12)
ndays <- 365; nyears <- 10; start_time <- c(2012,1)
n <- ndays *nyears
t <- 1 : (ndays*nyears)
w <- (2*pi)/ndays 
sea <- 5*sin(w*t) #seasonal 
tr <- 0.1*(t/ndays - 1)* (t/ndays- 6)* (t/ndays - 10)
#tr <- tr/sd(tr)
#err = rnorm(length(t))

#err = 3*rnorm(length(t))



err <- arima.sim(n =ndays*nyears, list(ar = c(0.5), ma = c(1, 0.5)))

# adding outliers
idx <- sample(n, round(0.1*n))
Xt <- sea + tr + err
#Xt[idx] <- rnorm(length(idx), 30, 2)

Xt <- ts(Xt, start=start_time , frequency = ndays) 

#tmp <- rbinom(length(idx), 1,0.5)
#Xt[idx] <- tmp*rnorm(length(idx), 25, 2) + (1-tmp)*rnorm(length(idx), -25, 2)


## stl_plus
res_stl_plus <- stlplus::stlplus(Xt, n.p=ndays, s.window="periodic", inner=2, outer=10)

## stl does not work "series is not periodic or has less than two periods"; check frequency first ... 
res_stl <- stl(Xt,  s.window="per", inner=2, outer=10)

## check results: 
plot(Xt, type = "p", col="grey")
lines(ts(res_stl_plus$data[,"seasonal"], start =start_time , frequency = ndays), lwd = 2, type = "l", col = "blue")
lines(ts(res_stl_plus$data[,"trend"], start = start_time , frequency = ndays), lwd = 2, type = "l", col = "blue") 

lines(ts(res_stl$time.series[,"seasonal"], start =start_time, frequency = ndays), lwd = 2, type = "l", col="green")
lines(ts(res_stl$time.series[,"trend"], start =start_time, frequency = ndays), lwd = 2, type = "l", col="green")

lines(ts(sea, start =start_time, frequency = ndays), lwd = 2, type = "l", col="red")
lines(ts(tr, start =start_time, frequency = ndays), lwd = 2, type = "l", col="red")


## try to mimic what stl_plus is doing ... 
cycle_indices <- rep(rep(1:ndays), nyears)
cs1 <- cs2 <- rep(1:ndays) # head and tail 
n.p <- ndays
c_tilde <- rep(NA, n)
for(i in 1:ndays){
  idx <- which(cycle_indices== i)
  c_tilde[idx] <- mean(Xt[which(cycle_indices==i)])
}
plot(c_tilde)
ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
c_tilde_ma <- ma(c_tilde, n.p) # this is a constant 
L <- c_tilde_ma
seasonal <- c_tilde - L   # not smooth; just the average of days per year   
D <- as.numeric(Xt) - seasonal # this is trend 

# smooth the trend to ignore things < 20 per days .... 
# Let's try to code Jason's method 

init.trend <- function(x, ndays, nyears){
  
  x_sub <- matrix(x, ncol = ndays)
  seasonal_med <- apply(x_sub, 2, median)
  #trend_rough <- ma(x-rep(seasonal_med, nyears), ndays)
  tmp <- as.numeric(x-rep(seasonal_med, nyears))
  time_pt <- 1:length(x)  
  trend_init <- loess(tmp~time_pt, family = "symmetric") # might be a little oversmoothed... so seasonal may have short term, but should be 
  #okay since the short term differ in years, and likely not be picked up by the next seasonal smoothing 
  return(predict(trend_init, time_pt))
}


cv.seasonal <- function(x_sub, candidate_bw){  
  
  time_pt <- 1:ncol(x_sub)  
  cv_err <- 0
  for(i in 1:nrow(x_sub)){
    tmp <- as.numeric(x_sub[-i,])
    loess_obj <- loess(tmp~rep(time_pt, nrow(x_sub)-1), span =  candidate_bw,  family = "symmetric")
    
    cv_err <- cv_err + sum( (x_sub[i,] - predict(loess_obj, time_pt))^2)
  }
  cv_err  <- cv_err/nrow(x_sub) # average prediction error per year	
  return(cv_err)
}

fit.seasonal <- function(x_sub, bw){
  
  time_pt <- 1:ncol(x_sub)  
  tmp <- as.numeric(x_sub)
  loess_obj <- loess(tmp~rep(time_pt, nrow(x_sub)), span =  bw,  family = "symmetric") 
  fitted <- predict(loess_obj, time_pt)

}

# if bw too large, trend flat; else trend rough 
decompose_trend_remainder <- function(x_trend_remainder, bw){
  
  time_pt <- 1:length(x_trend_remainder)  
  tmp <- as.numeric(x_trend_remainder)
  loess_obj <- loess(tmp~time_pt, span =  bw,  family = "symmetric") 
  fitted_trend <- predict(loess_obj, time_pt)
  return(list(fitted_trend = fitted_trend, fitted_remainder = x_trend_remainder - fitted_trend))
  
}


perform.decompose  <- function(x, ndays, nyears, niter){
  
  x_detrend <- x - init.trend(x, ndays, nyears) # first detrend
  
  for(i in 1:niter){
    
    x_sub <- matrix(x_detrend, ncol = ndays)
    opt_bw <- optim(par = 0.5, fn = function(candidate_bw){cv.seasonal(x_sub, candidate_bw)}, lower = 0.01, upper = 1, method = "Brent") 
    x_seasonal <- fit.seasonal(x_sub, opt_bw$par)
    x_seasonal_rep <- rep(x_seasonal, nyears)
    x_trend_remainder <- x - x_seasonal_rep
    
    par(mfrow = c(3,2))
    
    par(mfrow = c(1,1))
    
    for(bw in c(0.1, 0.5, 2)){
      
      tmp <- decompose_trend_remainder(x_trend_remainder, bw)
      plot(tmp$fitted_trend, main = var(tmp$fitted_trend), ylab = "fitted trend")
      plot(tmp$fitted_remainder, main = var(tmp$fitted_remainder), ylab = "fitted residual")
      
    }
    
    
    #library("funtimes")
    time_points <- 1:(ndays* nyears) 
    tmpp <- loess(tmp$fitted_remainder~time_points,  family = "symmetric")   # not super smooth but should be good enough .. 
    tmpp <- predict(tmpp, time_points)
    plot(tmp$fitted_trend)
    plot(tmpp)
    
    plot(abs(tmp$fitted_trend))
    plot(abs(tmpp))
    
    
    aa <- quantile(abs(tmp$fitted_trend), seq(0,1,0.1))
    bb <- quantile(abs(tmpp), seq(0,1,0.1))
    plot(bb~aa)
    abline(a = 0, b = 1, col = "red")
    
    
    t.test(abs(tmpp))
    
    tmpp <- tmpp[1:200]
    time_points1 <- 1:length(tmpp)
    
    wavk_test(tmpp ~ time_points1, factor.length = "adaptive.selection")
    
    tmpp  <- 0.5*time_points1 + rnorm(length(tmpp))
    wavk_test(tmpp ~ time_points1 + time_points1^2, factor.length = "adaptive.selection")
    
    
    
    var_bw <- NULL; var_remainder <- NULL
    mse_tr <-var_tr <- bias_tr <- NULL
    for(bw in seq(0.1,5,0.1)){
      tmp <- decompose_trend_remainder(x_trend_remainder, bw)
      var_bw <- c(var_bw, var(tmp$fitted_trend))
      mse_tr <- c(mse_tr, mean((tmp$fitted_trend -x_trend_remainder)^2))
      var_tr <- c(var_tr, var(tmp$fitted_trend -x_trend_remainder))
      bias_tr <- c(bias_tr, mean(tmp$fitted_trend -x_trend_remainder)^2)
      var_remainder <- c(var_remainder, mean( (tmp$fitted_remainder)^2))
    }
    
    
    #plot(seq(0.1,5,0.1),mse_tr)
    #plot(seq(0.1,5,0.1),bias_tr)
    
    
    par(mfrow = c(1,2))
    plot(seq(0.1,5,0.1),var_bw)
    plot(seq(0.1,5,0.1),var_remainder)
    
    #plot(var_bw+ var_remainder)
  }
  
}    


set.seed(25468)
n_var <- 2
fname <- ZDT3
lower <- rep(0, n_var)
upper <- rep(1, n_var)
res <- easyGParetoptim(fn=fname, lower=lower, upper=upper, budget=15,
                       control=list(method="EHI", inneroptim="pso", maxit=20))
par(mfrow=c(1,2))
plotGPareto(res)
title("Pareto Front")
plot(res$history$X, main="Pareto set", col = "red", pch = 20)
points(res$par, col="blue", pch = 17)







