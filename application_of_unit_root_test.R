rm(list = ls())
setwd("C:/Users/Neel/Desktop/2023-2024/Spring 2024/Global Economics/R Functions")
pwd = getwd() # getting the path for present working directory (pwd)

##########################################################
# Adding required packages and functions
##########################################################

source("model_selection_function.R") # function for model selection
source("ols_function.R") # function for OLS estimation
source("t_test_function.R") # fucntion for t test
source("expanding_window_forecast_function.R") # function for forecasting using expanding windows
source("rolling_window_forecast_function.R") # function for forecasting using rolling windows
library("quantmod") # add quantmod to the list of Packages
library("fBasics") # add fBasics to the list of Packages
library("fUnitRoots") # add fUnitRoots to the list of Packages

##########################################################
# Fetch and filter the data 
##########################################################

getSymbols(Symbols ="GDPC1",src = "FRED") # Download Quarterly data for real GDP
Y <- as.matrix(GDPC1[,1])
DY <- diff(Y)/Y[1:(dim(Y)[1]-1),]
head(Y)
tail(Y)
keep_data <- seq(from = as.Date("1947-04-01"), to = as.Date("2023-10-1"), by = "quarter")
Y_new = as.matrix(Y[as.Date(rownames(Y)) %in% keep_data,]) # precovid data
DY_new = as.matrix(DY[as.Date(rownames(DY)) %in% keep_data,]) # precovid data

n_obs = dim(DY_new)[1]
Y_new_date = as.Date(row.names(Y_new))
DY_new_date = Y_new_date
plot(Y_new_date,Y_new,ylab='Real GDP',type='l')

##########################################################
# Augmented Dickey Fuller Test
##########################################################
pacf(DY_new,lag=round(n_obs^(1/3))) # command to obtain sample PACF of the data

results <- model_selection(round(n_obs^(1/3)),DY_new)
aic_values = results$AIC
bic_values = results$BIC
num_lags_aic = results$op_lag_AIC  
num_lags_bic = results$op_lag_BIC

num_lags_aic
num_lags_bic

adfTest(Y_new,lags=num_lags_aic,type=c("ct")) # ADF test using AIC selected number of lags
adfTest(Y_new,lags=num_lags_bic,type=c("ct")) # ADF test using BIC selected number of lags

##########################################################
# Plot ACF and perform auto correlation test
##########################################################
acf(Y_new,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data

acf(DY_new,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data
Box.test(DY_new, lag = round(n_obs^(1/3)), type = "Ljung-Box") # applying Ljung and Box (1978) joint test of auto correlations

##########################################################
# Theoritical ACF and model checking using AIC 
##########################################################
num_lags = num_lags_aic
lags_DY_new = matrix(NA,nrow = n_obs, ncol = num_lags)
for (i in 1:num_lags) {
  lags_DY_new[(i+1):n_obs,i] = as.matrix(DY_new[1:(n_obs-i),1])
}
intercept = matrix(1,n_obs)
X = cbind(intercept,lags_DY_new)
y = DY_new
reg_result = ols(X[(num_lags+1):n_obs,],as.matrix(y[(num_lags+1):n_obs,1]))
beta_hat = reg_result$beta_hat
var_beta_hat = reg_result$var_beta_hat
test_result = t_test(beta_hat,var_beta_hat)
test_result$t_stat
test_result$p_value

ar_coeff <- as.numeric(beta_hat[2:(num_lags+1)])
ma_coeff <- 0
ACF = acf(DY_new,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data
TACF <- ARMAacf(ar_coeff, ma_coeff, lag.max = round(n_obs^(1/3))) # command to obtain theorical ACF
plot(c(0:round(n_obs^(1/3))),ACF$acf,type='l',xlab='Lag',ylab='ACF',ylim=c(-1,1))
lines(0:round(n_obs^(1/3)),TACF,lty=2)
grid(nx = 4, ny = 4)

residuals = reg_result$u_hat # get the AR model residuals
acf(residuals,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data
Box.test(residuals, lag = round(n_obs^(1/3)), type = "Ljung-Box") # applying Ljung and Box (1978) joint test of auto correlations

##########################################################
# Theoritical ACF and model checking using BIC 
##########################################################
num_lags = num_lags_bic
lags_DY_new = matrix(NA,nrow = n_obs, ncol = num_lags)
for (i in 1:num_lags) {
  lags_DY_new[(i+1):n_obs,i] = as.matrix(DY_new[1:(n_obs-i),1])
}
intercept = matrix(1,n_obs)
X = cbind(intercept,lags_DY_new)
y = DY_new
reg_result = ols(X[(num_lags+1):n_obs,],as.matrix(y[(num_lags+1):n_obs,1]))
beta_hat = reg_result$beta_hat
var_beta_hat = reg_result$var_beta_hat
test_result = t_test(beta_hat,var_beta_hat)
test_result$t_stat
test_result$p_value

ar_coeff <- as.numeric(beta_hat[2:(num_lags+1)])
ma_coeff <- 0
ACF = acf(DY_new,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data
TACF <- ARMAacf(ar_coeff, ma_coeff, lag.max = round(n_obs^(1/3))) # command to obtain theorical ACF
plot(c(0:round(n_obs^(1/3))),ACF$acf,type='l',xlab='Lag',ylab='ACF',ylim=c(-0.2,1))
lines(0:round(n_obs^(1/3)),TACF,lty=2)
grid(nx = 4, ny = 4)

residuals = reg_result$u_hat # get the AR model residuals
acf(residuals,lag=round(n_obs^(1/3))) # command to obtain sample ACF of the data
Box.test(residuals, lag = round(n_obs^(1/3)), type = "Ljung-Box") # applying Ljung and Box (1978) joint test of auto correlations

#########################################################################
# Forecasting using AR Model with aic selected number of lags
#########################################################################
lag_choice = NA
init_win_len = 120 # the first 30 years
num_step_ahead = 8 # 1 to 8 steps ahead forecastes 
prediction_results = expanding_window(y = DY_new, init_win_len = init_win_len, pre_sel_num_lags = lag_choice, num_step_ahead = num_step_ahead, sel_method = 'aic')
y_f_aic <- prediction_results$actual_value
yhat_f_aic <- prediction_results$forecast
selected_num_lags <- prediction_results$sel_num_lags

plot(x = DY_new_date[121:n_obs], y = y_f_aic,xlab='time',ylab='GDP growth',type='l',col="yellow")
lines(x = DY_new_date[121:n_obs],y = yhat_f_aic[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_aic[,8],lty=3, col = 2)

forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_aic) - yhat_f_aic
rmsfe_ar_aic = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))

#########################################################################
# Forecasting using AR Model with bic selected number of lags
#########################################################################
lag_choice = NA
init_win_len = 120 # the first 30 years
num_step_ahead = 8 # 1 to 8 steps ahead forecastes 
prediction_results = expanding_window(y = DY_new, init_win_len = init_win_len, pre_sel_num_lags = lag_choice, num_step_ahead = num_step_ahead, sel_method = 'bic')
y_f_bic <- prediction_results$actual_value
yhat_f_bic <- prediction_results$forecast
selected_num_lags <- prediction_results$sel_num_lags

plot(x = DY_new_date[121:n_obs], y = y_f_bic,xlab='time',ylab='GDP growth',type='l',col="yellow")
lines(x = DY_new_date[121:n_obs],y = yhat_f_bic[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_bic[,8],lty=3, col = 2)

forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_bic) - yhat_f_bic
rmsfe_ar_bic = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))

#########################################################################
# Forecasting by Model averaging and RMSFE comparison
#########################################################################

yhat_f_ave = (yhat_f_aic + yhat_f_bic)/2
forecast_error =  kronecker(matrix(1,ncol = num_step_ahead),y_f_bic) - yhat_f_ave
rmsfe_ave = sqrt(colMeans(forecast_error^2, na.rm = TRUE, dims = 1))

plot(x = DY_new_date[121:n_obs], y = y_f_bic,xlab='time',ylab='Change in number of new infected cases',type='l',col="yellow")
lines(x = DY_new_date[121:n_obs],y = yhat_f_ave[,1],lty=2, col = 4)
lines(x = DY_new_date[121:n_obs],y = yhat_f_ave[,7],lty=3, col = 2)

rmsfe_all = rbind(rmsfe_ar_aic,rmsfe_ar_bic,rmsfe_ave)
rmsfe_all
