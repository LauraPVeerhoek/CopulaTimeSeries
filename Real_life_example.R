library(astsa)
library(copula)
library(svines)
library(forecast)
library(vars)
library(tsDyn)
library(functClust)
library(tseries)


plot(prodn, ylab = "Federal Reserve Board Production Index", main = "Monthly Federal Reserve Board Production Index Time Series")

dprodn = diff(prodn)
plot(dprodn,  ylab = "Once Differenced time series", main = "Once Differenced Time Series")
acf2(dprodn)

ddprodn = diff(dprodn,12)
plot(ddprodn,  ylab = "Once differenced and seasonality removed", main = "Once Differenced and Seasonality Removed Time Series")

Box.test(ddprodn, type='Ljung-Box') # small p-value suggest stationarity
adf.test(ddprodn, alternative = 'stationary') # large p-value suggests stationarity
kpss.test(ddprodn) # large p-value suggests stationarity

train_data = ddprodn[1:287] 
test_data = ddprodn[288:359]


AICS_train = c()

for (i in 1:100){
  print(i)
  AICS_train = c(AICS_train, AIC(svine(train_data, p = i)))
}
plot(AICS_train)


fit_ar_train = auto.arima(train_data, max.q = 0, max.d = 0, max.P = 0, max.Q = 0, max.D = 0)
summary(fit_ar_train)
mean(abs(residuals(fit_ar_train)))

prediction_ar = forecast(fit_ar_train, h=72, level = c(95))

mean(abs((test_data-prediction_ar$mean)))

length(which((test_data<prediction_ar$upper)&(test_data>prediction_ar$lower)))/72

# determining the vine order 
AICS_train = c()

for (i in 1:100){
  print(i)
  AICS_train = c(AICS_train, AIC(svine(train_data, p = i)))
}
plot(AICS_train, xlab = 'order of the vine model', ylab = "AIC of the vine model", main = "AIC's of vine models with orders 1-100")

#fit vine model of order 2

vine_fit_2 = svine(train_data, p=2)

summary(vine_fit_2)

model_pred_med_2 = c()

for (i in 14:287){ #warning: takes long to run 
  print(i)
  vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = train_data[c((i-13):(i-1))])
  model_pred_med_2 = c(model_pred_med_2, quantile(vine_prediction,  probs = c(0.5)))
}

mean(abs(train_data[14:287]- model_pred_med_2))

logLik(vine_fit_2)
AIC(vine_fit_2)

vine_prediction_2 <- svine_sim(n = 72, rep = 1000, model = vine_fit_2, past = train_data)
vine_prediction_med_2 = c()
vine_prediction_lwr_2 = c()
vine_prediction_upr_2 = c()

for (j in 1:72){
  vine_prediction_med_2 = c(vine_prediction_med_2, quantile(vine_prediction_2[j, ,],  probs = c(0.5)))
  vine_prediction_lwr_2 = c(vine_prediction_lwr_2, quantile(vine_prediction_2[j, ,],  probs = c(0.025)))
  vine_prediction_upr_2 = c(vine_prediction_upr_2, quantile(vine_prediction_2[j, ,],  probs = c(0.975)))
}

mean(abs((test_data-vine_prediction_med_2)))

length(which((test_data<vine_prediction_upr_2)&(test_data>vine_prediction_lwr_2)))/72

vine_fit = svine(train_data, p=13)

summary(vine_fit)

model_pred_med = c()

for (i in 14:287){ #warning: takes long to run 
  print(i)
  vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = train_data[c(i-13, i-12, i-11, i-10, i-9, i-8, i-7, i-6, i-5, i-4, i-3, i-2, i-1),])
  model_pred_med = c(model_pred_med, quantile(vine_prediction,  probs = c(0.5)))
}

mean(abs(train_data[14:287]- model_pred_med))

logLik(vine_fit)
AIC(vine_fit)

vine_prediction <- svine_sim(n = 72, rep = 1000, model = vine_fit, past = train_data)
vine_prediction_med = c()
vine_prediction_lwr = c()
vine_prediction_upr = c()

for (j in 1:72){
  vine_prediction_med = c(vine_prediction_med, quantile(vine_prediction[j, ,],  probs = c(0.5)))
  vine_prediction_lwr = c(vine_prediction_lwr, quantile(vine_prediction[j, ,],  probs = c(0.025)))
  vine_prediction_upr = c(vine_prediction_upr, quantile(vine_prediction[j, ,],  probs = c(0.975)))
}

mean(abs((test_data-vine_prediction_med)))

length(which((test_data<vine_prediction_upr)&(test_data>vine_prediction_lwr)))/72


# plot of the data
plot(test_data, type='l', ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds")
lines(c(1:72), prediction_ar$mean, col='blue')
lines(c(1:72), prediction_ar$lower, col='blue')
lines(c(1:72), prediction_ar$upper, col='blue')
lines(c(1:72), vine_prediction_lwr, col='red')
lines(c(1:72), vine_prediction_med, col='red')
lines(c(1:72), vine_prediction_upr, col='red')
lines(c(1:72), vine_prediction_lwr_2, col='green')
lines(c(1:72), vine_prediction_med_2, col='green')
lines(c(1:72), vine_prediction_upr_2, col='green')
legend("bottomright", 
       legend = c("Data", "AR model(2)", "Vine model(2)", "Vine model(13)"), 
       col =  c("black", "blue", "green", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )



#plot of the sizes of the prediction interval
sizes_ar = c()
sizes_vine_2 =c()
sizes_vine_13 =c()

for(i in 1:72){
  sizes_ar = c(sizes_ar, abs(prediction_ar$upper[i] - prediction_ar$lower[i]))
  sizes_vine_2 = c(sizes_vine_2, abs(vine_prediction_upr_2[i] - vine_prediction_lwr_2[i]))
  sizes_vine_13 = c(sizes_vine_13, abs(vine_prediction_upr[i] - vine_prediction_lwr[i]))
}


plot(sizes_ar, type='l', ylim = c(min(sizes_vine_13), max(sizes_vine_2)), ylab = 'size of the prediction bounds at time t in the future', 
     xlab = "time t", main = 'Sizes of the prediction bounds', col='blue')
lines(sizes_vine_13, col='red')
lines(sizes_vine_2, col='green')

legend("bottomright", 
       legend = c( "AR model(2)", "Vine model(2)", "Vine model(13)"), 
       col =  c( "blue", "green", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )


######################################################################################################################
#MULTIVARIATE CASE

y = econ5[,c(1,5)]
plot(y, main = "Multivariate time series")
dy = diff(y)
plot(dy, main = "Once differenced multivariate time series")

Box.test(dy[,1], type='Ljung-Box')
Box.test(dy[,2], type='Ljung-Box')


adf.test(dy[,1], alternative = 'stationary')
adf.test(dy[,2], alternative = 'stationary')

# does not pass adf test 
kpss.test(dy[,1])
kpss.test(dy[,2])


train_data = dy[1:128,]
test_data = dy[128:160,]


var_model = VAR(train_data,p=VARselect(dy)$selection[1]) # if you want to get the coefficients nicely use ar(data_train)
AIC(var_model)
summary(var_model)
logLik(var_model)
(mean(abs(residuals(var_model)[,1]))+ mean(abs(residuals(var_model)[,2])))/2
prediction_ar = predict(var_model, n.ahead = 33, ci=.95)
( mean(abs(c(prediction_ar$fcst$unemp[,1] - test_data[,1]))) + mean(abs(c(prediction_ar$fcst$prinv[,1] - test_data[,2]))))/2
length(which((test_data[,1]<prediction_ar$fcst$unemp[,3])&
               (test_data[,1]>prediction_ar$fcst$unemp[,2])&
               (test_data[,2]<prediction_ar$fcst$prinv[,3])&
               (test_data[,2]>prediction_ar$fcst$prinv[,2])))/33



AICS_train = c()

for (i in 1:20){
  print(i)
  AICS_train = c(AICS_train, AIC(svine(train_data, p = i)))
}
plot(AICS_train)




vine_fit = svine(train_data, p=1)
AIC(vine_fit)
logLik(vine_fit)


summary(vine_fit)

# model_pred_med = array(c(NA,NA),dim = c(1,2))
# 
# for (i in 2:1000){ #warning: takes long to run 
#   vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = matrix(train_data[i-1,], ncol = 2))
#   model_pred_med = rbind(model_pred_med,c(quantile(vine_prediction[,1,],  probs = c(0.5)),quantile(vine_prediction[,2,],  probs = c(0.5))))
# }
# 
# model_pred_med = model_pred_med[-1,]
# ( mean(abs(train_data[2:1000,1]- model_pred_med[,1])) + mean(abs(train_data[2:1000,2]- model_pred_med[,2])))/2


vine_prediction <- svine_sim(n = 33, rep = 1000, model = vine_fit, past = train_data)
vine_prediction_med_1 = c()
vine_prediction_lwr_1 = c()
vine_prediction_upr_1 = c()
vine_prediction_med_2 = c()
vine_prediction_lwr_2 = c()
vine_prediction_upr_2 = c()

for (j in 1:33){
  vine_prediction_med_1 = c(vine_prediction_med_1, quantile(vine_prediction[j,1 ,],  probs = c(0.5)))
  vine_prediction_lwr_1 = c(vine_prediction_lwr_1, quantile(vine_prediction[j,1 ,],  probs = c(0.025)))
  vine_prediction_upr_1 = c(vine_prediction_upr_1, quantile(vine_prediction[j,1 ,],  probs = c(0.975)))
  vine_prediction_med_2 = c(vine_prediction_med_2, quantile(vine_prediction[j,2 ,],  probs = c(0.5)))
  vine_prediction_lwr_2 = c(vine_prediction_lwr_2, quantile(vine_prediction[j,2 ,],  probs = c(0.025)))
  vine_prediction_upr_2 = c(vine_prediction_upr_2, quantile(vine_prediction[j,2 ,],  probs = c(0.975)))
}

(mean(abs(c(vine_prediction_med_1 - test_data[,1]))) + mean(abs(c(vine_prediction_med_2 - test_data[,2]))))/2

length(which((test_data[,1]<vine_prediction_upr_1)&
               (test_data[,1]>vine_prediction_lwr_1)&
               (test_data[,2]<vine_prediction_upr_2)&
               (test_data[,2]>vine_prediction_lwr_2)))/33



