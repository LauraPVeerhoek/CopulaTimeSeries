library(copula)
library(svines)
library(forecast)
library(vars)
library(tsDyn)
library(functClust)

#################################     model 1    ##############################################

order_ar = c()
MAE_ar = c()
loglik_ar = c()
AIC_ar = c()
MAE_test_ar = c()
inside_conf_int_test_ar = c()
margin_vine = c()
MAE_vine = c()
loglik_vine = c()
AIC_vine = c()
MAE_test_vine = c()
inside_conf_int_test_vine = c()

for (i in 1:50){
  print(i)
  #simulate the data from ar(1) with parameter 0.5
  x <- arima.sim(n = 1100, model = list(order(1, 0, 0), ar = 0.5))
  data_train = x[1:1000]
  data_test = x[1001:1100]
  
  #fit AR model on the data
  fit_ar = auto.arima(data_train, max.q = 0, max.d = 0, max.P = 0, max.Q = 0, max.D = 0)

  order_ar = c(order_ar, fit_ar$arma[1])
  
  #MAE generally in the order of 0.1 
  MAE_ar = c(MAE_ar, mean(abs(residuals(fit_ar))))
  
  
  loglik_ar = c(loglik_ar, fit_ar$loglik)
  AIC_ar = c(AIC_ar, fit_ar$aic)


  prediction_ar = forecast(fit_ar, h=100, level = c(95))

  MAE_test_ar = c(MAE_test_ar, mean(abs((data_test-prediction_ar$mean))))
  
  inside_conf_int_test_ar = c(inside_conf_int_test_ar, length(which((data_test<prediction_ar$upper)&(data_test>prediction_ar$lower))))
  
  vine_fit = svine(data_train, p=1)
  
  margin_vine = c(margin_vine, lapply(summary(vine_fit)$margins[1,3][1], as.character)[[1]])
  
  model_pred_med = c()

  for (i in 2:1000){ #warning: takes long to run 
    vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = data_train[i-1])
    model_pred_med = c(model_pred_med, quantile(vine_prediction,  probs = c(0.5)))
  }
  
  MAE_vine = c(MAE_vine, mean(abs(data_train[2:1000]- model_pred_med)))
  
  
  loglik_vine = c(loglik_vine, logLik(vine_fit))
  AIC_vine = c(AIC_vine, AIC(vine_fit))

  vine_prediction <- svine_sim(n = 100, rep = 1000, model = vine_fit, past = data_train)
  vine_prediction_med = c()
  vine_prediction_lwr = c()
  vine_prediction_upr = c()
  
  for (j in 1:100){
    vine_prediction_med = c(vine_prediction_med, quantile(vine_prediction[j, ,],  probs = c(0.5)))
    vine_prediction_lwr = c(vine_prediction_lwr, quantile(vine_prediction[j, ,],  probs = c(0.025)))
    vine_prediction_upr = c(vine_prediction_upr, quantile(vine_prediction[j, ,],  probs = c(0.975)))
  }

  MAE_test_vine = c(MAE_test_vine, mean(abs((data_test-vine_prediction_med))))
  
  inside_conf_int_test_vine = c(inside_conf_int_test_vine, length(which((data_test<vine_prediction_upr)&(data_test>vine_prediction_lwr))))
}
data.frame( MAE = c(mean(MAE_ar), mean(MAE_vine)),
            MAE_sd = c(sd(MAE_ar), sd(MAE_vine)),
            LogLik = c(mean(loglik_ar), mean(loglik_vine)),
            LogLik_sd = c(sd(loglik_ar), sd(loglik_vine)),
            AIC = c(mean(AIC_ar), mean(AIC_vine)),
            AIC_sd = c(sd(AIC_ar), sd(AIC_vine)),
            MAE_test = c(mean(MAE_test_ar), mean(MAE_test_vine)),
            MAE_test_sd = c(sd(MAE_test_ar), sd(MAE_test_vine)),
            inside_conf_int_test = c(mean(inside_conf_int_test_ar), mean(inside_conf_int_test_vine)),
            inside_conf_int_test_sd = c(sd(inside_conf_int_test_ar), sd(inside_conf_int_test_vine)),
            row.names = c("AR model", "vine model"))


order_ar
margin_vine

# plot of the data
plot(x[1000:1100], type='l', ylim = c(min(x[1000:1100]),max(x[1000:1100])), ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds")
lines(c(1:100), prediction_ar$mean, col='blue')
lines(c(1:100), prediction_ar$lower, col='blue')
lines(c(1:100), prediction_ar$upper, col='blue')
lines(c(1:100), vine_prediction_lwr, col='red')
lines(c(1:100), vine_prediction_med, col='red')
lines(c(1:100), vine_prediction_upr, col='red')
legend("bottomright", 
       legend = c("Data", "AR model", "Vine model"), 
       col =  c("black", "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )


results_1 = data.frame(order_ar,
                       MAE_ar,
                       loglik_ar,
                       AIC_ar,
                       MAE_test_ar,
                       inside_conf_int_test_ar,
                       margin_vine,
                       MAE_vine,
                       loglik_vine,
                       AIC_vine,
                       MAE_test_vine,
                       inside_conf_int_test_vine)

#plot of the sizes of the prediction interval
sizes_ar = c()
sizes_vine =c()
for(i in 1:100){
  sizes_ar = c(sizes_ar, abs(prediction_ar$upper[i] - prediction_ar$lower[i]))
  sizes_vine = c(sizes_vine, abs(vine_prediction_upr[i] - vine_prediction_lwr[i]))
}


plot(sizes_ar, type='l', ylim = c(min(sizes_vine), max(sizes_vine)), ylab = 'size of the prediction bounds at time t in the future', 
     xlab = "time t", main = 'Sizes of the prediction bounds', col='blue')
lines(sizes_vine, col='red')
legend("bottomright", 
       legend = c("AR model", "Vine model"), 
       col =  c( "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )


length(which(results_1$AIC_ar>results_1$AIC_vine)) # 16%
length(which(results_1$margin_vine=="Normal")) # 62%
mean(results_1$order_ar) # 1.32
length(which(results_1$order_ar>1)) # 26%



#################################     model 2    ##############################################
order_ar = c()
MAE_ar = c()
loglik_ar = c()
AIC_ar = c()
MAE_test_ar = c()
inside_conf_int_test_ar = c()
margin_vine = c()
MAE_vine = c()
loglik_vine = c()
AIC_vine = c()
MAE_test_vine = c()
inside_conf_int_test_vine = c()

for (i in 1:50){
  print(i)
  #simulate the data from ar(1) with parameter 0.5
  x <- arima.sim(n = 1100, model = list(order(1, 0, 0), ar = 0.5))
  x = rank(x)/1101 # same as pobs(x)
  data_train = x[1:1000]
  data_test = x[1001:1100]
  
  #fit AR model on the data
  fit_ar = auto.arima(data_train, max.q = 0, max.d = 0, max.P = 0, max.Q = 0, max.D = 0)
  
  order_ar = c(order_ar, fit_ar$arma[1])
  
  #MAE generally in the order of 0.1 
  MAE_ar = c(MAE_ar, mean(abs(residuals(fit_ar))))
  
  
  loglik_ar = c(loglik_ar, fit_ar$loglik)
  AIC_ar = c(AIC_ar, fit_ar$aic)

  
  prediction_ar = forecast(fit_ar, h=100, level = c(95))
  
  MAE_test_ar = c(MAE_test_ar, mean(abs((data_test-prediction_ar$mean))))
  
  inside_conf_int_test_ar = c(inside_conf_int_test_ar, length(which((data_test<prediction_ar$upper)&(data_test>prediction_ar$lower))))
  

  

  
  vine_fit = svine(data_train, p=1)
  
  margin_vine = c(margin_vine, lapply(summary(vine_fit)$margins[1,3][1], as.character)[[1]])
  
  model_pred_med = c()
  
  for (i in 2:1000){ #warning: takes long to run 
    vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = data_train[i-1])
    model_pred_med = c(model_pred_med, quantile(vine_prediction,  probs = c(0.5)))
  }
  
  MAE_vine = c(MAE_vine, mean(abs(data_train[2:1000]- model_pred_med)))
  
  
  loglik_vine = c(loglik_vine, logLik(vine_fit))
  AIC_vine = c(AIC_vine, AIC(vine_fit))

  vine_prediction <- svine_sim(n = 100, rep = 1000, model = vine_fit, past = data_train)
  vine_prediction_med = c()
  vine_prediction_lwr = c()
  vine_prediction_upr = c()
  
  for (j in 1:100){
    vine_prediction_med = c(vine_prediction_med, quantile(vine_prediction[j, ,],  probs = c(0.5)))
    vine_prediction_lwr = c(vine_prediction_lwr, quantile(vine_prediction[j, ,],  probs = c(0.025)))
    vine_prediction_upr = c(vine_prediction_upr, quantile(vine_prediction[j, ,],  probs = c(0.975)))
  }
  
  MAE_test_vine = c(MAE_test_vine, mean(abs((data_test-vine_prediction_med))))
  
  inside_conf_int_test_vine = c(inside_conf_int_test_vine, length(which((data_test<vine_prediction_upr)&(data_test>vine_prediction_lwr))))
  
  
  

  
}
data.frame( MAE = c(mean(MAE_ar), mean(MAE_vine)),
            MAE_sd = c(sd(MAE_ar), sd(MAE_vine)),
            LogLik = c(mean(loglik_ar), mean(loglik_vine)),
            LogLik_sd = c(sd(loglik_ar), sd(loglik_vine)),
            AIC = c(mean(AIC_ar), mean(AIC_vine)),
            AIC_sd = c(sd(AIC_ar), sd(AIC_vine)),
            MAE_test = c(mean(MAE_test_ar), mean(MAE_test_vine)),
            MAE_test_sd = c(sd(MAE_test_ar), sd(MAE_test_vine)),
            inside_conf_int_test = c(mean(inside_conf_int_test_ar), mean(inside_conf_int_test_vine)),
            inside_conf_int_test_sd = c(sd(inside_conf_int_test_ar), sd(inside_conf_int_test_vine)),
            row.names = c("AR model", "vine model"))


order_ar
margin_vine


plot(x[1000:1100], type='l', ylim = c(min(prediction_ar$lower),max(prediction_ar$upper)), ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds")
lines(c(1:100), prediction_ar$mean, col='blue')
lines(c(1:100), prediction_ar$lower, col='blue')
lines(c(1:100), prediction_ar$upper, col='blue')
lines(c(1:100), vine_prediction_lwr, col='red')
lines(c(1:100), vine_prediction_med, col='red')
lines(c(1:100), vine_prediction_upr, col='red')
legend("bottomright", 
       legend = c("Data", "AR model", "Vine model"), 
       col =  c("black", "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )

#This does not work very well, probably due to the margin estimation. So, let us rank the data, and rerun the models. 

results_2 = data.frame(order_ar,
                       MAE_ar,
                       loglik_ar,
                       AIC_ar,
                       MAE_test_ar,
                       inside_conf_int_test_ar,
                       margin_vine,
                       MAE_vine,
                       loglik_vine,
                       AIC_vine,
                       MAE_test_vine,
                       inside_conf_int_test_vine)


length(which(results_2$AIC_ar>results_2$AIC_vine)) # 100%
length(which(results_2$margin_vine=="Uniform")) # 96%
mean(results_2$order_ar) # 1.3 
length(which(results_2$order_ar>1)) # 28%

#################################     model 3    ##############################################

order_ar = c()
MAE_ar = c()
loglik_ar = c()
AIC_ar = c()
MAE_test_ar = c()
inside_conf_int_test_ar = c()
margin_vine = c()
MAE_vine = c()
loglik_vine = c()
AIC_vine = c()
MAE_test_vine = c()
inside_conf_int_test_vine = c()

for (i in 1:50){
  print(i)
  #simulate the data from copula
  # create marginal model
  univ1 <- univariateML::mlnorm(rnorm(10))
  univ1[] <- c(5, 10) # define parameter
  
  cs_struct <- cvine_structure(1)
  
  pcs <- list(
    list(  # first tree
      bicop_dist("clayton", 0, 3)# the copula
    )
  )
  cop <- svinecop_dist(
    pcs, cs_struct, p = 1, out_vertices = 1, in_vertices = 1)
  
  model <- svine_dist(margins = list(univ1), copula = cop)
  x <- svine_sim(n = 1100, rep = 1, model = model)
  
  data_train = x[1:1000]
  data_test = x[1001:1100]
  
  #fit AR model on the data
  fit_ar = auto.arima(data_train, max.q = 0, max.d = 0, max.P = 0, max.Q = 0, max.D = 0)
  
  order_ar = c(order_ar, fit_ar$arma[1])
  
  #MAE generally in the order of 0.1 
  MAE_ar = c(MAE_ar, mean(abs(residuals(fit_ar))))
  
  
  loglik_ar = c(loglik_ar, fit_ar$loglik)
  AIC_ar = c(AIC_ar, fit_ar$aic)

  
  prediction_ar = forecast(fit_ar, h=100, level = c(95))
  
  MAE_test_ar = c(MAE_test_ar, mean(abs((data_test-prediction_ar$mean))))
  
  inside_conf_int_test_ar = c(inside_conf_int_test_ar, length(which((data_test<prediction_ar$upper)&(data_test>prediction_ar$lower))))
  
  vine_fit = svine(data_train, p=1)
  
  margin_vine = c(margin_vine, lapply(summary(vine_fit)$margins[1,3][1], as.character)[[1]])
  
  model_pred_med = c()
  
  for (i in 2:1000){ #warning: takes long to run 
    vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = data_train[i-1])
    model_pred_med = c(model_pred_med, quantile(vine_prediction,  probs = c(0.5)))
  }
  
  MAE_vine = c(MAE_vine, mean(abs(data_train[2:1000]- model_pred_med)))
  
  
  loglik_vine = c(loglik_vine, logLik(vine_fit))
  AIC_vine = c(AIC_vine, AIC(vine_fit))

  vine_prediction <- svine_sim(n = 100, rep = 1000, model = vine_fit, past = data_train)
  vine_prediction_med = c()
  vine_prediction_lwr = c()
  vine_prediction_upr = c()
  
  for (j in 1:100){
    vine_prediction_med = c(vine_prediction_med, quantile(vine_prediction[j, ,],  probs = c(0.5)))
    vine_prediction_lwr = c(vine_prediction_lwr, quantile(vine_prediction[j, ,],  probs = c(0.025)))
    vine_prediction_upr = c(vine_prediction_upr, quantile(vine_prediction[j, ,],  probs = c(0.975)))
  }
  
  MAE_test_vine = c(MAE_test_vine, mean(abs((data_test-vine_prediction_med))))
  
  inside_conf_int_test_vine = c(inside_conf_int_test_vine, length(which((data_test<vine_prediction_upr)&(data_test>vine_prediction_lwr))))
}
data.frame( MAE = c(mean(MAE_ar), mean(MAE_vine)),
            MAE_sd = c(sd(MAE_ar), sd(MAE_vine)),
            LogLik = c(mean(loglik_ar), mean(loglik_vine)),
            LogLik_sd = c(sd(loglik_ar), sd(loglik_vine)),
            AIC = c(mean(AIC_ar), mean(AIC_vine)),
            AIC_sd = c(sd(AIC_ar), sd(AIC_vine)),
            MAE_test = c(mean(MAE_test_ar), mean(MAE_test_vine)),
            MAE_test_sd = c(sd(MAE_test_ar), sd(MAE_test_vine)),
            inside_conf_int_test = c(mean(inside_conf_int_test_ar), mean(inside_conf_int_test_vine)),
            inside_conf_int_test_sd = c(sd(inside_conf_int_test_ar), sd(inside_conf_int_test_vine)),
            row.names = c("AR model", "vine model"))


order_ar
margin_vine

plot(x[1000:1100], type='l', ylim = c(min(prediction_ar$lower),max(x[1000:1100])), ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds")
lines(c(1:100), prediction_ar$mean, col='blue')
lines(c(1:100), prediction_ar$lower, col='blue')
lines(c(1:100), prediction_ar$upper, col='blue')
lines(c(1:100), vine_prediction_lwr, col='red')
lines(c(1:100), vine_prediction_med, col='red')
lines(c(1:100), vine_prediction_upr, col='red')
legend("bottomright", 
       legend = c("Data", "AR model", "Vine model"), 
       col =  c("black", "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )

#This does not work very well, probably due to the margin estimation. So, let us rank the data, and rerun the models. 

results_3 = data.frame(order_ar,
                       MAE_ar,
                       loglik_ar,
                       AIC_ar,
                       MAE_test_ar,
                       inside_conf_int_test_ar,
                       margin_vine,
                       MAE_vine,
                       loglik_vine,
                       AIC_vine,
                       MAE_test_vine,
                       inside_conf_int_test_vine)

length(which(results_3$AIC_ar>results_3$AIC_vine)) # 100%
length(which(results_3$margin_vine=="Normal")) # 12% normal rest: as before 
mean(results_3$order_ar) #2.7
sd(results_3$order_ar) # 1.28
length(which(results_3$order_ar>1)) #80%


#################################     model 4    ##############################################

order_ar = c()
MAE_ar = c()
loglik_ar = c()
AIC_ar = c()
MAE_test_ar = c()
inside_conf_int_test_ar = c()
margin_vine = c()
MAE_vine = c()
loglik_vine = c()
AIC_vine = c()
MAE_test_vine = c()
inside_conf_int_test_vine = c()

for (i in 1:50){
  print(i)
  #simulate the data from copula
  # create marginal model
  univ1 <- univariateML::mlunif(runif(10))
  univ1[] <- c(0, 1) # define parameter
  
  cs_struct <- cvine_structure(1)
  
  pcs <- list(
    list(  # first tree
      bicop_dist("clayton", 0, 3)# the copula
    )
  )
  cop <- svinecop_dist(
    pcs, cs_struct, p = 1, out_vertices = 1, in_vertices = 1)
  
  model <- svine_dist(margins = list(univ1), copula = cop)
  x <- svine_sim(n = 1100, rep = 1, model = model)
  
  data_train = x[1:1000]
  data_test = x[1001:1100]
  
  #fit AR model on the data
  fit_ar = auto.arima(data_train, max.q = 0, max.d = 0, max.P = 0, max.Q = 0, max.D = 0)
  
  order_ar = c(order_ar, fit_ar$arma[1])
  
  #MAE generally in the order of 0.1 
  MAE_ar = c(MAE_ar, mean(abs(residuals(fit_ar))))
  
  
  loglik_ar = c(loglik_ar, fit_ar$loglik)
  AIC_ar = c(AIC_ar, fit_ar$aic)
  
  
  prediction_ar = forecast(fit_ar, h=100, level = c(95))
  
  MAE_test_ar = c(MAE_test_ar, mean(abs((data_test-prediction_ar$mean))))
  
  inside_conf_int_test_ar = c(inside_conf_int_test_ar, length(which((data_test<prediction_ar$upper)&(data_test>prediction_ar$lower))))
  
  vine_fit = svine(data_train, p=1)
  
  margin_vine = c(margin_vine, lapply(summary(vine_fit)$margins[1,3][1], as.character)[[1]])
  
  model_pred_med = c()
  
  for (i in 2:1000){ #warning: takes long to run 
    vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = data_train[i-1])
    model_pred_med = c(model_pred_med, quantile(vine_prediction,  probs = c(0.5)))
  }
  
  MAE_vine = c(MAE_vine, mean(abs(data_train[2:1000]- model_pred_med)))
  
  
  loglik_vine = c(loglik_vine, logLik(vine_fit))
  AIC_vine = c(AIC_vine, AIC(vine_fit))
  
  vine_prediction <- svine_sim(n = 100, rep = 1000, model = vine_fit, past = data_train)
  vine_prediction_med = c()
  vine_prediction_lwr = c()
  vine_prediction_upr = c()
  
  for (j in 1:100){
    vine_prediction_med = c(vine_prediction_med, quantile(vine_prediction[j, ,],  probs = c(0.5)))
    vine_prediction_lwr = c(vine_prediction_lwr, quantile(vine_prediction[j, ,],  probs = c(0.025)))
    vine_prediction_upr = c(vine_prediction_upr, quantile(vine_prediction[j, ,],  probs = c(0.975)))
  }
  
  MAE_test_vine = c(MAE_test_vine, mean(abs((data_test-vine_prediction_med))))
  
  inside_conf_int_test_vine = c(inside_conf_int_test_vine, length(which((data_test<vine_prediction_upr)&(data_test>vine_prediction_lwr))))
}
data.frame( MAE = c(mean(MAE_ar), mean(MAE_vine)),
            MAE_sd = c(sd(MAE_ar), sd(MAE_vine)),
            LogLik = c(mean(loglik_ar), mean(loglik_vine)),
            LogLik_sd = c(sd(loglik_ar), sd(loglik_vine)),
            AIC = c(mean(AIC_ar), mean(AIC_vine)),
            AIC_sd = c(sd(AIC_ar), sd(AIC_vine)),
            MAE_test = c(mean(MAE_test_ar), mean(MAE_test_vine)),
            MAE_test_sd = c(sd(MAE_test_ar), sd(MAE_test_vine)),
            inside_conf_int_test = c(mean(inside_conf_int_test_ar), mean(inside_conf_int_test_vine)),
            inside_conf_int_test_sd = c(sd(inside_conf_int_test_ar), sd(inside_conf_int_test_vine)),
            row.names = c("AR model", "vine model"))


order_ar
margin_vine


plot(x[1000:1100], type='l', ylim = c(min(prediction_ar$lower),max(prediction_ar$upper)), ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds")
lines(c(1:100), prediction_ar$mean, col='blue')
lines(c(1:100), prediction_ar$lower, col='blue')
lines(c(1:100), prediction_ar$upper, col='blue')
lines(c(1:100), vine_prediction_lwr, col='red')
lines(c(1:100), vine_prediction_med, col='red')
lines(c(1:100), vine_prediction_upr, col='red')
legend("bottomright", 
       legend = c("Data", "AR model", "Vine model"), 
       col =  c("black", "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )

results_4 = data.frame(order_ar,
                       MAE_ar,
                       loglik_ar,
                       AIC_ar,
                       MAE_test_ar,
                       inside_conf_int_test_ar,
                       margin_vine,
                       MAE_vine,
                       loglik_vine,
                       AIC_vine,
                       MAE_test_vine,
                       inside_conf_int_test_vine)


data.frame( MAE = c(mean(results_4$MAE_ar), mean(results_4$MAE_vine)),
            MAE_sd = c(sd(results_4$MAE_ar), sd(results_4$MAE_vine)),
            LogLik = c(mean(results_4$loglik_ar), mean(results_4$loglik_vine)),
            LogLik_sd = c(sd(results_4$loglik_ar), sd(results_4$loglik_vine)),
            AIC = c(mean(results_4$AIC_ar), mean(results_4$AIC_vine)),
            AIC_sd = c(sd(results_4$AIC_ar), sd(results_4$AIC_vine)),
            MAE_test = c(mean(results_4$MAE_test_ar), mean(results_4$MAE_test_vine)),
            MAE_test_sd = c(sd(results_4$MAE_test_ar), sd(results_4$MAE_test_vine)),
            inside_conf_int_test = c(mean(results_4$inside_conf_int_test_ar), mean(results_4$inside_conf_int_test_vine)),
            inside_conf_int_test_sd = c(sd(results_4$inside_conf_int_test_ar), sd(results_4$inside_conf_int_test_vine)),
            row.names = c("AR model", "vine model"))


results_4$order_ar
results_4$margin_vine

length(which(results_4$AIC_ar>results_4$AIC_vine)) # 100%
length(which(results_4$margin_vine=="Uniform")) # 48% uniform rest: beta, logit-normal, power distribution, kumaraswamy

mean(results_4$order_ar) #2.7
sd(results_4$order_ar) # 1.28
length(which(results_4$order_ar>1)) #80%



#################################     model 5    ##############################################
order_ar = c()
MAE_ar_1 = c()
MAE_ar_2 = c()
loglik_ar = c()
AIC_ar = c() 
MAE_test_ar_1 = c()
MAE_test_ar_2 = c()
inside_conf_int_test_ar_1 = c()
inside_conf_int_test_ar_2 = c()
inside_conf_int_test_ar_both = c()
margin_vine_1 = c()
margin_vine_2 = c()
MAE_vine_1 = c()
MAE_vine_2 = c()
loglik_vine = c()
AIC_vine = c()
MAE_test_vine_1 = c()
MAE_test_vine_2 = c()
inside_conf_int_test_vine_both = c()
inside_conf_int_test_vine_1 = c()
inside_conf_int_test_vine_2 = c()


for (i in 1:50){
  print(i)
  #simulate the data from VAR process
  A <- matrix(c(0.5,0.3, 0.02, 0.8), byrow = TRUE, nrow = 2, ncol = 2)
  
  x <- VAR.sim(B = A, lag = 1, include = "none", n = 1100)
  data_train = x[1:1000,]
  data_test = x[1001:1100,]
  
  
  #fit VAR model on the data
  var_model = VAR(data_train,p=VARselect(data_train)$selection[1]) # if you want to get the coefficients nicely use ar(data_train)
  
  #summary(var_model) to access the coefficients
  
  
  order_ar = c(order_ar, var_model$p)
  
  MAE_ar_1 = c(MAE_ar_1, mean(abs(residuals(var_model)[,1])))
  MAE_ar_2 = c(MAE_ar_2, mean(abs(residuals(var_model)[,2])))
  
  loglik_ar = c(loglik_ar, logLik(var_model))
  AIC_ar = c(AIC_ar, AIC(var_model))
  

  prediction_ar = predict(var_model, n.ahead = 100, ci=.95)
  
  MAE_test_ar_1 = c(MAE_test_ar_1, mean(abs(c(prediction_ar$fcst$y1[,1] - data_test[,1]))))
  MAE_test_ar_2 = c(MAE_test_ar_2, mean(abs(c(prediction_ar$fcst$y1[,2] - data_test[,2]))))
  
  
  inside_conf_int_test_ar_both = c(inside_conf_int_test_ar_both, length(which((data_test[,1]<prediction_ar$fcst$y1[,3])&
                                                                      (data_test[,1]>prediction_ar$fcst$y1[,2])&
                                                                      (data_test[,2]<prediction_ar$fcst$y2[,3])&
                                                                      (data_test[,2]>prediction_ar$fcst$y2[,2]))))
  inside_conf_int_test_ar_1 = c(inside_conf_int_test_ar_1, length(which((data_test[,1]<prediction_ar$fcst$y1[,3])&
                                                                                (data_test[,1]>prediction_ar$fcst$y1[,2]))))
  inside_conf_int_test_ar_2 = c(inside_conf_int_test_ar_2, length(which((data_test[,2]<prediction_ar$fcst$y2[,3])&
                                                                          (data_test[,2]>prediction_ar$fcst$y2[,2]))))
  
  vine_fit = svine(data_train, p=1)
  
  margin_vine_1 = c(margin_vine_1, lapply(summary(vine_fit)$margins[,3], as.character)[[1]])
  margin_vine_2 = c(margin_vine_2, lapply(summary(vine_fit)$margins[,3], as.character)[[2]])
  
  model_pred_med = array(c(NA,NA),dim = c(1,2))
  
  for (i in 2:1000){ #warning: takes long to run 
    vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = matrix(data_train[i,], ncol = 2))
    model_pred_med = rbind(model_pred_med,c(quantile(vine_prediction[,1,],  probs = c(0.5)),quantile(vine_prediction[,2,],  probs = c(0.5))))
  }
  
  model_pred_med = model_pred_med[-1,]
  
  MAE_vine_1 = c(MAE_vine_1, mean(abs(data_train[2:1000,1]- model_pred_med[,1])))
  MAE_vine_2 = c(MAE_vine_2, mean(abs(data_train[2:1000,2]- model_pred_med[,2])))
  
  loglik_vine = c(loglik_vine, logLik(vine_fit))
  AIC_vine = c(AIC_vine, AIC(vine_fit))
  
  vine_prediction <- svine_sim(n = 100, rep = 1000, model = vine_fit, past = data_train)
  vine_prediction_med_1 = c()
  vine_prediction_lwr_1 = c()
  vine_prediction_upr_1 = c()
  vine_prediction_med_2 = c()
  vine_prediction_lwr_2 = c()
  vine_prediction_upr_2 = c()
  
  for (j in 1:100){
    vine_prediction_med_1 = c(vine_prediction_med_1, quantile(vine_prediction[j,1 ,],  probs = c(0.5)))
    vine_prediction_lwr_1 = c(vine_prediction_lwr_1, quantile(vine_prediction[j,1 ,],  probs = c(0.025)))
    vine_prediction_upr_1 = c(vine_prediction_upr_1, quantile(vine_prediction[j,1 ,],  probs = c(0.975)))
    vine_prediction_med_2 = c(vine_prediction_med_2, quantile(vine_prediction[j,2 ,],  probs = c(0.5)))
    vine_prediction_lwr_2 = c(vine_prediction_lwr_2, quantile(vine_prediction[j,2 ,],  probs = c(0.025)))
    vine_prediction_upr_2 = c(vine_prediction_upr_2, quantile(vine_prediction[j,2 ,],  probs = c(0.975)))
  }
  
  MAE_test_vine_1 = c(MAE_test_vine_1, mean(abs(c(vine_prediction_med_1 - data_test[,1]))))
  MAE_test_vine_2 = c(MAE_test_vine_2, mean(abs(c(vine_prediction_med_2 - data_test[,2]))))

  inside_conf_int_test_vine_both = c(inside_conf_int_test_vine_both, length(which((data_test[,1]<vine_prediction_upr_1)&
                                                                          (data_test[,1]>vine_prediction_lwr_1)&
                                                                          (data_test[,2]<vine_prediction_upr_2)&
                                                                          (data_test[,2]>vine_prediction_lwr_2))))

  
  inside_conf_int_test_vine_1 = c(inside_conf_int_test_vine_1, length(which((data_test[,1]<vine_prediction_upr_1)&
                                                                                    (data_test[,1]>vine_prediction_lwr_1))))

  inside_conf_int_test_vine_2 = c(inside_conf_int_test_vine_2, length(which((data_test[,2]<vine_prediction_upr_2)&
                                                                              (data_test[,2]>vine_prediction_lwr_2))))
}


data.frame( MAE_1 = c(mean(MAE_ar_1), mean(MAE_vine_1)),
            MAE_1_sd = c(sd(MAE_ar_1), sd(MAE_vine_1)),
            MAE_2 = c(mean(MAE_ar_2), mean(MAE_vine_2)),
            MAE_2_sd = c(sd(MAE_ar_2), sd(MAE_vine_2)),
            LogLik = c(mean(loglik_ar), mean(loglik_vine)),
            LogLik_sd = c(sd(loglik_ar), sd(loglik_vine)),
            AIC = c(mean(AIC_ar), mean(AIC_vine)),
            AIC_sd = c(sd(AIC_ar), sd(AIC_vine)),
            MAE_test_1 = c(mean(MAE_test_ar_1), mean(MAE_test_vine_1)),
            MAE_test_1_sd = c(sd(MAE_test_ar_1), sd(MAE_test_vine_1)),
            MAE_test_2 = c(mean(MAE_test_ar_2), mean(MAE_test_vine_2)),
            MAE_test_2_sd = c(sd(MAE_test_ar_2), sd(MAE_test_vine_2)),
            inside_conf_int_test_both = c(mean(inside_conf_int_test_ar_both), mean(inside_conf_int_test_vine_both)),
            inside_conf_int_test_both_sd = c(sd(inside_conf_int_test_ar_both), sd(inside_conf_int_test_vine_both)),
            inside_conf_int_test_1 = c(mean(inside_conf_int_test_ar_1), mean(inside_conf_int_test_vine_1)),
            inside_conf_int_test_1_sd = c(sd(inside_conf_int_test_ar_1), sd(inside_conf_int_test_vine_1)),
            inside_conf_int_test_2 = c(mean(inside_conf_int_test_ar_2), mean(inside_conf_int_test_vine_2)),
            inside_conf_int_test_2_sd = c(sd(inside_conf_int_test_ar_2), sd(inside_conf_int_test_vine_2)),
            row.names = c("AR model", "vine model"))


order_ar



results_5 = data.frame(order_ar, 
                       MAE_ar_1,
                       MAE_ar_2,
                       loglik_ar,
                       AIC_ar, 
                       MAE_test_ar_1,
                       MAE_test_ar_2,
                       inside_conf_int_test_ar_1,
                       inside_conf_int_test_ar_2,
                       inside_conf_int_test_ar_both,
                       margin_vine_1,
                       margin_vine_2,
                       MAE_vine_1,
                       MAE_vine_2,
                       loglik_vine,
                       AIC_vine,
                       MAE_test_vine_1,
                       MAE_test_vine_2,
                       inside_conf_int_test_vine_both,
                       inside_conf_int_test_vine_1,
                       inside_conf_int_test_vine_2)


data.frame( MAE_1 = c(mean(results_5$MAE_ar_1), mean(results_5$MAE_vine_1)),
            MAE_1_sd = c(sd(results_5$MAE_ar_1), sd(results_5$MAE_vine_1)),
            MAE_2 = c(mean(results_5$MAE_ar_2), mean(results_5$MAE_vine_2)),
            MAE_2_sd = c(sd(results_5$MAE_ar_2), sd(results_5$MAE_vine_2)),
            LogLik = c(mean(results_5$loglik_ar), mean(results_5$loglik_vine)),
            LogLik_sd = c(sd(results_5$loglik_ar), sd(results_5$loglik_vine)),
            AIC = c(mean(results_5$AIC_ar), mean(results_5$AIC_vine)),
            AIC_sd = c(sd(results_5$AIC_ar), sd(results_5$AIC_vine)),
            MAE_test_1 = c(mean(results_5$MAE_test_ar_1), mean(results_5$MAE_test_vine_1)),
            MAE_test_1_sd = c(sd(results_5$MAE_test_ar_1), sd(results_5$MAE_test_vine_1)),
            MAE_test_2 = c(mean(results_5$MAE_test_ar_2), mean(results_5$MAE_test_vine_2)),
            MAE_test_2_sd = c(sd(results_5$MAE_test_ar_2), sd(results_5$MAE_test_vine_2)),
            inside_conf_int_test_both = c(mean(results_5$inside_conf_int_test_ar_both), mean(results_5$inside_conf_int_test_vine_both)),
            inside_conf_int_test_both_sd = c(sd(results_5$inside_conf_int_test_ar_both), sd(results_5$inside_conf_int_test_vine_both)),
            inside_conf_int_test_1 = c(mean(results_5$inside_conf_int_test_ar_1), mean(results_5$inside_conf_int_test_vine_1)),
            inside_conf_int_test_1_sd = c(sd(results_5$inside_conf_int_test_ar_1), sd(results_5$inside_conf_int_test_vine_1)),
            inside_conf_int_test_2 = c(mean(results_5$inside_conf_int_test_ar_2), mean(results_5$inside_conf_int_test_vine_2)),
            inside_conf_int_test_2_sd = c(sd(results_5$inside_conf_int_test_ar_2), sd(results_5$inside_conf_int_test_vine_2)),
            row.names = c("AR model", "vine model"))



results_5$order_ar

length(which(results_5$AIC_ar>results_5$AIC_vine)) # 100%
length(which(results_5$loglik_ar>results_5$loglik_vine)) # 90%

length(which(results_5$margin_vine_1=="Normal")) # 48% 
length(which(results_5$margin_vine_2=="Normal")) # 30% 
length(which((results_5$margin_vine_1=="Normal")&(results_5$margin_vine_2=="Normal"))) #12%

mean(results_5$order_ar) #2.7
sd(results_5$order_ar) # 0.274
length(which(results_5$order_ar>1)) #8%


#################################     model 6    ##############################################


order_ar = c()
MAE_ar_1 = c()
MAE_ar_2 = c()
loglik_ar = c()
AIC_ar = c() 
MAE_test_ar_1 = c()
MAE_test_ar_2 = c()
inside_conf_int_test_ar_1 = c()
inside_conf_int_test_ar_2 = c()
inside_conf_int_test_ar_both = c()
margin_vine_1 = c()
margin_vine_2 = c()
MAE_vine_1 = c()
MAE_vine_2 = c()
loglik_vine = c()
AIC_vine = c()
MAE_test_vine_1 = c()
MAE_test_vine_2 = c()
inside_conf_int_test_vine_both = c()
inside_conf_int_test_vine_1 = c()
inside_conf_int_test_vine_2 = c()


for (i in 1:50){
  print(i)
  #simulate the data from vine process
 
  univ1 <- univ2 <- univariateML::mlnorm(rnorm(10))
  
  # modify the parameters to N(5, 10) and N(0, 2) distributions
  univ1[] <- c(5, 10)
  univ2[] <- c(0, 2)
  
  cs_struct <- cvine_structure(1:2)
  
  pcs <- list(
    list(  # first tree
      bicop_dist("clayton", 0, 3), # cross sectional copula
      bicop_dist("gumbel", 0, 3)  # serial copula
    ),
    list(  # second tree
      bicop_dist("frank", 0, 15), 
      bicop_dist("indep")  
    ),
    list( # third tree
      bicop_dist("clayton", 0, 5)
    )
  )
  
  cop <- svinecop_dist(pcs, cs_struct, p = 1, out_vertices = 1:2, in_vertices = 1:2)
  
  model <- svine_dist(margins = list(univ1, univ2), copula = cop)
  x <- svine_sim(n = 1100, rep = 1, model = model)
  
  data_train = x[1:1000,]
  data_test = x[1001:1100,]
  
  
  #fit VAR model on the data
  var_model = VAR(data_train,p=VARselect(data_train)$selection[1]) # if you want to get the coefficients nicely use ar(data_train)
  
  #summary(var_model) to access the coefficients
  
  
  order_ar = c(order_ar, var_model$p)
  
  MAE_ar_1 = c(MAE_ar_1, mean(abs(residuals(var_model)[,1])))
  MAE_ar_2 = c(MAE_ar_2, mean(abs(residuals(var_model)[,2])))
  
  loglik_ar = c(loglik_ar, logLik(var_model))
  AIC_ar = c(AIC_ar, AIC(var_model))
  
  
  prediction_ar = predict(var_model, n.ahead = 100, ci=.95)
  
  MAE_test_ar_1 = c(MAE_test_ar_1, mean(abs(c(prediction_ar$fcst$y1[,1] - data_test[,1]))))
  MAE_test_ar_2 = c(MAE_test_ar_2, mean(abs(c(prediction_ar$fcst$y1[,2] - data_test[,2]))))
  
  
  inside_conf_int_test_ar_both = c(inside_conf_int_test_ar_both, length(which((data_test[,1]<prediction_ar$fcst$y1[,3])&
                                                                                (data_test[,1]>prediction_ar$fcst$y1[,2])&
                                                                                (data_test[,2]<prediction_ar$fcst$y2[,3])&
                                                                                (data_test[,2]>prediction_ar$fcst$y2[,2]))))
  inside_conf_int_test_ar_1 = c(inside_conf_int_test_ar_1, length(which((data_test[,1]<prediction_ar$fcst$y1[,3])&
                                                                          (data_test[,1]>prediction_ar$fcst$y1[,2]))))
  inside_conf_int_test_ar_2 = c(inside_conf_int_test_ar_2, length(which((data_test[,2]<prediction_ar$fcst$y2[,3])&
                                                                          (data_test[,2]>prediction_ar$fcst$y2[,2]))))
  
  vine_fit = svine(data_train, p=1)
  
  margin_vine_1 = c(margin_vine_1, lapply(summary(vine_fit)$margins[,3], as.character)[[1]])
  margin_vine_2 = c(margin_vine_2, lapply(summary(vine_fit)$margins[,3], as.character)[[2]])
  
  model_pred_med = array(c(NA,NA),dim = c(1,2))
  
  for (i in 2:1000){ #warning: takes long to run 
    vine_prediction <- svine_sim(n = 1, rep = 1000, model = vine_fit, past = matrix(data_train[i,], ncol = 2))
    model_pred_med = rbind(model_pred_med,c(quantile(vine_prediction[,1,],  probs = c(0.5)),quantile(vine_prediction[,2,],  probs = c(0.5))))
  }
  
  model_pred_med = model_pred_med[-1,]
  
  MAE_vine_1 = c(MAE_vine_1, mean(abs(data_train[2:1000,1]- model_pred_med[,1])))
  MAE_vine_2 = c(MAE_vine_2, mean(abs(data_train[2:1000,2]- model_pred_med[,2])))
  
  loglik_vine = c(loglik_vine, logLik(vine_fit))
  AIC_vine = c(AIC_vine, AIC(vine_fit))
  
  vine_prediction <- svine_sim(n = 100, rep = 1000, model = vine_fit, past = data_train)
  vine_prediction_med_1 = c()
  vine_prediction_lwr_1 = c()
  vine_prediction_upr_1 = c()
  vine_prediction_med_2 = c()
  vine_prediction_lwr_2 = c()
  vine_prediction_upr_2 = c()
  
  for (j in 1:100){
    vine_prediction_med_1 = c(vine_prediction_med_1, quantile(vine_prediction[j,1 ,],  probs = c(0.5)))
    vine_prediction_lwr_1 = c(vine_prediction_lwr_1, quantile(vine_prediction[j,1 ,],  probs = c(0.025)))
    vine_prediction_upr_1 = c(vine_prediction_upr_1, quantile(vine_prediction[j,1 ,],  probs = c(0.975)))
    vine_prediction_med_2 = c(vine_prediction_med_2, quantile(vine_prediction[j,2 ,],  probs = c(0.5)))
    vine_prediction_lwr_2 = c(vine_prediction_lwr_2, quantile(vine_prediction[j,2 ,],  probs = c(0.025)))
    vine_prediction_upr_2 = c(vine_prediction_upr_2, quantile(vine_prediction[j,2 ,],  probs = c(0.975)))
  }
  
  MAE_test_vine_1 = c(MAE_test_vine_1, mean(abs(c(vine_prediction_med_1 - data_test[,1]))))
  MAE_test_vine_2 = c(MAE_test_vine_2, mean(abs(c(vine_prediction_med_2 - data_test[,2]))))
  
  inside_conf_int_test_vine_both = c(inside_conf_int_test_vine_both, length(which((data_test[,1]<vine_prediction_upr_1)&
                                                                                    (data_test[,1]>vine_prediction_lwr_1)&
                                                                                    (data_test[,2]<vine_prediction_upr_2)&
                                                                                    (data_test[,2]>vine_prediction_lwr_2))))
  
  
  inside_conf_int_test_vine_1 = c(inside_conf_int_test_vine_1, length(which((data_test[,1]<vine_prediction_upr_1)&
                                                                              (data_test[,1]>vine_prediction_lwr_1))))
  
  inside_conf_int_test_vine_2 = c(inside_conf_int_test_vine_2, length(which((data_test[,2]<vine_prediction_upr_2)&
                                                                              (data_test[,2]>vine_prediction_lwr_2))))
  
}


data.frame( MAE_1 = c(mean(MAE_ar_1), mean(MAE_vine_1)),
            MAE_1_sd = c(sd(MAE_ar_1), sd(MAE_vine_1)),
            MAE_2 = c(mean(MAE_ar_2), mean(MAE_vine_2)),
            MAE_2_sd = c(sd(MAE_ar_2), sd(MAE_vine_2)),
            LogLik = c(mean(loglik_ar), mean(loglik_vine)),
            LogLik_sd = c(sd(loglik_ar), sd(loglik_vine)),
            AIC = c(mean(AIC_ar), mean(AIC_vine)),
            AIC_sd = c(sd(AIC_ar), sd(AIC_vine)),
            MAE_test_1 = c(mean(MAE_test_ar_1), mean(MAE_test_vine_1)),
            MAE_test_1_sd = c(sd(MAE_test_ar_1), sd(MAE_test_vine_1)),
            MAE_test_2 = c(mean(MAE_test_ar_2), mean(MAE_test_vine_2)),
            MAE_test_2_sd = c(sd(MAE_test_ar_2), sd(MAE_test_vine_2)),
            inside_conf_int_test_both = c(mean(inside_conf_int_test_ar_both), mean(inside_conf_int_test_vine_both)),
            inside_conf_int_test_both_sd = c(sd(inside_conf_int_test_ar_both), sd(inside_conf_int_test_vine_both)),
            inside_conf_int_test_1 = c(mean(inside_conf_int_test_ar_1), mean(inside_conf_int_test_vine_1)),
            inside_conf_int_test_1_sd = c(sd(inside_conf_int_test_ar_1), sd(inside_conf_int_test_vine_1)),
            inside_conf_int_test_2 = c(mean(inside_conf_int_test_ar_2), mean(inside_conf_int_test_vine_2)),
            inside_conf_int_test_2_sd = c(sd(inside_conf_int_test_ar_2), sd(inside_conf_int_test_vine_2)),
            row.names = c("AR model", "vine model"))


order_ar



results_6 = data.frame(order_ar, 
                       MAE_ar_1,
                       MAE_ar_2,
                       loglik_ar,
                       AIC_ar, 
                       MAE_test_ar_1,
                       MAE_test_ar_2,
                       inside_conf_int_test_ar_1,
                       inside_conf_int_test_ar_2,
                       inside_conf_int_test_ar_both,
                       margin_vine_1,
                       margin_vine_2,
                       MAE_vine_1,
                       MAE_vine_2,
                       loglik_vine,
                       AIC_vine,
                       MAE_test_vine_1,
                       MAE_test_vine_2,
                       inside_conf_int_test_vine_both,
                       inside_conf_int_test_vine_1,
                       inside_conf_int_test_vine_2)

write_xlsx(results_6,"C:\\Users\\laura\\Dropbox\\TU Delft Msc\\Uncertainty and sensitivity analysis WI4050\\results_6.xlsx")


plot(x[1000:1100,1], type='l', ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds", ylim = c(-30,30))
lines(c(1:100), prediction_ar$fcst$y1[,1], col='blue')
lines(c(1:100), prediction_ar$fcst$y1[,2], col='blue')
lines(c(1:100), prediction_ar$fcst$y1[,3], col='blue')
lines(c(1:100), vine_prediction_lwr_1, col='red')
lines(c(1:100), vine_prediction_med_1, col='red')
lines(c(1:100), vine_prediction_upr_1, col='red')
legend("bottomright", 
       legend = c("Data", "AR model", "Vine model"), 
       col =  c("black", "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )


plot(x[1000:1100,2], type='l', ylab = expression('x'[t]), xlab = "Time", 
     main = "Prediction on the test set with 95% confidence bounds", ylim = c(-7,6))
lines(c(1:100), prediction_ar$fcst$y2[,1], col='blue')
lines(c(1:100), prediction_ar$fcst$y2[,2], col='blue')
lines(c(1:100), prediction_ar$fcst$y2[,3], col='blue')
lines(c(1:100), vine_prediction_lwr_2, col='red')
lines(c(1:100), vine_prediction_med_2, col='red')
lines(c(1:100), vine_prediction_upr_2, col='red')
legend("bottomright", 
       legend = c("Data", "AR model", "Vine model"), 
       col =  c("black", "blue", "red"), 
       pt.cex = 1, 
       lty=c(1,1),
       cex = 0.8, 
       text.col = "black", 
       horiz = F )

data.frame( MAE_1 = c(mean(results_6$MAE_ar_1), mean(results_6$MAE_vine_1)),
            MAE_1_sd = c(sd(results_6$MAE_ar_1), sd(results_6$MAE_vine_1)),
            MAE_2 = c(mean(results_6$MAE_ar_2), mean(results_6$MAE_vine_2)),
            MAE_2_sd = c(sd(results_6$MAE_ar_2), sd(results_6$MAE_vine_2)),
            LogLik = c(mean(results_6$loglik_ar), mean(results_6$loglik_vine)),
            LogLik_sd = c(sd(results_6$loglik_ar), sd(results_6$loglik_vine)),
            AIC = c(mean(results_6$AIC_ar), mean(results_6$AIC_vine)),
            AIC_sd = c(sd(results_6$AIC_ar), sd(results_6$AIC_vine)),
            MAE_test_1 = c(mean(results_6$MAE_test_ar_1), mean(results_6$MAE_test_vine_1)),
            MAE_test_1_sd = c(sd(results_6$MAE_test_ar_1), sd(results_6$MAE_test_vine_1)),
            MAE_test_2 = c(mean(results_6$MAE_test_ar_2), mean(results_6$MAE_test_vine_2)),
            MAE_test_2_sd = c(sd(results_6$MAE_test_ar_2), sd(results_6$MAE_test_vine_2)),
            inside_conf_int_test_both = c(mean(results_6$inside_conf_int_test_ar_both), mean(results_6$inside_conf_int_test_vine_both)),
            inside_conf_int_test_both_sd = c(sd(results_6$inside_conf_int_test_ar_both), sd(results_6$inside_conf_int_test_vine_both)),
            inside_conf_int_test_1 = c(mean(results_6$inside_conf_int_test_ar_1), mean(results_6$inside_conf_int_test_vine_1)),
            inside_conf_int_test_1_sd = c(sd(results_6$inside_conf_int_test_ar_1), sd(results_6$inside_conf_int_test_vine_1)),
            inside_conf_int_test_2 = c(mean(results_6$inside_conf_int_test_ar_2), mean(results_6$inside_conf_int_test_vine_2)),
            inside_conf_int_test_2_sd = c(sd(results_6$inside_conf_int_test_ar_2), sd(results_6$inside_conf_int_test_vine_2)),
            row.names = c("AR model", "vine model"))


data.frame(MAE_test = c(mean(results_6$MAE_test_ar_1), mean(results_6$MAE_test_ar_2), mean(results_6$MAE_test_vine_1), mean(results_6$MAE_test_vine_2)),
           MAE_test_sd = c(sd(results_6$MAE_test_ar_1), sd(results_6$MAE_test_ar_2), sd(results_6$MAE_test_vine_1), sd(results_6$MAE_test_vine_2)))

data.frame(MAE_train = c(mean(results_6$MAE_ar_1), mean(results_6$MAE_ar_2), mean(results_6$MAE_vine_1), mean(results_6$MAE_vine_2)),
           MAE_train_sd = c(sd(results_6$MAE_ar_1), sd(results_6$MAE_ar_2), sd(results_6$MAE_vine_1), sd(results_6$MAE_vine_2)))



length(which(results_6$AIC_ar>results_6$AIC_vine)) # 100%
length(which(results_6$loglik_ar>results_6$loglik_vine)) # 90%

length(which(results_6$margin_vine_1=="Normal")) # 28% 
length(which(results_6$margin_vine_2=="Normal")) # 40% 
length(which((results_6$margin_vine_1=="Normal")&(results_6$margin_vine_2=="Normal"))) #20%

mean(results_6$order_ar) #2.7
sd(results_6$order_ar) # 0.274
length(which(results_6$order_ar>1)) #74%



