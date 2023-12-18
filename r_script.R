
######################################## Appendix Script ########################################

# R-Script for:
# Interrupted Time Series: A Comparison Of Different Methods With Theoretical And Practical Considerations
# Uni Assignment For -Data Science In Clinical Practice-
# Henrik Godmann (13802453)

######################################## Prepare data ########################################

source("/Users/henrikgodmann/Desktop/workspace/GitHub/functions/colors/Rcolors.R")

dat <- read.csv("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/ESMdata/ESMdata.csv")

dat$date <- as.Date(dat$date, format = "%d/%m/%y") 

library(missForest)

dat1 <- subset(dat, dat$phase < 4)
dat2 <- subset(dat, dat$phase >3)


numeric_vars1 <- dat1[sapply(dat1, function(x) is.numeric(x) || is.integer(x))]
set.seed(42)
imputed <- missForest(numeric_vars1)
imputed2 <- round(imputed$ximp,4)
imputed2_cut <- imputed2[c(2,82)]
dat1 <- cbind(dat1[2], imputed2_cut)


numeric_vars2 <- dat2[sapply(dat2, function(x) is.numeric(x) || is.integer(x))]
set.seed(42)
imputed <- missForest(numeric_vars2)
imputed3 <- round(imputed$ximp,4)
imputed3_cut <- imputed3[c(2,82)]
dat2 <- cbind(dat2[2], imputed3_cut)


dat3 <- rbind(dat1,dat2)


# which phases do we have?
table(dat3$phase)
# 1 = baseline
# 2 = double blind before reducing medication
# 3 = double blind during medication reduction
# 4 = phase after medication reduction
# 5 = phase after experiment


######################################## ARIMA ########################################

# Steps following Schaffer et al. (2021)

#### Plot data to understand patterns

# Plot the time series to understand the patterns, specifically pre-existing trends, seasonal effects
library(ggplot2)
library(tseries)
library(gridExtra)

dat <- dat[c(2,86)]
dat <- na.omit(dat)


# setwd("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/plots")
# png("imputation.png", width=12*300, height=8*300, res=300)

# 'dat3' contains imputed data and 'dat' contains original data.

y_range <- range(dat3$dep, dat$dep, na.rm = TRUE)

plot(dat3$date, dat3$dep, type = "l", col = my_red, lwd = 2,
     xlim = range(dat3$date, dat$date), ylim = y_range,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")


points(dat$date, dat$dep, col = my_blue, pch = 16, cex = 1.5) 

title(main = "Depression Score Over Time", font.main = 2, cex.main = 2)
title(xlab = "Time", cex.lab = 1.5)
title(ylab = "Depression Score", cex.lab = 1.5)

month_dates <- seq(from = min(dat3$date), to = max(dat3$date), by = "month")
axis(1, at = month_dates, labels = format(month_dates, "%b %Y"), cex.axis = 1.2, las = 1)

axis(2, las = 1, cex.axis = 1.2)

box(lwd=2)

# dev.off()

# looks non-stationairy with a change in mean and variance (and also autocorrelation, see Wichers et al., 2016)
# likely no seasonal effects in this time window due to data structure and non seasonal length of data collection.
# However, longer data is needed to actually account for seasonal effects, which in depression might actually occur
# (for instance in winter)

# from now on, we will model the data before interruption to allow for a later comparison


library(dplyr)


# Group the data by date and phase, then summarize it by taking the mean of dep; run only for second analysis
dat1 <- dat1 %>%
  group_by(date, phase) %>%
  summarise(mean_dep = mean(dep, na.rm = TRUE))

colnames(dat1)[3] <- "dep"


dat2 <- dat2 %>%
  group_by(date, phase) %>%
  summarise(mean_dep = mean(dep, na.rm = TRUE))

colnames(dat2)[3] <- "dep"

# We can also test stationarity with a hypothesis test: (Dickey-Fuller test)
adf.test(dat1$dep, alternative = "stationary")
# significant with all measures, not significant with daily measures, so the time series rather needs differencing
# -> first order difference to account for the trend (d = 1)

# check auto correlation plot
acf(dat1$dep, main="ACF Plot")

# positive autocorrelation, so AR term is needed. but as this is still the non-stationairy time series, let's use
# an automated algorithm to be sure to select an appropriate model

library(forecast)

# Fit Arima
fit <- auto.arima(dat1$dep, d = 1, seasonal = FALSE)
summary(fit)

forecasts <- forecast(fit, h = nrow(dat2))


# Prepare data for plotting all data (skip if wanting daily data)
date <- c(dat1$date[450:671], dat2$date[1:150])
actual_data <- c(dat1$dep[450:671], dat2$dep[1:150])
predicted <- c(rep(NA, length(dat1$dep[450:671])), forecasts$mean[1:150])
lower_bound <- c(rep(NA, length(dat1$dep[450:671])), forecasts$lower[, "95%"][1:150])
upper_bound <- c(rep(NA, length(dat1$dep[450:671])), forecasts$upper[, "95%"][1:150])

# prepare to plot daily data
date <- c(dat1$date[60:98], dat2$date[1:20])
actual_data <- c(dat1$dep[60:98], dat2$dep[1:20])
predicted <- c(rep(NA, length(dat1$dep[60:98])), forecasts$mean[1:20])
lower_bound <- c(rep(NA, length(dat1$dep[60:98])), forecasts$lower[, "95%"][1:20])
upper_bound <- c(rep(NA, length(dat1$dep[60:98])), forecasts$upper[, "95%"][1:20])


df <- data.frame(date, actual_data, predicted, lower_bound, upper_bound)

# setwd("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/plots")
# png("arima3.png", width=12*300, height=8*300, res=300)
# Plot
plot(df$date, df$actual_data, type = "l", col = my_blue, lwd = 2,
     xlab = "", ylab = "",
     main = "",
     ylim = c(0, 4),  # Set y-axis limits
     cex.axis = 1.2, las = 1,cex.lab = 1.5)

title(main = "Actual vs. Forecasted Depression Score ARIMA", font.main = 2, cex.main = 1.5)
title(xlab = "Time", cex.lab = 1.5)
title(ylab = "Depression Score", cex.lab = 1.5)

# Add predicted data with thicker lines
lines(df$date, df$predicted, col = my_red, lwd = 3)  # Thicker line for predicted values

# Add confidence intervals
lines(df$date, df$lower_bound, col = my_red_strong, lwd = 2, lty = 2)  # Lower bound of confidence interval
lines(df$date, df$upper_bound, col = my_red_strong, lwd = 2, lty = 2)  # Upper bound of confidence interval

# Add a legend
legend("topright", legend = c("Actual Data", "Predicted", "80% Confidence Interval"),
       col = c(my_blue, my_red, my_red_strong), lty = c(1, 1, 2), cex = 0.8, lwd = c(2, 3, 2))

# Draw a thicker box around the plot
box(lwd = 2)

# Add a horizontal dashed line at the date 2012-11-19
abline(v=as.Date("2012-11-19"), lty=2, lwd=1.5)

# dev.off()


####### mood predition (see Appendix)

imputed2$mood_mean <- rowMeans(imputed2[c(8,9,11,12,14,16,17)])
names(imputed2)

imputedM1 <- imputed2[c(83)]
datM1 <- cbind(dat1[1], imputedM1)


imputed3$mood_mean <- rowMeans(imputed3[c(8,9,11,12,14,16,17)])
names(imputed3)

imputedM2 <- imputed3[c(83)]
datM2 <- cbind(dat2[1], imputedM2)

datM3 <- rbind(datM1,datM2)



datM1 <- datM1 %>%
  group_by(date) %>%
  summarise(mood_mean = mean(mood_mean, na.rm = TRUE))

datM2 <- datM2 %>%
  group_by(date) %>%
  summarise(mood_mean = mean(mood_mean, na.rm = TRUE))


plot(datM3$date, datM3$mood_mean, type = "l", col = my_blue, lwd = 2,
     xlim = range(datM3$date), ylim = range(datM3$mood_mean),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")

title(main = "Negative Mood Score Over Time", font.main = 2, cex.main = 2)
title(xlab = "Time", cex.lab = 1.5)
title(ylab = "Mood Score", cex.lab = 1.5)

month_dates <- seq(from = min(datM3$date), to = max(datM3$date), by = "month")
axis(1, at = month_dates, labels = format(month_dates, "%b %Y"), cex.axis = 1.2, las = 1)

axis(2, las = 1, cex.axis = 1.2)
box(lwd=2)


#fit autom- arima to mood data
fit2 <- auto.arima(datM1$mood_mean, d = 0, seasonal = FALSE)
summary(fit2)

forecasts <- forecast(fit2, h = nrow(datM2))

# Prepare data for plotting all (skip if you want to plot daily)
date <- c(datM1$date[500:671], datM2$date[1:30])
actual_data <- c(datM1$mood_mean[500:671], datM2$mood_mean[1:30])
predicted <- c(rep(NA, length(datM1$mood_mean[500:671])), forecasts$mean[1:30])
lower_bound <- c(rep(NA, length(datM1$mood_mean[500:671])), forecasts$lower[, "80%"][1:30])
upper_bound <- c(rep(NA, length(datM1$mood_mean[500:671])), forecasts$upper[, "80%"][1:30])

# Prepare data for plotting daily
date <- c(datM1$date[60:98], datM2$date[1:10])
actual_data <- c(datM1$mood_mean[60:98], datM2$mood_mean[1:10])
predicted <- c(rep(NA, length(datM1$mood_mean[60:98])), forecasts$mean[1:10])
lower_bound <- c(rep(NA, length(datM1$mood_mean[60:98])), forecasts$lower[, "80%"][1:10])
upper_bound <- c(rep(NA, length(datM1$mood_mean[60:98])), forecasts$upper[, "80%"][1:10])



df <- data.frame(date, actual_data, predicted, lower_bound, upper_bound)

# setwd("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/plots")
# png("arima_mood.png", width=12*300, height=8*300, res=300)
# Plot
plot(df$date, df$actual_data, type = "l", col = my_blue, lwd = 2,
     xlab = "", ylab = "",
     main = "",
     ylim = c(0, 4),  # Set y-axis limits
     cex.axis = 1.2, las = 1,cex.lab = 1.5)

title(main = "Actual vs. Forecasted Mood Score ARIMA", font.main = 2, cex.main = 1.5)
title(xlab = "Time", cex.lab = 1.5)
title(ylab = "Mood Score", cex.lab = 1.5)

lines(df$date, df$predicted, col = my_red, lwd = 3) 

lines(df$date, df$lower_bound, col = my_red_strong, lwd = 2, lty = 2)  
lines(df$date, df$upper_bound, col = my_red_strong, lwd = 2, lty = 2)  

legend("topleft", legend = c("Actual Data", "Predicted", "80% Confidence Interval"),
       col = c(my_blue, my_red, my_red_strong), lty = c(1, 1, 2), cex = 0.8, lwd = c(2, 3, 2))
box(lwd = 2)
abline(v=as.Date("2012-11-19"), lty=2, lwd=1.5)

# dev.off()



#################################### Network Analysis ####################################

####### CLEAR ENVIRONMENT AND RUN PREPARE DATA SECTION AGAIN FIRST


library(mlVAR)
library(graphicalVAR)
library(psychonetrics)
library(dplyr)
library(qgraph)



# I will now estimate a saturated GVAR model with psychonetrics. I will plot the estimated temporal (partial
# directed correlations) and contemporaneous networks

imputed2_cut <- imputed2[c(8,11,12,15,18,82)]
datNT1 <- cbind(dat1[1], imputed2_cut)

imputed3_cut <- imputed3[c(8,11,12,15,18,82)]
datNT2 <- cbind(dat2[1], imputed3_cut)

vars_all <- colnames(datNT1[c(2:7)])


# check assumption of stationairity and difference data
for (var in vars_all) {
  print(adf.test(datNT1[[var]], alternative = "stationary"))
}


for (var in vars_all) {
  print(adf.test(datNT1[[var]], alternative = "stationary"))
}


datNT1 <- datNT1 %>%
  mutate(across(starts_with("mood"), ~c(NA, diff(.)))) %>%
  mutate(dep_diff = c(NA, diff(dep)))

datNT1 <- datNT1[-1, ]

datNT2 <- datNT2 %>%
  mutate(across(starts_with("mood"), ~c(NA, diff(.)))) %>%
  mutate(dep_diff = c(NA, diff(dep)))

datNT2 <- datNT2[-1, ]



# saturated gaussian graphical VAR model
sat_model1 <- gvar(datNT1,  dayvar = "date", vars = vars_all, estimator = "FIML")
sat_model1 <- sat_model1 %>% runmodel()# estimation
cont_sat_model1 <- sat_model1 %>% getmatrix("omega_zeta")# Get contemporaneous network
temp_sat_model1 <- sat_model1 %>% getmatrix("PDC")# Get temporal network


# saturated gaussian graphical VAR model
sat_model2 <- gvar(datNT2,  dayvar = "date", vars = vars_all, estimator = "FIML")
sat_model2 <- sat_model2 %>% runmodel()# estimation
cont_sat_model2 <- sat_model2 %>% getmatrix("omega_zeta")# Get contemporaneous network
temp_sat_model2 <- sat_model2 %>% getmatrix("PDC")# Get temporal network


# Write own function to plot networks and specify desired outcomes
plot_network1 <- function(model_temp, model_contemp, data, dayvar, vars, maximum_edge_weight, cut_value, layout=NULL) {

  max_edge_weight <- maximum_edge_weight

  cut_value <- cut_value
  qgraph(model_temp, labels = vars, theme = "colorblind", layout = layout,
         maximum = max_edge_weight, cut = cut_value,vsize = 15,  label.cex = 1.5)

  mtext("A", side = 3, line = 1, adj = 0, cex = 1, font = 2.5)

  qgraph(model_contemp, labels = vars, theme = "colorblind", layout = layout,
         maximum = max_edge_weight, cut = cut_value,vsize = 15,  label.cex = 1.5)

  mtext("B", side = 3, line = 1, adj = 0, cex = 1, font = 2.5)

  # sat_model %>% parameters
}

plot_network2 <- function(model_temp, model_contemp, data, dayvar, vars, maximum_edge_weight, cut_value, layout=NULL) {

  max_edge_weight <- maximum_edge_weight

  cut_value <- cut_value
  qgraph(model_temp, labels = vars, theme = "colorblind", layout = layout,
         maximum = max_edge_weight, cut = cut_value,vsize = 15,  label.cex = 1.5)

  mtext("C", side = 3, line = 1, adj = 0, cex = 1, font = 2.5)

  qgraph(model_contemp, labels = vars, theme = "colorblind", layout = layout,
         maximum = max_edge_weight, cut = cut_value,vsize = 15,  label.cex = 1.5)

  mtext("D", side = 3, line = 1, adj = 0, cex = 1, font = 2.5)
  
  # sat_model %>% parameters
}

# for comparability
max_edge_weight <- 0.5
cut_value <- 0.5


# setwd("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/plots")
# png("networks2.png", width=8*300, height=8*300, res=300)
par(mfrow = c(2,2))
plot_network1(cont_sat_model1,temp_sat_model1, datNT1, dayvar =datNT1$date , vars = vars_all, max_edge_weight, cut_value)

plot_network2(cont_sat_model2,temp_sat_model2,datNT2, dayvar = datNT2$date, vars = vars_all, max_edge_weight, cut_value)
# dev.off()



#################################### Random Effects Multilevel Modeling ####################################

####### CLEAR ENVIRONMENT AND RUN PREPARE DATA SECTION AGAIN FIRST

### prepare data for e-clip

library(lme4)

imputed2_cut <- imputed2[c(2,82)]
datML1 <- cbind(dat1[1], imputed2_cut)
datML1$phase <- 0

imputed3_cut <- imputed3[c(2,82)]
datML2 <- cbind(dat2[1], imputed3_cut)
datML2$phase <- 1

datML <- rbind(datML1,datML2)

write.csv(datML,"/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/data/data.csv")


#################################### State Space Models ####################################

####### CLEAR ENVIRONMENT AND RUN PREPARE DATA SECTION AGAIN FIRST

imputed2_cut <- imputed2[c(2,69:82)]
datSSM1 <- cbind(dat1[1], imputed2_cut)
datSSM1$phase <- 0

imputed3_cut <- imputed3[c(2,69:82)]
datSSM2 <- cbind(dat2[1], imputed3_cut)
datSSM2$phase <- 1

datSSM <- rbind(datSSM1,datSSM2)

head(datSSM)

library(dplyr)
library(KFAS)

# For State Space Models, we will aggregate the data to daily measures
daily_data <- datSSM %>%
  group_by(date) %>%
  summarise(across(starts_with("dep"), mean, na.rm = TRUE),
            phase = mean(phase, na.rm = TRUE))


ssm_model <- SSModel(daily_data$dep ~ SSMtrend(1, Q = list(matrix(NA))) + 
                       SSMregression(~ phase, data = daily_data, Q = matrix(NA)), 
                     H = matrix(NA))


fit <- fitSSM(ssm_model, inits = rep(0, length = 3), method = "BFGS")# Fit the model
fit

optimized_parameters <- fit$optim.out$par
print(optimized_parameters)

kfs_results <- KFS(fit$model)
smoothed_states <- kfs_results$a
print(smoothed_states)

df <- as.data.frame(smoothed_states)


# setwd("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/plots")
# png("phase_ssm.png", width=12*300, height=8*300, res=300)
plot(daily_data$date, df$phase[-c(1)], type = "l", col = my_red, lwd = 3,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")

title(xlab = "Time", cex.lab = 1)
title(ylab = "State Estimate", cex.lab = 1)
month_dates <- seq(from = min(dat3$date), to = max(dat3$date), by = "month")
axis(1, at = month_dates, labels = format(month_dates, "%b %Y"), cex.axis = 1.2, las = 1)
axis(2, las = 1, cex.axis = 1.2)
box(lwd=2)
# dev.off()


# setwd("/Users/henrikgodmann/Desktop/workspace/semester_5/Clinical_Data_Science/final_assignment/plots")
# png("ssm_results.png", width=8*300, height=8*300, res=300)
par(mfrow = c(2,1))

plot(daily_data$date, daily_data$dep, type = "l", col = my_red, lwd = 2,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")

title(main = "Depression Score An Estimated States Over Time", font.main = 2, cex.main = 1)
title(xlab = "Time", cex.lab = 1)
title(ylab = "Depression Score", cex.lab = 1)
month_dates <- seq(from = min(dat3$date), to = max(dat3$date), by = "month")
axis(1, at = month_dates, labels = format(month_dates, "%b %Y"), cex.axis = 1.2, las = 1)
axis(2, las = 1, cex.axis = 1.2)
box(lwd=2)

plot(daily_data$date, smoothed_states[-1], type = "l", col = my_red, lwd = 2,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")

# title(main = "Depression Score Over Time", font.main = 2, cex.main = 1)
title(xlab = "Time", cex.lab = 1)
title(ylab = "State Estimate", cex.lab = 1)
month_dates <- seq(from = min(dat3$date), to = max(dat3$date), by = "month")
axis(1, at = month_dates, labels = format(month_dates, "%b %Y"), cex.axis = 1.2, las = 1)
axis(2, las = 1, cex.axis = 1.2)
box(lwd=2)


abline(v = daily_data$date[99], col = my_blue, lty = 2, lwd=2) 


# dev.off()








