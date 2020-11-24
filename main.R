library(fda.usc)
library(ftsa)
library(fda)
library(rainbow)
library(VIM)
library(tidyverse)
library(magrittr)
library(dplyr)
library(base)
library(ggplot2)
library(glmnet)
library(gridExtra)

#--------------- READING THE DATA -----------------------#

Mortality <- as.data.frame(read.table("Data/7052eng_UntypedDataSet_23032020_165127.csv", header = TRUE, sep = ";"))
Population <- as.data.frame(read.table("Data/7461eng_UntypedDataSet_23032020_165214.csv", header = TRUE, sep = ";"))

#Adjust the names of some variables
colnames(Population)[1] <- "ID"
colnames(Population)[4] <- "Year"
colnames(Population)[5] <- "Population"


colnames(Mortality)[1] <- "ID"
colnames(Mortality)[4] <- "Year"
colnames(Mortality)[5] <- "Total_deaths"

## Keep only the variables of interest
Population <- Population %>%
  select(c(Sex, Age, Year, Population))

Mortality <- Mortality %>%
  select(c(Sex, Age, Year, Total_deaths))

## The variable Year in both datasets needs to be converted to numeric values

# first we have to convert the variable Year from factor to character
Population$Year <- as.character(Population$Year)
Mortality$Year <- as.character(Mortality$Year)

to_number <- function(x) as.numeric(strsplit(x, "J")[[1]][1])

Population$Year <- unlist(lapply(Population$Year, FUN = to_number))
Mortality$Year <- unlist(lapply(Mortality$Year, FUN = to_number))

## Change the encoding of the Sex variable to M, F and T
Population$Sex <- as.factor(Population$Sex)
Mortality$Sex <- as.factor(Mortality$Sex)

levels(Population$Sex) <- c("M", "F", "T")
levels(Mortality$Sex) <- c("M", "F", "T")


## Fix Age variable encodings which are a bit different in the two datasets.
# We need to create a new age group in the Population dataset (51300), 
# containing the 1 to 5 ages like in the Mortality data.
# So, we need to substract the population of newborns from the 0-5 age group.
# We want to keep the newborns as a separate group.


#levels(as.factor(Mortality$Age)) #10010: 0, 51300: 1-5
#levels(as.factor(Population$Age))#10010: 0, 70100: 0-5

Population <- Population %>%
  spread(Age, value = Population) %>%
  mutate(`51300` = `70100` - `10010`) %>%
  gather(key = "Age", value = "Population", `10000`:`51300`)



## Merge the two datasets 
# Now in the Ratio_data we have the 0-5 age group (70100) but also the 0 age group (10010)
Ratio_data <- merge(Population, Mortality, by = c("Sex", "Year", "Age"))


#Plots
p1 <-  ggplot(Population %>%
              filter(Age == 10000), aes(x = Year, y = Population)) +
              geom_point(aes(col=Sex)) +
              geom_line(aes(col=Sex)) +
              labs(title = "Population in the Netherlands", subtitle = "For males and females") +
              theme(legend.position="bottom") +
              xlab("Year")

p2 <- ggplot(Mortality %>%
             filter(Age == 10000), aes(x = Year, y = Total_deaths)) +
             geom_point(aes(col=Sex)) +
             geom_line(aes(col=Sex)) +
             labs(title = "Mortality in the Netherlands", subtitle = "For males and females") +
             theme(legend.position="bottom") +
             xlab("Year")

grid.arrange(p1, p2, ncol = 2)


# No missing values
summary(Ratio_data)

#--------------------- CALCULATING MORTALITY RATIOS ------------------

# Mortality ratios
Ratio_data$Mortality_ratio <- round(Ratio_data$Total_deaths/Ratio_data$Population, digits = 5)

# checking for weird values
summary(Ratio_data$Mortality_ratio)

# log ratios
Ratio_data$log_ratio <- log(Ratio_data$Mortality_ratio)

# Plot sex-specific log ratios 
ggplot(Ratio_data %>%
      filter(Age == 10000, Sex != "T"), aes(x = Year, y = log_ratio, group = Sex)) +
  geom_point(aes(col = Sex)) +
  geom_line(aes(col = Sex)) +
  labs(title = "Log-Mortality rates in the Netherlands", subtitle = "For males and females separately") +
  ylab("log-mortality rate")


# Plot population of two age groups
p3 <- ggplot(Ratio_data %>%
            filter(Age == 22000, Sex == "T" ), aes(x = Year, y = Population, group = 1)) +
  geom_point() +
  geom_line() +
  labs(title = "Population of 95+ year old group") +
  xlab("Year")

p4 <- ggplot(Ratio_data %>%
            filter(Age == 10010, Sex == "T" ), aes(x = Year, y = Population, group = 1)) +
  geom_point() +
  geom_line() +
  labs(title = "Population of <1 year old group") +
  xlab("Year")

grid.arrange(p3, p4, ncol = 2)


#Plot Life expectancies of the two genders
Life_exp <- as.data.frame(read.table("Data/Health_expectancy__since_1981_29032020_194215.csv", header = FALSE, sep = ","))
Life_exp <- as.data.frame(Life_exp[-c(1:5, 20), ])
Life_exp <- Life_exp %>%
  separate(col = 1, into = c("Sex", "Age", "Periods", "Life.expectancy",
                             "Life_exp_in_good_health", "Life_exp_no_physical_limit",
                             "Life_exp_no_chr", "mental_health", "GALI"), sep = ";") %>%
  select(c(Sex, Age, Periods, Life.expectancy))



ggplot(Life_exp,
       aes(x = as.factor(Periods), y = Life.expectancy, fill = Sex)) +
  geom_bar(position="dodge", stat="identity") +
  labs(title = "Life expectancy") +
  ylab("Life expectancy") +
  theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = "bottom") +
  xlab("Year") +
  scale_fill_manual(values=c("red", "orange"))

levels(as.factor(Ratio_data$Age))

#---------------------------- PREPARATION FOR SMOOTHING ---------------------

# First we need to sort our data by age groups.
Ratio_data <- Ratio_data %>%
  filter(Age != 10000)

# Since the age group 22000 is for 95+ year old people, we assign it to the largest number
Ratio_data[Ratio_data$Age == 22000, 3] <- 72000
# Sort by age
Ratio_data <- Ratio_data %>% 
  arrange(Age)

# For fitting least squares we need wide format
Ratio_wide <- Ratio_data %>%
  select(-c(Population, Mortality_ratio, Total_deaths)) %>%
  spread(Year, log_ratio)

### Bsplines smoothing 

#Creating the basis function decomposition using B-splines
# We have 20 age groups
#0 corresponds to newborns (10010) and 100 corresponds 95+ age group
age_points <- seq(from = 0, to = 100, by = 5)


## ---------------- Test the impact of K, at the smoothed curve (only for visualization)

#Now we need to estimate the coefficients in order to produce the smooth function.
#We are doing that using OLS
#Y is the mortality ratios of females of every age, at 2000. We want to fit a smooth curve to its data points.
Y <- Ratio_wide %>% filter(Sex == "F") %>% select(`2000`)
Y <- unlist(Y)


### K = 4
basis1 <- create.bspline.basis(c(0, 100), nbasis = 4)
basis_eval1 <- eval.basis(age_points, basis1)


model1 <- lm(Y ~ 0 + basis_eval1)
Yhat1 <- model1$fitted.values

### K = 
basis2 <- create.bspline.basis(c(0, 100), nbasis = 50)
basis_eval2 <- eval.basis(age_points, basis2)


model2 <- lm(Y ~ 0 + basis_eval2)
Yhat2 <- model2$fitted.values


### Plot
plot(age_points, Y, type="n",lwd=4, col="black",
     xlab="Age", ylab="Mortality ratio", 
     main = "B-splines curves")

points(age_points, Y, pch=1, cex=.5, col="blue", lwd=1)
lines(age_points, Yhat1, lwd=1, col="red")
lines(age_points, Yhat2, lwd=1, col="black")
legend("bottomright", legend = c("K=4", "K=50"), col = c("red", "black"),lty=1:2)


#---------------------------- PENALIZED B-SPLINES SMOOTHING ---------------------

set.seed(1993) 
# Discrete mortality observations for males and females
ym <- as.matrix(Ratio_wide %>% filter(Sex == "M") %>% select(-c(Sex, Age)))
yf <- as.matrix(Ratio_wide %>% filter(Sex == "F") %>% select(-c(Sex, Age)))



### In order to choose the best value for lambda, we compute the RMSE, GCV and degrees of freedom
# for a sequence of values for lambda.
summarise_penalties <- function(loglambdas, basis_expansion, data, argvalues, fdParobj){
  #Initializing an empty dataframe, based on the length of lambdas sequence
  results <- numeric(0)
  
  #Looping for every values of lambda
  for (i in 1:length(loglambdas)){
    loglambda <- loglambdas[i]
    #Fitting the penalized B-splines basis
    fdParobj$lambda <- 10^(loglambda)
    smooth_basis <- smooth.basis(argvals = argvalues, y = data, fdParobj)
    
    #Extracting predicted curves, df, GCV and RMSE
    smooth_curve <- smooth_basis$fd
    df <- round(smooth_basis$df, digits = 5)
    gcv <- round(mean(smooth_basis$gcv), digits = 5)
    error_set <- eval.fd(argvalues, smooth_curve) - data
    RMSE <- mean(apply(error_set, 2, function(x) round(sqrt(mean(x^2)), digits = 3)))
    #Adding a row in the results table
    results <- rbind(results, c(loglambda, df, gcv, RMSE))
  }
  
  results <- as.data.frame(results)
  colnames(results) <- c("log-lambda", "df", "GCV", "RMSE")
  return(results)
}


loglambdas <- seq(-20, -3, 0.25)
basis_exppen <- create.bspline.basis(c(0,100), norder = 9, breaks = c(0,100))
fdParobj <- fdPar(basis_exppen, Lfdobj = NULL, lambda = 10^(-5))
## RESULTS
penalty_summary <- summarise_penalties(loglambdas, basis_expansion = basis_exppen, data = ym, argvalues = age_points, fdParobj)


#Plot the change of GCV value
plot(penalty_summary$`log-lambda`, penalty_summary$GCV, xlab = "Log-lamnda",
     ylab = "GCV value", main = "GCV criterion", type = "b")

#Optimal log lambda is equal to -9. Has the lowest GCV and RMSE values.
lambda <- 10^(-15)
fdParobj <- fdPar(basis_exppen, Lfdobj = NULL, lambda)

#Smooth curves for males
mortality_smooth_m <- smooth.basis(argvals = age_points, y = ym, fdParobj)
mortalityfd_m <- mortality_smooth_m$fd  

#Smooth curves for females
mortality_smooth_f <- smooth.basis(argvals = age_points, y = yf, fdParobj)
mortalityfd_f <- mortality_smooth_f$fd 

## Unsmoothed data
dataM <- fts(age_points, Ratio_wide %>%
               filter(Sex == "M") %>%
               select(-c(Sex, Age)), xname = "Age", yname = "Log-Mortality")


dataF <- fts(age_points, Ratio_wide %>%
               filter(Sex == "F") %>%
               select(-c(Sex, Age)), xname = "Age", yname = "Log-Mortality")


### Plots for males
par(mfrow = c(1,2))
plot(dataM, plot.type = "functions", 
     plotlegend = TRUE, legendpos = "bottomright", 
     main = "Original data for males")

plot.fd(mortalityfd_m, col = rainbow(80), 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Smoothed data for males")
legend("bottomright", legend = c("1950", "1984", "2018"), col = c("red", "green", "purple"), lty = 1)



### Plots for females
par(mfrow = c(1,2))
plot(dataF, plot.type = "functions", 
     plotlegend = TRUE, legendpos = "bottomright", 
     main = "Original data for females")

plot.fd(mortalityfd_f, col = rainbow(80), 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Smoothed data for females")
legend("bottomright", legend = c("1950", "1984", "2018"), col = c("red", "green", "purple"), lty = 1)

### Functional HDR for males
#functional
#fboxplot(dataM, plot.type = "functional", type = "hdr", projmethod = "PCAproj", legendpos = "bottomright")
#bivariate
#fboxplot(dataM, plot.type = "bivariate", type = "hdr", projmethod = "PCAproj")


### Functional HDR for females
#functional
#fboxplot(dataF, plot.type = "functional", type = "hdr", projmethod = "PCAproj")
#bivariate
#fboxplot(dataF, plot.type = "bivariate", type = "hdr", projmethod = "PCAproj")

# Smooth functional time series
smoothM <- fts(age_points, eval.fd(age_points, mortalityfd_m), xname = "Age", yname = "Log-Mortality")
smoothF <- fts(age_points, eval.fd(age_points, mortalityfd_f), xname = "Age", yname = "Log-Mortality")



###----------------------- FPCA MODELS ON MALE POPULATION-----------------------

##First two principal components for visualization
plot(forecast(ftsm(smoothM, order = 2), h = 20), "components")

###Training set
y_train_m <- fts(age_points, as.data.frame(eval.fd(age_points, mortalityfd_m)) %>%
                   select(`1950`:`1999`), xname = "Age", yname = "Log-Mortality ratios")

###Forecasts for different values of K
fpca_m_4 <- forecast(ftsm(y_train_m, order = 4), h = 19)
fpca_m_6 <- forecast(ftsm(y_train_m, order = 6), h = 19)
fpca_m_8 <- forecast(ftsm(y_train_m, order = 8), h = 19)

###Goodness of fit 
#For each choice of K we get slightly different errors.
summary(ftsm(y_train_m, order = 4))
summary(ftsm(y_train_m, order = 6))
summary(ftsm(y_train_m, order = 8))

###RESIDUALS
plot.fmres(residuals.fm(ftsm(y_train_m, order = 6)),
           type = "fts", xlab = "Age", main = "Residuals for males (FPCA)")

#Smooth training data
train_smooth_m <- smooth.basis(argvals = age_points,
                               y = as.matrix(Ratio_wide %>%
                                               filter(Sex == "M") %>%
                                               select(-c(Sex, Age)) %>%
                                               select(`1950`:`1999`)),
                               fdParobj = fdParobj)$fd
#Smooth forecasts
fpca_forecasts_m_6 <- smooth.basis(argvals = age_points, y = fpca_m_6$mean$y, 
                                   fdParobj)$fd

### Plot data with forecasts
#Smooth data
plot.fd(train_smooth_m, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for males")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
#Smooth forecasts
plot.fd(fpca_forecasts_m_6, col = rainbow(19), add = TRUE)





###--------------------------- Weighted FPCA models -----------------------------

###Forecasts for males
wfpca_m_4 <- forecast(ftsm(y_train_m, order = 4, weight = TRUE), h = 19)
wfpca_m_6 <- forecast(ftsm(y_train_m, order = 6, weight = TRUE), h = 19)
wfpca_m_8 <- forecast(ftsm(y_train_m, order = 8, weight = TRUE), h = 19)

###Goodness of fit 
summary(ftsm(y_train_m, order = 4, weight = TRUE))
summary(ftsm(y_train_m, order = 6, weight = TRUE))
summary(ftsm(y_train_m, order = 8, weight = TRUE))

### RESIDUALS
plot.fmres(residuals.fm(ftsm(y_train_m, order = 6, weight = TRUE)),
           type = "fts", xlab = "Age", main = "Residuals for males (Weighted FPCA)")


#The difference between the weighted FPCA and the regular FPCA can be clearly seen from the residuals
#The weighted FPCA focuses on most recent data which have low residuals, and not so much on older data
#which have clearly larger residuals.



#Smoothing the forecasts
wfpca_forecasts_m_6 <- smooth.basis(argvals = age_points, y = wfpca_m_6$mean$y, 
                                    fdParobj = fdParobj)$fd


### Plot data with forecasts
plot.fd(train_smooth_m, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for males")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(wfpca_forecasts_m_6, col = rainbow(19), add = TRUE)



###----------------------- FPCA MODELS ON FEMALE POPULATION-----------------------


##First two principal components for visualization
plot(forecast(ftsm(smoothF, order = 2), h = 20), "components")

###Training set
y_train_f <- fts(age_points, as.data.frame(eval.fd(age_points, mortalityfd_f)) %>%
                   select(`1950`:`1999`), xname = "Age", yname = "Log-Mortality ratios")

###Forecasts for different values of K
fpca_f_4 <- forecast(ftsm(y_train_f, order = 4), h = 19)
fpca_f_6 <- forecast(ftsm(y_train_f, order = 6), h = 19)
fpca_f_8 <- forecast(ftsm(y_train_f, order = 8), h = 19)

###Goodness of fit 
#For each choice of K we get slightly different errors.
summary(ftsm(y_train_f, order = 4))
summary(ftsm(y_train_f, order = 6))
summary(ftsm(y_train_f, order = 8))

### RESIDUALS
plot.fmres(residuals.fm(ftsm(y_train_f, order = 6)),
           type = "fts", xlab = "Age", main = "Residuals for females (FPCA)")


#Smoothing the forecasts
fpca_forecasts_f_6 <- smooth.basis(argvals = age_points, y = fpca_f_6$mean$y, 
                                   fdParobj = fdParobj)$fd

#Smooth training data
train_smooth_f <- smooth.basis(argvals = age_points,
                               y = as.matrix(Ratio_wide %>%
                                               filter(Sex == "F") %>%
                                               select(-c(Sex, Age)) %>%
                                               select(`1950`:`1999`)),
                               fdParobj = fdParobj)$fd

### Plot data with forecasts
plot.fd(train_smooth_f, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for females")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(fpca_forecasts_f_6, col = rainbow(19), add = TRUE)





###--------------------------- Weighted FPCA models -----------------------------

###Forecasts for females
wfpca_f_4 <- forecast(ftsm(y_train_f, order = 4, weight = TRUE), h = 19)
wfpca_f_6 <- forecast(ftsm(y_train_f, order = 6, weight = TRUE), h = 19)
wfpca_f_8 <- forecast(ftsm(y_train_f, order = 8, weight = TRUE), h = 19)




###Goodness of fit 
summary(ftsm(y_train_f, order = 4, weight = TRUE))
summary(ftsm(y_train_f, order = 6, weight = TRUE))
summary(ftsm(y_train_f, order = 8, weight = TRUE))



### RESIDUALS
plot.fmres(residuals.fm(ftsm(y_train_f, order = 6, weight = TRUE)),
           type = "fts", xlab = "Age", main = "Residuals for females (Weighted FPCA)")

#Smoothing the forecasts
wfpca_forecasts_f_6 <- smooth.basis(argvals = age_points, y = wfpca_f_6$mean$y, 
                                    fdParobj = fdParobj)$fd

### Plot data with forecasts
plot.fd(train_smooth_f, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for females")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(wfpca_forecasts_f_6, col = rainbow(19), add = TRUE)



###------------------Functional Partial Least Squares on male population ------------------------

##Forecasts for different values of K
fplsr_m_4 <- forecastfplsr(y_train_m, components = 4, h = 19)
fplsr_m_6 <- forecastfplsr(y_train_m, components = 6, h = 19)
fplsr_m_8 <- forecastfplsr(y_train_m, components = 8, h = 19)
fplsr_m_10 <- forecastfplsr(y_train_m, components = 10, h = 19)



###Goodness of fit 
#For each choice of K we get slightly different errors.
summary(fplsr(y_train_m, order = 4))
summary(fplsr(y_train_m, order = 6))
summary(fplsr(y_train_m, order = 8))
summary(fplsr(y_train_m, order = 10))

### RESIDUALS
#Plot smooth residuals (K = 6)

fplsr_res_m <- fts(age_points, t(residuals.fm(fplsr(y_train_m, order = 6))$z), xname = "Age", yname = "Residuals")

plot(fplsr_res_m , plot.type = "functions", 
     plotlegend = TRUE, legendpos = "bottomright", 
     main = "Residuals for males (FPLSR)")


# Smooth forecasts
fplsr_forecasts_m_6 <- smooth.basis(argvals = age_points, y = fplsr_m_6$y, 
                                    fdParobj = fdParobj)$fd

#Plot data with forecasts for K=6
plot.fd(train_smooth_m, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for males")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(fplsr_forecasts_m_6, col = rainbow(19), add = TRUE)



###------------------Functional Partial Least Squares on female population ------------------------

##Forecasts for different values of K
fplsr_f_4 <- forecastfplsr(y_train_f, components = 4, h = 19)
fplsr_f_6 <- forecastfplsr(y_train_f, components = 6, h = 19)
fplsr_f_8 <- forecastfplsr(y_train_f, components = 8, h = 19)

###Goodness of fit 
#For each choice of K we get slightly different errors.
summary(fplsr(y_train_f, order = 4))
summary(fplsr(y_train_f, order = 6))
summary(fplsr(y_train_f, order = 8))

### RESIDUALS
#Plot smooth residuals (K = 6)
fplsr_res_f <- fts(age_points, t(residuals.fm(fplsr(y_train_f, order = 6))$z), xname = "Age", yname = "Residuals")

plot(fplsr_res_f , plot.type = "functions", 
     plotlegend = TRUE, legendpos = "bottomright", 
     main = "Residuals for females (FPLSR)")

# Smooth forecasts
fplsr_forecasts_f_6 <- smooth.basis(argvals = age_points, y = fplsr_f_6$y, 
                                    fdParobj = fdParobj)$fd

#Plot data with forecasts for K=6
plot.fd(train_smooth_f, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for females")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(fplsr_forecasts_f_6, col = rainbow(19), add = TRUE)



### --------------------- Produce Plots ------------------------
###RESIDUALS FPCA
par(mfrow = c(1,2))

plot.fmres(residuals.fm(ftsm(y_train_m, order = 6)),
           type = "fts", xlab = "Age", ylim = c(-0.15, 0.15), main = "Residuals for males")

plot.fmres(residuals.fm(ftsm(y_train_f, order = 6)),
           type = "fts", xlab = "Age", ylim = c(-0.15, 0.15), main = "Residuals for females")


par(mfrow = c(1,2))
### Plot data with forecasts FPCA
#Smooth data
plot.fd(train_smooth_m, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for males")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
#Smooth forecasts
plot.fd(fpca_forecasts_m_6, col = rainbow(19), add = TRUE)

plot.fd(train_smooth_f, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for females")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(fpca_forecasts_f_6, col = rainbow(19), add = TRUE)


par(mfrow = c(1,2))  
### RESIDUALS WFPCA
plot.fmres(residuals.fm(ftsm(y_train_m, order = 6, weight = TRUE)),
           type = "fts", xlab = "Age", ylim = c(-0.2, 0.2), main = "Residuals for males")

plot.fmres(residuals.fm(ftsm(y_train_f, order = 6, weight = TRUE)),
           type = "fts", xlab = "Age", ylim = c(-0.2, 0.2), main = "Residuals for females")


par(mfrow = c(1,2))
### Plot data with forecasts WFPCA
plot.fd(train_smooth_m, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for males")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(wfpca_forecasts_m_6, col = rainbow(19), add = TRUE)


plot.fd(train_smooth_f, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for females")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(wfpca_forecasts_f_6, col = rainbow(19), add = TRUE)



# Residuals FPLSR
par(mfrow = c(1,2))
plot(fplsr_res_m , plot.type = "functions", legendpos = "bottomright", 
     main = "Residuals for males")

plot(fplsr_res_f , plot.type = "functions", legendpos = "bottomright", 
     main = "Residuals for females")


# FORECASTS FPLSR
par(mfrow = c(1,2))
plot.fd(train_smooth_m, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for males")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(fplsr_forecasts_m_6, col = rainbow(19), add = TRUE)

plot.fd(train_smooth_f, col = "grey", 
        xlab = "Age", ylab = "Log-Mortality", 
        main = "Forecasts for females")
legend("topleft", legend = c("2000", "2018"), col = c("red", "purple"), lty = 1:2)
plot.fd(fplsr_forecasts_f_6, col = rainbow(19), add = TRUE)


#---------------------- MODEL EVALUATION FOR MALES -------------------------

### First we are going to create a function that takes the forecasts and test data as input.
## And then computes some errors
##Test set 
y_test_m <- as.data.frame(eval.fd(age_points, mortalityfd_m)) %>% select(`2000`:`2018`)             


### Compute the RMSE for each curve (2000-2018)
## For every model

# FPCA models
fpca_m_4_rmse <- apply(fpca_m_4$mean$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fpca_m_6_rmse <- apply(fpca_m_6$mean$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fpca_m_8_rmse <- apply(fpca_m_8$mean$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))

# Weighted FPCA models
wfpca_m_4_rmse <- apply(wfpca_m_4$mean$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
wfpca_m_6_rmse <- apply(wfpca_m_6$mean$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
wfpca_m_8_rmse <- apply(wfpca_m_8$mean$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))

# FPLSR models
fplsr_m_4_rmse <- apply(fplsr_m_4$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fplsr_m_6_rmse <- apply(fplsr_m_6$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fplsr_m_8_rmse <- apply(fplsr_m_8$y - y_test_m, 2, function(x) round(sqrt(mean(x^2)), digits = 3))


eval_males <- as.data.frame(rbind(fpca_m_6_rmse, wfpca_m_6_rmse, fplsr_m_6_rmse))

eval_males$model <- c("FPCA", "WFPCA", "FPLSR")

apply(eval_males[ ,-20], 1, FUN = mean)


#---------------------- MODEL EVALUATION FOR FEMALES -------------------------

### First we are going to create a function that takes the forecasts and test data as input.
## And then computes some errors
##Test set 
y_test_f <- as.data.frame(eval.fd(age_points, mortalityfd_f)) %>% select(`2000`:`2018`) 

### Compute the RMSE for each curve (2000-2018)
## For every model
fpca_f_4_rmse <- apply(fpca_f_4$mean$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fpca_f_6_rmse <- apply(fpca_f_6$mean$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fpca_f_8_rmse <- apply(fpca_f_8$mean$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))

wfpca_f_4_rmse <- apply(wfpca_f_4$mean$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
wfpca_f_6_rmse <- apply(wfpca_f_6$mean$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
wfpca_f_8_rmse <- apply(wfpca_f_8$mean$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))

# FPLSR models
fplsr_f_4_rmse <- apply(fplsr_f_4$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fplsr_f_6_rmse <- apply(fplsr_f_6$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))
fplsr_f_8_rmse <- apply(fplsr_f_8$y - y_test_f, 2, function(x) round(sqrt(mean(x^2)), digits = 3))

eval_females <- as.data.frame(rbind(fpca_f_6_rmse, wfpca_f_6_rmse, fplsr_f_6_rmse))

eval_females$model <- c("FPCA", "WFPCA", "FPLSR")


apply(eval_females[ ,-20], 1, FUN = mean)

p5 <- ggplot( eval_males %>%
                gather(key = "Year", value = "RMSE", `2000`:`2018`),
              aes(x = as.numeric(Year), y = RMSE, group = model)) +
  geom_line(aes(col = model)) +
  geom_point(aes(col = model)) +
  labs(title = "Model evaluation for males") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
  xlab("Year") +
  theme(legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 6))


# FEMALES
p6 <- ggplot( eval_females %>%
                gather(key = "Year", value = "RMSE", `2000`:`2018`),
              aes(x = as.numeric(Year), y = RMSE, group = model)) +
  geom_line(aes(col = model)) +
  geom_point(aes(col = model)) +
  labs(title = "Model evaluation for females") +
  ylab("RMSE") +
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
  xlab("Year") +
  theme(legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 6))


grid.arrange(p5, p6, ncol = 2)



### Compute the RMSE for each age group
## For every model


### MALE POPULATION

# FPCA models
fpca_m_6_rmse2 <- apply(fpca_m_6$mean$y - y_test_m, 1, function(x) round(sqrt(mean(x^2)), digits = 3))

# Weighted FPCA models
wfpca_m_6_rmse2 <- apply(wfpca_m_6$mean$y - y_test_m, 1, function(x) round(sqrt(mean(x^2)), digits = 3))

# FPLSR models
fplsr_m_6_rmse2 <- apply(fplsr_m_6$y - y_test_m, 1, function(x) round(sqrt(mean(x^2)), digits = 3))

eval_males2 <- as.data.frame(rbind(fpca_m_6_rmse2, wfpca_m_6_rmse2, fplsr_m_6_rmse2))
colnames(eval_males2) <- as.character(age_points)
eval_males2$model <- c("FPCA", "WFPCA", "FPLSR")

apply(eval_males2[ ,-22], 1, FUN = mean)

### FEMALE POPULATION

# FPCA models
fpca_f_6_rmse2 <- apply(fpca_f_6$mean$y - y_test_f, 1, function(x) round(sqrt(mean(x^2)), digits = 3))

# Weighted FPCA models
wfpca_f_6_rmse2 <- apply(wfpca_f_6$mean$y - y_test_f, 1, function(x) round(sqrt(mean(x^2)), digits = 3))

# FPLSR models
fplsr_f_6_rmse2 <- apply(fplsr_f_6$y - y_test_f, 1, function(x) round(sqrt(mean(x^2)), digits = 3))

eval_females2 <- as.data.frame(rbind(fpca_f_6_rmse2, wfpca_f_6_rmse2, fplsr_f_6_rmse2))
colnames(eval_females2) <- as.character(age_points)
eval_females2$model <- c("FPCA", "WFPCA", "FPLSR")

apply(eval_females2[ ,-22], 1, FUN = mean)


# MALES
p7 <- ggplot( eval_males2 %>%
                gather(key = "Age", value = "RMSE", `0`:`100`),
              aes(x = as.numeric(Age), y = RMSE, group = model)) +
  geom_line(aes(col = model)) +
  geom_point(aes(col = model)) +
  labs(title = "Model evaluation for males") +
  ylab("RMSE") +
  xlab("Age") +
  theme(legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 6))


# FEMALES
p8 <- ggplot( eval_females2 %>%
                gather(key = "Age", value = "RMSE", `0`:`100`),
              aes(x = as.numeric(Age), y = RMSE, group = model)) +
  geom_line(aes(col = model)) +
  geom_point(aes(col = model)) +
  labs(title = "Model evaluation for females") +
  ylab("RMSE") +
  xlab("Age") +
  theme(legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 6))


grid.arrange(p7, p8, ncol = 2)