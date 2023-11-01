## Step 1 Prob distribution for disease onset to hospital #####
## Disease onset to hospitalization; calculating shape and scale from mean and IQR
f = function(par, targetmedian, targetIQR) {
  k = par[1]
  theta = par[2]
  
  #median
  median = qgamma(0.5,shape = k,scale = theta)
  #IQR 
  IQR = qgamma(0.75,shape = k,scale = theta) - qgamma(0.25,shape = k,scale = theta)
  
  error = (median-targetmedian)^2+(IQR-targetIQR)^2
  return(error)
}

### optim function below gives result k = shape theta = scale
p2delta = optim(par  = c(1,1), f, targetmedian = 5.1, targetIQR = 2.4)$par
## Delta = Shape= 8.51, scale= 0.62

##optim function below gives result k = shape theta = scale
p2omi = optim(par  = c(1,1), f, targetmedian = 2.8, targetIQR = 1.6)$par
## Omicron = Shape= 5.873, scale= 0.505

###Step 2: Prob distribution for Incubation period ###
## Incubation period (assuming normal distribution)

meandel = 4.41
sedel = (5.50-3.76)/3.92

meanomi = 3.42
seomi = (3.69-2.88)/3.92

###Step 3 Convulation of 2 probabilities ####

#generate 1000 random values that uses
# a gamma distribution shape
#install.packages("survival")
install.packages("fitdistrplus")
library(fitdistrplus)

##Delta variant
datadelta = rgamma(10000, p2delta[1], scale = p2delta[2]) + rnorm(10000, meandel, sedel)
#fit the dataset to a gamma dist using mle
ansdelta = fitdistrplus::fitdist(datadelta, distr = "gamma", method = "mle")
summary(ansdelta)

tiff("distribution-delta.tiff", res = 300, units = "in", width = 8, height = 7, compression = "lzw")
#Display the summary of ans
plot(ansdelta)
dev.off()


##Omicron variant
dataom = rgamma(1000, p2omi[1], scale = p2omi[2]) + rnorm(1000, meanomi, seomi)
ansomi = fitdistrplus::fitdist(dataom, distr = "gamma", method = "mle")
summary(ansomi)

tiff("distribution-omicron.tiff", res = 300, units = "in", width = 8, height = 7, compression = "lzw")
plot(ansomi)
dev.off()

