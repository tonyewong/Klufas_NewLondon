# Alex Klufas
# optimization_sf.R
# Written on June 14, 2017 
# Modified on June 14, 2017
#
# Script in order to determine best fit parameters using optimization 

#---------OPTIMIZATION---------------------------------------------------------------------------------------- 
#library(DEoptim)

log.like.calc <- function(data, location1, scale1, shape1){
  ll<-sum(devd(data, loc=location1, scale=scale1, shape=shape1, type=c('GEV'), log=TRUE))
  return(ll)
}

neg.log.like.calc <- function(p, data){
  mu <- p[1]
  sigma <- p[2]
  xi <- p[3]
  
  nll <- -1*sum(devd(data, loc=mu, scale=sigma, shape=xi, type=c('GEV'), log=TRUE))
  
  return(nll)
}


upper.bound <- c(3000,3000,5)
lower.bound <- c(0,0,-5)

p.names<- c('mu', 'sigma', 'xi')

de.optim.val <- DEoptim.control(VTR = -Inf, strategy = 2, bs = FALSE, NP = 200,
                                itermax = 200, CR = 0.5, F = 0.8, trace = TRUE, initialpop = NULL,
                                storepopfrom = 101, storepopfreq = 1, p = 0.2, c = 0,
                                parallelType = 0, cluster = NULL, packages = c(), parVar = c(),
                                foreachArgs = list())

optim.like <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data= lsl.max)

value.to.use <- pevd(x.hgt, loc = optim.like$optim$bestmem[1], scale=optim.like$optim$bestmem[2], shape=optim.like$optim$bestmem[3], type=c('GEV'), lower.tail = FALSE)

#plotting sf's (both LHS and DEoptim)
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT')

points(lsl.sorted.vals, log10(esf.vals)) 
lines(x.hgt, log10(sf.hgt2), col='red')

for (s in 1:10){
  sf.hgt.test <- 1-pevd(x.hgt, loc=sf.location.top10[s], scale=sf.scale.top10[s], shape=sf.shapes.top10[s], type=c("GEV"), lower.tail=TRUE)
  lines(x.hgt, log10(sf.hgt.test), col='blue')
}

lines(x.hgt, log10(value.to.use), col='green')

#-------------non-stationary work 

#read temp data 
temp.data <- read.table('noaa-temp-1880-2017.csv', header = TRUE, sep=',')
temp.years <- temp.data$Year
temp.values <- temp.data$Value[59:134]     #limit temp data to just the years I am working with 

#log likelihood
log.like.calc <- function(data, location1, scale1, shape1){
  ll<-sum(devd(data, loc=location1, scale=scale1, shape=shape1, type=c('GEV'), log=TRUE))
  return(ll)
}

#negative log likelihood calculator 
neg.log.like.calc <- function(p, parnames, data, temps){
  n.parnames <- length(parnames)
  #checks length of parameters, which will be used to determine which parameters are non stationary 
  
  #if all parameters are stationary 
  if (n.parnames == 3){
    mu <- p[1]
    sigma <- p[2]
    xi <- p[3]
  }
  
  #one stationary parameter
  else if(n.parnames ==4){
    if (parnames[1] == 'mu0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi <- p[4]
      
      mu <- mu0 + mu1*temps
      sigma <- exp(sigma)
    }
    else if (parnames[2] == 'sigma0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi <- p[4]
      
      sigma <- exp(sigma0+ sigma1*temps)
    }
    else{
      mu <- p[1]
      sigma <- p[2]
      xi0 <- p[3]
      xi1 <- p[4]
      
      sigma <- exp(sigma)
      xi <- xi0 + xi1*temps
    }
  }
  
  #two stationary parameters
  else if (n.parnames == 5){
    if (parnames[1] == 'mu0' & parnames[3]=='sigma0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma0 <- p[3]
      sigma1 <- p[4]
      xi <- p[5]
      
      mu <- mu0 + mu1*temps
      sigma <- exp(sigma0 + sigma1*temps)
    }
    else if(parnames[2] == 'sigma0' & parnames[4] == 'xi0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      
      sigma <- exp(sigma0 + sigma1*temps)
      xi <- xi0 + xi1*temps
      
    }
    else{
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      
      mu <- mu0 + mu1*temps
      sigma <- exp(sigma)
      xi <- xi0 + xi1*temps
    }
  }
  #all parameters non stationary 
  else if (n.parnames ==6){
    mu0 <- p[1]
    mu1 <- p[2]
    sigma0 <- p[3]
    sigma1 <- p[4]
    xi0 <- p[5]
    xi1 <- p[6]
    
    mu <- mu0 + mu1*temps
    sigma <- exp(sigma0 + sigma1*temps)
    xi <- xi0 + xi1*temps
    
  }
  nll <- -1*sum(devd(data, loc=mu, scale=sigma, shape=xi, type=c('GEV'), log=TRUE))
  
  return(nll)
}

#---------Non Stationary Plots Using DE Optim-------------------------------------------------------------------

#----non stationary mu only---------------------------------works!!!!!----------------------
p.names <- c('mu0', 'mu1', 'sigma', 'xi')

upper.bound <- c(3000,1000, 1000, 5)
lower.bound <- c(0,0, -100, -5)

optim.like.temp.mu <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)

plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT, Non Stationary Mu')
points(lsl.sorted.vals, log10(esf.vals))

fit.vals.optim.mu <- rep(0, length(lsl.max))

for (j in 1:length(temp.values)){
  sf.optim.mu <- 1-pevd(lsl.max[j], loc = (optim.like.temp.mu$optim$bestmem[1] + optim.like.temp.mu$optim$bestmem[2] * temp.values[j]), scale = exp(optim.like.temp.mu$optim$bestmem[3]) , shape = optim.like.temp.mu$optim$bestmem[4],type=c('GEV'))
  fit.vals.optim.mu[j] <- sf.optim.mu
  points(lsl.max[j], log10(sf.optim.mu), pch = 2, col='red')
}
#-----non sationary mu and sigma ----------- works ---------------------------------------- --------------------------------------------------- --------------------------------------------------- 

p.names <- c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi')

upper.bound <- c(3000,100, 1000, 10 , 5)
lower.bound <- c(0,-100, 0, 0, -5)

optim.like.temp.mu.sigma <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values,parnames = p.names)

fit.vals.optim.mu.sigma <- rep(0, length(lsl.max))

#-------run for plot of non sationary mu and sigma 
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT, non stationary mu and sigma')
points(lsl.sorted.vals, log10(esf.vals))
for (j in 1:length(temp.values)){
  sf.optim.mu.sigma <- 1 - pevd(lsl.max[j], loc = (optim.like.temp.mu.sigma$optim$bestmem[1] + optim.like.temp.mu.sigma$optim$bestmem[2] * temp.values[j]), scale = exp(optim.like.temp.mu.sigma$optim$bestmem[3] + optim.like.temp.mu.sigma$optim$bestmem[4]*temp.values[j]), shape = optim.like.temp.mu.sigma$optim$bestmem[5],type=c('GEV'))
  fit.vals.optim.mu.sigma[j] <- sf.optim.mu.sigma
  points(lsl.max[j], log10(sf.optim.mu.sigma), pch = 2, col='green')
  
}
#------all variables non sationary --------------------------------------this one does work ------------- --------------------------------------------------- --------------------------------------------------- 
p.names <- c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi0', 'xi1')

upper.bound <- c(3000,100, 1000, 10 , 1, 1)
lower.bound <- c(0,-100, 0, 0, -1, -1)

optim.like.temp.mu.sigma.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)

fit.vals.optim.mu.sigma.xi <- rep(0, length(lsl.max))

#---------run for plot of non sationary mu and sigma and xi
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT, mu, sigma, xi all nonstationary')
points(lsl.sorted.vals, log10(esf.vals))
for (j in 1:length(temp.values)){
  sf.optim.mu.sigma.xi <-1-pevd(lsl.max[j], loc = optim.like.temp.mu.sigma.xi$optim$bestmem[1] + optim.like.temp.mu.sigma.xi$optim$bestmem[2] * temp.values[j], scale = exp(optim.like.temp.mu.sigma.xi$optim$bestmem[3] + optim.like.temp.mu.sigma.xi$optim$bestmem[4]*temp.values[j]), shape = optim.like.temp.mu.sigma.xi$optim$bestmem[5] + optim.like.temp3$optim$bestmem[6]*temp.values[j],type=c('GEV'))
  fit.vals.optim.mu.sigma.xi[j] <- sf.optim.mu.sigma.xi
  points(lsl.max[j], log10(sf.optim.mu.sigma.xi), pch = 2, col='red')
}



#------sigma and xi non stationary----------------------works now! ----------------------- --------------------------------------------------- --------------------------------------------------- 
p.names <- c('mu', 'sigma0', 'sigma1' , 'xi0', 'xi1')

upper.bound <- c(3000, 1000, 10 , 1, 1)
lower.bound <- c(0, -100, 0, -1, -1)

optim.like.temp.sigma.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)

fit.vals.optim.sigma.xi <- rep(0, length(lsl.max))


#---------run for plot of non sationary sigma and xi
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT, non stationary sigma and xi ')
points(lsl.sorted.vals, log10(esf.vals))
for (j in 1:length(temp.values)){
  sf.optim5.sigma.xi <-1-pevd(lsl.max[j], loc = optim.like.temp.sigma.xi$optim$bestmem[1], scale = exp(optim.like.temp.sigma.xi$optim$bestmem[2] + optim.like.temp.sigma.xi$optim$bestmem[3]*temp.values[j]), shape = optim.like.temp.sigma.xi$optim$bestmem[4] + optim.like.temp.sigma.xi$optim$bestmem[5]*temp.values[j],type=c('GEV'))
  fit.vals.optim.sigma.xi[j] <- sf.optim5.sigma.xi
  points(lsl.max[j], log10(sf.optim5.sigma.xi), pch = 2, col='red')
}



#------non stationary mu and xi----------this one works----------------------------------------- --------------------------------------------------- --------------------------------------------------- 
p.names <- c('mu0', 'mu1', 'sigma', 'xi0', 'xi1')

upper.bound <- c(3000, 100, 1000, 10, 10)
lower.bound <- c(0, -100, 0, -1, -10)

optim.like.temp.mu.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)

fit.vals.optim.mu.xi <- rep(0, length(lsl.max))
#---------run for plot of non sationary mu and xi
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT, nonstationary mu and xi')
points(lsl.sorted.vals, log10(esf.vals))
for (j in 1:length(temp.values)){
  sf.optim.mu.xi<-1-pevd(lsl.max[j], loc = optim.like.temp.mu.xi5$optim$bestmem[1] +optim.like.temp.mu.xi$optim$bestmem[2]*temp.values[j], scale = exp(optim.like.temp.mu.xi$optim$bestmem[3]), shape = optim.like.temp.mu.xi$optim$bestmem[4] + optim.like.temp.mu.xi$optim$bestmem[5]*temp.values[j],type=c('GEV'))
  fit.vals.optim.mu.xi[j] <- sf.optim.mu.xi
  points(lsl.max[j], log10(sf.optim.mu.xi), pch = 2, col='red')
}

#------only sigma non stationary -------------------------this works!!-------------------------- --------------------------------------------------- --------------------------------------------------- 

p.names <- c('mu', 'sigma0', 'sigma1', 'xi')

upper.bound <- c(3000, 1000, 10, 1)
lower.bound <- c(0, 0, -10, -1)

optim.like.temp.sigma <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)

fit.vals.optim.sigma <- rep(0, length(lsl.max))
#---------run for plot of non sationary sigma 
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT, sigma non stationary')
points(lsl.sorted.vals, log10(esf.vals))
for (j in 1:length(temp.values)){
  sf.optim.sigma <-1-pevd(lsl.max[j], loc = optim.like.temp.sigma$optim$bestmem[1], scale = exp(optim.like.temp.sigma$optim$bestmem[2] + optim.like.temp.sigma$optim$bestmem[3]*temp.values[j]), shape = optim.like.temp.sigma$optim$bestmem[4],type=c('GEV'))
  fit.vals.optim.sigma[j] <- sf.optim.sigma
  points(lsl.max[j], log10(sf.optim.sigma), pch = 2, col='red')
}

#------ only xi non stationary---------------------------------this works! ------------------ --------------------------------------------------- ---------------------------------------------------  
p.names <- c('mu', 'sigma', 'xi0', 'xi1')

upper.bound <- c(3000, 500 , 1, 1)
lower.bound <- c(0, 0, -1, -1)

optim.like.temp.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)

fit.vals.optim.xi <- rep(0, length(lsl.max))
#---------run for plot of non sationary xi
plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
title('Survival Function of height of sea [mm], New London CT')
points(lsl.sorted.vals, log10(esf.vals))
for (j in 1:length(temp.values)){
  sf.optim.xi <-1-pevd(lsl.max[j], loc = optim.like.temp.xi$optim$bestmem[1], scale = exp(optim.like.temp.xi$optim$bestmem[2]), shape = optim.like.temp.xi$optim$bestmem[3] + optim.like.temp.xi$optim$bestmem[4]*temp.values[j],type=c('GEV'))
  fit.vals.optim.xi[j] <- sf.optim.xi
  points(lsl.max[j], log10(sf.optim.xi), pch = 2, col='red')
}


#--------------functions for calculating fit values using RSME, BIC, AIC-------------------#
rmse.calc <- function(data, fit.vals, esf.vals){
  val <- sqrt(1/(length(data))*(sum(esf.vals - fit.vals)^2))
  return(val)
}

aic.calc <- function(n.param, optim.best.val){
  val <- 2*n.param - 2*log(optim.best.val)
  return(val)
}

bic.calc <- function(n.param, optim.best.val, data){
  val <- n.param*log(length(data)) - 2*log(optim.best.val)
  return(val)
}
#------------------------------------------------------------------------------------------#

#--------------------------RMSE For Optim Fits------------------------------------------------------------------------------------------------

#only mu non stationary 
rmse.mu <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.mu)^2))

#only sigma non stationary 
rmse.sigma<- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.sigma)^2))

#xi non stationary 
rmse.xi <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.xi))^2)

#mu and sigma non stationary 
rmse.mu.sigma <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.mu.sigma)^2))

#mu and xi non stationary 
rmse.mu.xi<- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.mu.xi))^2)

#sigma and xi non stationary 
rmse.sigma.xi <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.sigma.xi))^2)

#mu, sigma, xi non stationary 
rmse.mu.sigma.xi <- sqrt(1/(length(lsl.max))*(sum(esf.vals - sf.optim.mu.sigma.xi))^2)

#----------------AIC For Optim Fits------------------------------------------------------------------------------------------------------------------------
#2k - 2 ln (L)

aic.calc <- function(n.param, optim.best.val){
  val <- 2*n.param - 2*log(optim.best.val)
  return(val)
}

aic.mu <- 2*4 - 2*log(optim.like.temp.mu$optim$bestval)

aic.sigma <- 2*4 - 2*log(optim.like.temp.sigma$optim$bestval)

aic.xi <- 2*4 - 2*log(optim.like.temp.xi$optim$bestval)

aic.mu.sigma <- 2*5 - 2*log(optim.like.temp.mu.sigma$optim$bestval)

aic.mu.xi<- 2*5 - 2*log(optim.like.temp.mu.xi$optim$bestval)

aic.sigma.xi<- 2*5 - 2*log(optim.like.temp.sigma.xi$optim$bestval)

aic.mu.sigma.xi <- 2*6 - 2*log(optim.like.temp.mu.sigma.xi$optim$bestval)

#------------------BIC For Optim Fits------------------------------------------------------------------------------------------------------------------------ 
#kln(n) - 2 ln(L)

bic.calc <- function(n.param, optim.best.val, data){
  val <- n.param*log(length(data)) - 2*log(optim.best.val)
  return(val)
}

bic.mu<- 4*log(length(lsl.max)) - 2*log(optim.like.temp.mu$optim$bestval)

bic.sigma <- 4*log(length(lsl.max))- 2*log(optim.like.temp.sigma$optim$bestval)

bic.xi <- 4*log(length(lsl.max))- 2*log(optim.like.temp.xi$optim$bestval)

bic.mu.sigma<- 5*log(length(lsl.max))- 2*log(optim.like.temp.mu.sigma$optim$bestval)

bic.mu.xi <- 5*log(length(lsl.max))- 2*log(optim.like.temp.mu.xi$optim$bestval)

bic.sigma.xi <- 5*log(length(lsl.max))- 2*log(optim.like.temp.sigma.xi$optim$bestval)

bic.mu.sigma.xi <- 6*log(length(lsl.max))- 2*log(optim.like.temp.mu.sigma.xi$optim$bestval)