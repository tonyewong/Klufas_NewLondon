# Alex Klufas
# optimization_sf.R
# Written on June 14, 2017 
# Modified on June 14, 2017
#
# Script in order to determine best fit parameters using optimization 
library(DEoptim)

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data(1939, 2014)

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
tide.data <- read.tide.data()

#log likelihood
log.like.calc <- function(data, location1, scale1, shape1){
  ll<-sum(devd(data, loc=location1, scale=scale1, shape=shape1, type=c('GEV'), log=TRUE))
  return(ll)
}

#negative log likelihood calculator 
neg.log.like.calc <- function(p, parnames, data, temps){
  #checks length of parameters, which will be used to determine which parameters are non stationary 
  n.parnames <- length(parnames)
  
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
  if (is.na(nll)){
    nll <- 0 
  }
  return(nll)
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
