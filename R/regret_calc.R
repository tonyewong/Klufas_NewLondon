# Alex Klufas 
# June 27 2017 
# regret_calc.R

library(mvtnorm)
library(extRemes)
library(coda)
library(ncdf4)

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data()

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
tide.data <- read.tide.data()

setwd('~/codes/Klufas_NewLondon/R/')
source('gev_nonstationary_MCMC.R')


#add temp trends stuff 

setwd('~/codes/Klufas_NewLondon/R/')
ncdata <- nc_open('../global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc')
temperature_proj <- ncvar_get(ncdata, 'tas')
time_proj <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

time_proj <- floor(time_proj/10000)

ind_norm <- which(time_proj==1901):which(time_proj==2000)
temperature_proj <- temperature_proj - mean(temperature_proj[ind_norm])

#------------------------some regret stuff--------------------------
#build or not to build calculator 
niter <- 1e5
build <- rep(0,niter)
regret.sow <- rep(0,niter)
initial.build.cost <- 18e6 

i <- seq(from = 1, to = 20, by = 1)

upkeep.cost <- sum(75e3 / ((1 +.035)^(i-1)))

#calculates pevd with given parameters 
pevd.calculator <- function(p, param.names, temperature, data.pt){
  n.parnames <- length(param.names)
  #checks length of parameters, which will be used to determine which parameters are non stationary 
  
  #if all parameters are stationary 
  if (n.parnames == 3){
    mu <- p[1]
    sigma <- p[2]
    xi <- p[3]
    #print('here1')
  }
  
  #one stationary parameter
  else if(n.parnames ==4){
    if (param.names[1] == 'mu0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi <- p[4]
     # print('here1')
      mu <- mu0 + mu1*temperature
      
    }
    else if (param.names[2] == 'sigma0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi <- p[4]
     # print('here2')
      sigma <- exp(sigma0+ sigma1*temperature)
    }
    else{
      mu <- p[1]
      sigma <- p[2]
      xi0 <- p[3]
      xi1 <- p[4]
     #print('here2')
     
      xi <- xi0 + xi1*temperature
    }
  }
  
  #two stationary parameters
  else if (n.parnames == 5){
    if (param.names[1] == 'mu0' & param.names[3]=='sigma0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma0 <- p[3]
      sigma1 <- p[4]
      xi <- p[5]
     #print('here4')
      mu <- mu0 + mu1*temperature
      sigma <- exp(sigma0 + sigma1*temperature)
    }
    else if(parnames[2] == 'sigma0' &param.names[4] == 'xi0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      sigma <- exp(sigma0 + sigma1*temperature)
      xi <- xi0 + xi1*temperature
      
    }
    else{
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      
      mu <- mu0 + mu1*temperature
      xi <- xi0 + xi1*temperature
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
    mu <- mu0 + mu1*temperature # I think it is freaking out because it doesn't know what temperature to use here well....
    
    sigma <- exp(sigma0 + sigma1*temperature)
    
    xi <- xi0 + xi1*temperature
  }
  val <- 1 - pevd(data.pt, loc = mu , scale = sigma, shape = xi, type=c('GEV'))
  return(val)
}

devd.calculator <- function(p, param.names, temperature, data){
  n.parnames <- length( param.names)
  #checks length of parameters, which will be used to determine which parameters are non stationary 
  
  #if all parameters are stationary 
  if (n.parnames == 3){
    mu <- p[1]
    sigma <- p[2]
    xi <- p[3]
    #print('here')
  }
  
  #one stationary parameter
  else if(n.parnames ==4){
    if (param.names[1] == 'mu0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi <- p[4]
      
      mu <- mu0 + mu1*temperature
      
    }
    else if (param.names[2] == 'sigma0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi <- p[4]
      
      sigma <- exp(sigma0+ sigma1*temperature)
    }
    else{
      mu <- p[1]
      sigma <- p[2]
      xi0 <- p[3]
      xi1 <- p[4]
      
     
      xi <- xi0 + xi1*temperature
    }
  }
  
  #two stationary parameters
  else if (n.parnames == 5){
    if (param.names[1] == 'mu0' & param.names[3]=='sigma0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma0 <- p[3]
      sigma1 <- p[4]
      xi <- p[5]
      
      mu <- mu0 + mu1*temperature
      sigma <- exp(sigma0 + sigma1*temperature)
    }
    else if(param.names[2] == 'sigma0' & param.names[4] == 'xi0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      
      sigma <- exp(sigma0 + sigma1*temperature)
      xi <- xi0 + xi1*temperature
      
    }
    else{
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      
      mu <- mu0 + mu1*temperature
      xi <- xi0 + xi1*temperature
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
    
    mu <- mu0 + mu1*temperature # I think it is freaking out because it doesn't know what temperature to use here well....
    
    sigma <- exp(sigma0 + sigma1*temperature)
    
    xi <- xi0 + xi1*temperature
  }
  val <- devd(data, loc = mu , scale = sigma, shape = xi, type=c('GEV'))
  return(val)
}

#calculates damage when no barrier 
damage.calc <- function(x){
  if (x > 2.850){
    val <- 1e6*1.951*exp(.4116*x) #/ upkeep.cost #in meters, in millions of dollars 
  }
  else{
    val <- 0 
  }
  return(val)
}


#calculates damage with barrier present 
damage.calc.w.barrier <- function(x){
  if (x > 3.28){
    val <- 1e6*1.951*exp(.4116*x) #/ upkeep.cost #in m, in millions of dollars
  }
  else{
    val <- 0 
  }
  return(val)
}

#will ultimately figure out whether to build or not and then calculate the regret depending on the reponse 
regret.calculator <- function(niter, 
                              years, 
                              effect.height.no.barrier, 
                              effect.height.barrier, 
                              slr, 
                              mcmc.chain, 
                              n.params, 
                              param.names, 
                              temps
                              ){ #niter = number iterations (int), years - sequence of years to cover, 
#prob of a certain flood happening                    #effect.height.no.barrier (height in mm of no damage area)
                                                      #effect.height.barrier (height in mm of no damage area w barrier present)
                                                       #slr, rate at which sea level is changing in the particular area 
  #change in effective height over course of years 
  #lower.level <- effect.height.no.barrier - 21*slr
  #lower.level.barrier <- effect.height.barrier - 21*slr
 # hgt.changing <- effect.height.no.barrier - length(years)*slr
  #hgt.barrier.changing <- effect.height.barrier - length(years)*slr
  hgt.changing <- seq(from = 2850, to = 2808, by = -2 ) #this is hard coded at the moment to try to get it to work 
  hgt.barrier.changing <- seq(from = 3280, to = 3238, by = -2)
  
  n.years <- length(years)
  
  #calculate cdf of flooding w no barrier and sea level  
  #temp included if non stationary 
  prob.of.flood <- matrix(ncol = length(years), nrow = niter)
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for(i in 5e4:niter){ #cutting off the start before convergence , need temp prediction of the future 
    for (t in 1:n.years){
      prob.of.flood[i,t] <- pevd.calculator(p=mcmc.chain[i,], param.names=param.names, temps[t], data.pt=hgt.changing[t]) 
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #calculate cdf of flooding w barrier present and sea level rise 
  #temp included if non stationary 
  prob.of.flood.barrier <- matrix(ncol = length(years), nrow = niter)
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for(i in 5e4:niter){
    for (t in 1:n.years){
      prob.of.flood.barrier[i,t] <- pevd.calculator(p=mcmc.chain[i,], param.names=param.names, temps[t], data.pt=hgt.barrier.changing[t]) 
    }
    setTxtProgressBar(pb, i)
  }
  
  #calculating the gev of each set of parameters - w barrier 
  new.vals.barrier <- matrix(ncol = length(years), nrow = niter)
  for(i in 5e4:niter){
    for (t in 1:n.years){
      new.vals.barrier[i,t] <- devd.calculator(p= mcmc.chain[i,], param.names = param.names, temps[t], data = hgt.barrier.changing)
    }
  }
  #no barrier 
  new.vals.barrier <- matrix(ncol = length(years), nrow = niter)
  for (t in 1:n.years){
    new.vals.no.barrier[t] <- devd.calculator( p= mcmc.chain[t,], param.names = param.names, temps[t], data = hgt.changing)
  }
  
  
  #calculating damage at each height on the spectrum (0 to 4.75m)
  sea.level <- seq(from=0, to = 10.0, by = .05) #len = 20 
  
  damage.values <- rep(0, length(sea.level)) #len = 20 
  
  #also need to calculate expected flood damage in each year 
  expect.flood.damage.no.barrier <- rep(0, n.years)
  expect.flood.damage.w.barrier <- rep(0, n.years)
  
  for (i in 1:length(sea.level)){
    damage.values[i] <- damage.calc(sea.level[i])
    #expect.flood.damage.no.barrier[i] <- sum(damage.values[i]* new.vals.no.barrier) / sum(new.vals.no.barrier)
  }
  
  expect.value.no.barrier <- sum(damage.values * new.vals.no.barrier) / sum(new.vals.no.barrier)
  
  damage.values.w.barrier <- rep(0, length(sea.level)) #len = 21 
  for (i in 1:length(sea.level)){
    damage.values.w.barrier[i] <- damage.calc.w.barrier(sea.level[i])
    #expect.flood.damage.w.barrier[i] <- sum(damage.values.w.barrier[i]* new.vals.barrier) / sum(new.vals.barrier)
  }
  
  expect.value.w.barrier <- sum(damage.values.w.barrier*new.vals.barrier) / sum(new.vals.barrier)
  
  #steps between the years 
  dx <- .05
  
  #calculating the avergae value of damage w and w/o barrier 
 # avg.w.barrier <- sum(damage.values.w.barrier*new.vals.barrier) / sum(new.vals.barrier) #took out dx temporarily 
 # avg.wo.barrier <- sum(damage.values*new.vals.no.barrier) / sum(new.vals.no.barrier)
  avg.w.barrier <- mean(expect.value.w.barrier)
  avg.wo.barrier <- mean(expect.value.no.barrier)
  #calculating flood damage overall w and w/o barrier 
  flood.damage.no.barrier <- avg.wo.barrier * prob.of.flood[5e4:niter,]
  flood.damage.w.barrier <- avg.w.barrier * prob.of.flood.barrier[5e4:niter,] + initial.build.cost
  
  #find the average values of damage costs for each iteration 
  avg.damage.no.barrier <- apply(flood.damage.no.barrier, 1, mean)
  avg.damage.w.barrier <- apply(flood.damage.w.barrier, 1, mean)
  
  #use which to find the number of places where we should and shouldnt build 
  should.build <- which(avg.damage.w.barrier < avg.damage.no.barrier)
  shouldnt.build <- which(avg.damage.no.barrier < avg.damage.w.barrier)
  
  regret.vals <- rep(0, length(avg.damage.no.barrier))
  if (length(should.build) > length(shouldnt.build)){
    #we will be building - now need to calculate regret of the no build values 
    for (i in 1:length(shouldnt.build)){
      regret.vals[i] <- abs(avg.damage.no.barrier[shouldnt.build[i]] - avg.damage.w.barrier[shouldnt.build[i]])
    }
  }else{
    for (i in 1:length(should.build)){
     regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
    }
    #regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
  }
  
  regret.calc.out <- vector('list', 9)
  
  names(regret.calc.out) <- c('no.barrier', 'w.barrier', 'regret.array', 'should.build', 
                              'shouldnt.build', 'avg.no.barrier','avg.w.barrier', 'flood.prob', 'flood.w.barrier.prob')
  
  regret.calc.out$no.barrier <- flood.damage.no.barrier[,1:20]
  
  regret.calc.out$w.barrier <- flood.damage.w.barrier[,1:20]
  
  regret.calc.out$regret.array <- regret.vals
  
  regret.calc.out$should.build <- should.build
  
  regret.calc.out$shouldnt.build <- shouldnt.build 
  
  regret.calc.out$avg.no.barrier <- avg.damage.no.barrier
  
  regret.calc.out$avg.w.barrier <- avg.damage.w.barrier
  
  regret.calc.out$flood.prob <- prob.of.flood[5e4:niter,]
  
  regret.calc.out$flood.w.barrier.prob <- prob.of.flood.barrier[5e4:niter,]
  
  return(regret.calc.out)
  
}

setwd('~/codes/Klufas_NewLondon/')
#stuff to call for this :
temp.vals <- temperature_proj[166:186]
niter <- 1e5
years <- seq(from=2015, to = 2035, by =1)
effect.height.barrier <- 3280 
effect.height.no.barrier <- 2850
slr <- 2
#source('adaptive_mcmc_nonstat.R')
load('mcmc.test.stationary.parallel.RData')
param.names <- c('mu', 'sigma', 'xi')
test1 <- mcmc.test.stationary.parallel[[1]]$samples

regret.stationary <- regret.calculator(niter = niter, 
                                       years = years, 
                                       effect.height.no.barrier = effect.height.no.barrier,
                                       effect.height.barrier=effect.height.barrier,  
                                       slr=slr, 
                                       mcmc.chain = test1 , 
                                       n.params = 3, 
                                       param.names = param.names,  
                                       temps = temp.vals)



avg.damage.no.barrier <- apply(regret.stationary$no.barrier, 1, mean)
avg.damage.w.barrier <- apply(regret.stationary$w.barrier, 1, mean)

#use which to find the number of places where we should and shouldnt build 
should.build <- which(avg.damage.w.barrier < avg.damage.no.barrier)
shouldnt.build <- which(avg.damage.no.barrier < avg.damage.w.barrier)

regret.vals <- rep(0, length(avg.damage.no.barrier))
if (length(should.build) > length(shouldnt.build)){
  #we will be building - now need to calculate regret of the no build values 
  for (i in 1:length(shouldnt.build)){
    regret.vals[i] <- abs(avg.damage.no.barrier[shouldnt.build[i]] - avg.damage.w.barrier[shouldnt.build[i]])
  }
}else{
  for (i in 1:length(should.build)){
    regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
  }
  #regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
}




#sigma nonstationary 
sigma.nonstat <- mcmc.sigma.nonstat.parallel[[1]]$samples
param.names <- c('mu', 'sigma0', 'sigma1', 'xi')

regret.sigma.nonstationary <- regret.calculator(niter = niter, 
                                              years = years, 
                                              effect.height.no.barrier = effect.height.no.barrier,
                                              effect.height.barrier=effect.height.barrier,  
                                              slr=slr, 
                                              mcmc.chain = sigma.nonstat, 
                                              n.params = 4, 
                                              param.names = param.names,  
                                              temps = temp.vals)


avg.damage.no.barrier <- apply(regret.sigma.nonstationary$no.barrier, 1, mean)
avg.damage.w.barrier <- apply(regret.sigma.nonstationary$w.barrier, 1, mean)

#use which to find the number of places where we should and shouldnt build 
should.build <- which(avg.damage.w.barrier < avg.damage.no.barrier)
shouldnt.build <- which(avg.damage.no.barrier < avg.damage.w.barrier)

regret.vals <- abs(avg.damage.no.barrier - avg.damage.w.barrier)
  
  #rep(0, length(avg.damage.no.barrier))
if (length(should.build) > length(shouldnt.build)){
  #we will be building - now need to calculate regret of the no build values 
  for (i in 1:length(shouldnt.build)){
    regret.vals[i] <- abs(avg.damage.no.barrier[shouldnt.build[i]] - avg.damage.w.barrier[shouldnt.build[i]])
  }
}else{
  for (i in 1:length(should.build)){
    regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
  }
  #regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
}



#xi non stationary 
xi.nonstat <- mcmc.xi.nonstat.parallel[[1]]$samples
param.names <- c('mu', 'sigma', 'xi0', 'xi1')
regret.xi.nonstationary <- regret.calculator(niter = niter, 
                                              years = years, 
                                              effect.height.no.barrier = effect.height.no.barrier,
                                              effect.height.barrier=effect.height.barrier,  
                                              slr=slr, 
                                              mcmc.chain = xi.nonstat, 
                                              n.params = 4, 
                                              param.names = param.names,  
                                              temps = temp.vals)

#mu non stationary
mu.nonstat <- mcmc.mu.nonstat.parallel[[1]]$samples
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
regret.mu.nonstationary <- regret.calculator(niter = niter, 
                                              years = years, 
                                              effect.height.no.barrier = effect.height.no.barrier,
                                              effect.height.barrier=effect.height.barrier,  
                                              slr=slr, 
                                              mcmc.chain = mu.nonstat, 
                                              n.params = 4, 
                                              param.names = param.names,  
                                              temps = temp.vals)

#sigma and xi non stationary 
sigma.xi.nonstat <- mcmc.sigma.xi.nonstat.parallel[[1]]$samples
param.names <- c('mu', 'sigma0', 'sigma1', 'xi0', 'xi1')
regret.sigma.xi.nonstationary <- regret.calculator(niter = niter, 
                                              years = years, 
                                              effect.height.no.barrier = effect.height.no.barrier,
                                              effect.height.barrier=effect.height.barrier,  
                                              slr=slr, 
                                              mcmc.chain = sigma.xi.nonstat, 
                                              n.params = 5, 
                                              param.names = param.names,  
                                              temps = temp.vals)
#mu and xi non stationary 
mu.xi.nonstat <- mcmc.mu.xi.nonstat.parallel[[1]]$samples
param.names <- c('mu0', 'mu1', 'sigma', 'xi0', 'xi1')
regret.mu.xi.nonstationary <- regret.calculator(niter = niter, 
                                              years = years, 
                                              effect.height.no.barrier = effect.height.no.barrier,
                                              effect.height.barrier=effect.height.barrier,  
                                              slr=slr, 
                                              mcmc.chain = mu.xi.nonstat, 
                                              n.params = 5, 
                                              param.names = param.names,  
                                              temps = temp.vals)

#mu and sigma nonstationary 
mu.sigma.nonstat <- mcmc.mu.sigma.nonstat.parallel[[1]]$samples
param.names <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi')
regret.mu.sigma.nonstationary <- regret.calculator(niter = niter, 
                                              years = years, 
                                              effect.height.no.barrier = effect.height.no.barrier,
                                              effect.height.barrier=effect.height.barrier,  
                                              slr=slr, 
                                              mcmc.chain = mu.sigma.nonstat, 
                                              n.params = 5, 
                                              param.names = param.names,  
                                              temps = temp.vals)


#all non stationary 
test2 <- mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples
param.names <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi0', 'xi1')

regret.all.nonstationary <- regret.calculator(niter = niter, 
                                             years = years, 
                                             effect.height.no.barrier = effect.height.no.barrier,
                                             effect.height.barrier=effect.height.barrier,  
                                             slr=slr, 
                                             mcmc.chain = test2, 
                                             n.params = 6, 
                                             param.names = param.names,  
                                             temps = temp.vals)


hist(regret.stationary$regret.array[,2])
hist(regret.mu.nonstationary$regret.array[,2])
hist(regret.mu.xi.nonstationary$regret.array[,2])
hist(regret.mu.sigma.nonstationary$regret.array[,2])




for (i in 1:niter){
  if (mean(flood.damage.no.barrier[i,]) > mean(flood.damage.w.barrier[i,])){
    #build the barrier 
  }
  else if (mean(flood.damage.no.barrier[i,]) < mean(flood.damage.w.barrier[i,])){
    #flood damage w barrier greater than flood damage w/o barrier 
    # then we don't build 
  }
}








#2.85 effective height 
#making synthetic slr data 
slr <- seq (from = 0, to = 50, by = 2)

effect.height <- seq(from = 2850 , to = 2810, by = -2)  #in mm b/c data

non.stat.loc <- run1.mu.nonstat$chains[,1] + run1.mu.nonstat$chains[,2]*temps$values

mcmc.scale <- run1.mu.nonstat$chains[,3]
  
mcmc.shape <- run1.mu.nonstat$chains[,4]

temp.vals <- temperature_proj[166:186]

prob.of.flood <- matrix(ncol = 20, nrow = niter)
#still trying to figure out temperature in this all 
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for(i in 1e4:niter){ #cutting off the start before convergence , need temp prediction of the future 
    for (t in 1:20){
        prob.of.flood[i,t] <- 1 - pevd(effect.height[t], 
                                       loc = run1.mu.xi.nonstat$chains[i,1] + run1.mu.xi.nonstat$chains[i,2]*temp.vals[t], 
                                       scale = exp(run1.mu.xi.nonstat$chains[i,3] + run1.mu.xi.nonstat$chains[i,4]*temp.vals[t]), 
                                       shape = run1.mu.xi.nonstat$chains[i,5] + run1.mu.xi.nonstat$chains[i,6]*temp.vals[t],
                                       type=c('GEV'))
    }
  setTxtProgressBar(pb, i)
}
close(pb)

#3.28 effective height 
prob.of.flood.barrier <- matrix(ncol = 20, nrow = niter)
effect.height.barrier <- seq(from = 3280, to = 3240, by = -2)#in mm b/c data
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for(i in 1e4:niter){
  for (t in 1:20){
      prob.of.flood.barrier[i,t] <- 1 - pevd(effect.height.barrier[t], loc = run1.mu.xi.nonstat$chains[i,1] + run1.mu.xi.nonstat$chains[i,2]*temp.vals[t], 
                                         scale = exp(run1.mu.xi.nonstat$chains[i,3] + run1.mu.xi.nonstat$chains[i,4]*temp.vals[t]), 
                                         shape = run1.mu.xi.nonstat$chains[i,5] + run1.mu.xi.nonstat$chains[i,6]*temp.vals[t],
                                          type=c('GEV'))
  }
  setTxtProgressBar(pb, i)
}
close(pb)

hist(prob.of.flood.barrier, freq = FALSE)
hist(prob.of.flood, freq = FALSE)

#------flood damage costs------------------------------------------

#no decrease guture rate thing 
damage.calc <- function(x){
  if (x > 2.850){
     val <- (1.951e6 * exp(.4116*x)) / upkeep.cost #in m, in millions of dollars
  }
  else{
    val <- 0 
  }
  return(val)
}

damage.calc.w.barrier <- function(x){
  if (x > 3.28){
    val <- (1.951e6 * exp(.4116*x)) / upkeep.cost #in m, in millions of dollars
  }
  else{
    val <- 0 
  }
  return(val)
}

dx <- seq(from = 0 , to = 20, by = 1)

years <- seq(from = 2015, to = 2035, by = 1)



#calculating the damages at certain heights w/o barrier
ss.stuff <- seq(from=0, to=4.75, by=.250)

damage.values <- rep(0, length(ss.stuff))

for (i in 1:length(ss.stuff)){
  damage.values[i] <- damage.calc(ss.stuff[i])*1e6
}

#calculating the damages at certain heights w/ barrier
damage.values.w.barrier <- rep(0, length(ss.stuff))

for (i in 1:length(ss.stuff)){
  damage.values.w.barrier[i] <- damage.calc.w.barrier(ss.stuff[i])*1e6
}

#trying to calculate overall avg value of damage 
full.damage.calc <- function(barrier, hgt){
  if (barrier == TRUE){
    damage.cost <- damage.calc.w.barrier(hgt)
  }
  else{
    damage.cost <- damage.calc(hgt)
  }
  return(damage.cost)
}

new.vals.barrier <- rep(0, 20)
new.vals.no.barrier <- rep(0,20)

dx <- seq(from = 0 , to = 19, by = 1)

#calculating the gev of each set of parameters ? I think - with barrier
for (t in 1:20){
  new.vals.barrier[t] <- sum(devd(effect.height.barrier, loc = run1.mu.xi.nonstat$chains[1e4:niter,1] + run1.mu.xi.nonstat$chains[1e4:niter,2]*temp.vals[t], 
                          scale = exp(run1.mu.xi.nonstat$chains[1e4:niter,3] + run1.mu.xi.nonstat$chains[1e4:niter,4]*temp.vals[t]), 
                          shape = run1.mu.xi.nonstat$chains[1e4:niter,5] + run1.mu.xi.nonstat$chains[1e4:niter,6]*temp.vals[t],
                          type=c('GEV')))
}
#calculating the gev of each set of parameters ? I think - w/o barrier
for (t in 1:20){
  new.vals.no.barrier[t] <- sum(devd(effect.height, loc = run1.mu.xi.nonstat$chains[1e4:niter,1] + run1.mu.xi.nonstat$chains[1e4:niter,2]*temp.vals[t], 
                          scale = exp(run1.mu.xi.nonstat$chains[1e4:niter,3] + run1.mu.xi.nonstat$chains[1e4:niter,4]*temp.vals[t]), 
                          shape = run1.mu.xi.nonstat$chains[1e4:niter,5] + run1.mu.xi.nonstat$chains[1e4:niter,6]*temp.vals[t],
                          type=c('GEV')))
}

#trying to calculate costs w/ and w/o barrier over all with this set of params 

avg.w.barrier <- sum(damage.values.w.barrier*new.vals.barrier*dx)

avg.wo.barrier <- sum(damage.values*new.vals.no.barrier*dx)

#multiply these values by the cdfs I think , I think something weird is going on here 

flood.damage.no.barrier <- avg.wo.barrier * prob.of.flood[1e4:niter,]
  
flood.damage.w.barrier <- avg.w.barrier * prob.of.flood.barrier[1e4:niter,] + initial.build.cost + upkeep.cost

#cost of upkeep only needed for part w barrier 

hist(flood.damage.no.barrier)
abline(v=mean(flood.damage.no.barrier), col='red')
mean1 <- mean(flood.damage.no.barrier)

hist(flood.damage.w.barrier)
abline(v=mean(flood.damage.w.barrier), col='red')
mean2 <- mean(flood.damage.w.barrier)


#-------------------------------------------------------------------
# pb <- txtProgressBar(min=0,max=time,initial=0,style=3)
#for (i in 1:time){

#  if (mean(flood.damage.no.barrier[i,]) > mean(flood.damage.w.barrier[i,])){
#build barrier 
#   array[i, 1] <- 1
# }

#  else if (mean(flood.damage.no.barrier[i,]) < mean(flood.damage.w.barrier[i,])){
#don't build barrier
#   array[i, 1] <- 0
# }
#setTxtProgressBar(pb, i)
#}

# #to build or not to build decider 
#  if ((sum(array[,1]) / time) >= .5){
# #going to build 
#need to calculate regret for just the sow where you didn't want to build 
#  pb <- txtProgressBar(min=0,max=time,initial=0,style=3)
# for(i in 1:time){
#   if (array[i,1] == 0){
#    array[i,2] <- 0 
#  }
#  else{
#     array[i,2] <-  mean(flood.damage.no.barrier[i,]) - mean(flood.damage.w.barrier[i,]) 
#    }
#    setTxtProgressBar(pb, i)
#  }
# }else{
#  for(i in 1:time){
#     if (array[i,1] == 1){
#      array[i,2] <- 0 
#    }
#    else{
#     array[i,2] <-   mean(flood.damage.w.barrier[i,]) - mean(flood.damage.no.barrier[i,])
##    }
#    setTxtProgressBar(pb, i)
}
#do the opposite
# }

