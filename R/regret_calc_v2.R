# Alex Klufas 
# June 27 2017 
# regret_calc_v2.R

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
damage.calc <- function(x, year){
  if (x > 2.850 - (year-1)*(.00257)){
    val <- 1.951*exp(.4116*x) / ((1.035)^(year-1)) #/ upkeep.cost #in meters, in millions of dollars 
  }
  else{
    val <- 0 
  }
  return(val)
}

yr <- seq(from = 0, to = 10 , by= 1)
sl <- seq(from = 2.5, to =3, by =.01)
test <- matrix(ncol = length(sl), nrow =length(yr) )
for (j in 1:length(yr)){
  for (i in 1:length(sl)){
    test[j, i] <- damage.calc(sl[i], yr[j])
  }
}
#calculates damage with barrier present 
damage.calc.w.barrier <- function(x, year){
  if (x > (3.28 - (year-1)*.00257) ){ #slr hard coded at the moment 
    val <- 1.951*exp(.4116*x) /((1.035)^(year-1)) #/ upkeep.cost #in m, in millions of dollars
  }
  else{
    val <- 0 
  }
  return(val)
}

sea.level.m <- seq(from=0, to = 10.0, by = .05)

sea.level.mm <- seq(from=0, to = 10000, by = 50)

#calculate damage w no barrier, make matrix 
damage.values <- rep(0, length(sea.level.m)) 
for (i in 1:length(sea.level.m)){
  damage.values[i] <- damage.calc(sea.level.m[i])*1e6
}

#calculate damage w barrier, make matrix 
damage.values.w.barrier <- rep(0, length(sea.level.m)) #len = 21 
for (i in 1:length(sea.level.m)){
  damage.values.w.barrier[i] <- damage.calc.w.barrier(sea.level.m[i])*1e6
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
  #hgt.changing <- seq(from = 2850, to = 2796 , by = -2.57 ) #this is hard coded at the moment to try to get it to work 
  #hgt.barrier.changing <- seq(from = 3280, to = 3227, by = -2.57)
  
  n.years <- length(years)
  
  #calculate cdf of flooding w no barrier and sea level  
  #temp included if non stationary 
  #prob.of.flood <- matrix(ncol = length(years), nrow = niter)
  #new.vals.barrier <- matrix(ncol = length(years), nrow = niter)
  
  #prob.of.flood.barrier <- matrix(ncol = length(years), nrow = niter)
  #new.vals.barrier <- matrix(ncol = length(years), nrow = niter)
  #prob.of.flood <- rep(0,length(sea.level))
  #prob.of.flood.barrier <- rep(0,length(sea.level))
  
  total.expected.no.barrier <- rep(0, niter)
  total.expected.w.barrier <- rep(0,niter)
  
  #new.vals.no.barrier <- rep(0,length(sea.level))
  #new.vals.barrer <- rep(0,length(sea.level))

  #new.vals.no.barrier <- rep(0, niter)
  #new.vals.barrier <- rep(0, niter)
  
  #expected.damage.no.barrier <- rep(0,n.years)#0,n.years)
  #expected.damage.w.barrier <- rep(0,n.years)
  expected.total.costs.no.barrier <- rep(0,niter)
  expected.total.costs.w.barrier <- rep(0,niter)
  
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for(i in 1:niter){ #cutting off the start before convergence , need temp prediction of the future , going to make it shorter to see how long this takes to do 
    
    #uncomment when doing more than one iteration 
    expected.damage.no.barrier <- rep (0,n.years)
    expected.damage.w.barrier <- rep(0,n.years)
    total.costs.no.barrier <- rep (0,n.years)
    total.costs.w.barrier <- rep (0,n.years)
     for (t in 1:n.years){
       #new.vals.no.barrier[i] <- devd.calculator(p= mcmc.chain[5e4+i,], param.names = param.names, temps[t], data = sea.level) #hgt.changing
       #new.vals.barrier[i] <- devd.calculator(p= mcmc.chain[5e4+i,], param.names = param.names, temps[t], data = sea.level) #hgt.changing
     # prob.of.flood[i,t] <- pevd.calculator(p=mcmc.chain[5e4+i,], param.names=param.names, temperature = temps[t], data.pt=hgt.changing[t]) 
     # new.vals.no.barrier[t] <- devd.calculator( p= mcmc.chain[(i + 5e4),], param.names = param.names, temps[t], data = hgt.changing)
      
      #prob.of.flood.barrier[i,t] <- pevd.calculator(p=mcmc.chain[5e4+i,], param.names=param.names, temperature =temps[t], data.pt=hgt.barrier.changing[t]) 
     # new.vals.barrier[t] <- devd.calculator(p= mcmc.chain[(i + 5e4),], param.names = param.names, temps[t], data = hgt.barrier.changing)
      
      #uncomment when doing more than one iteration 
      #new.vals.no.barrier <- rep (0,length(sea.level))
      #new.vals.barrer <- rep(0,length(sea.level))
      #for (q in 1:length(sea.level.mm)){
        #pevd ? 
      new.vals.no.barrier <- devd.calculator(p= mcmc.chain[5e4+i,], param.names = param.names, temperature =temps[t], data = sea.level.mm) #hgt.changing
      new.vals.barrier <- devd.calculator(p= mcmc.chain[5e4+i,], param.names = param.names, temperature =temps[t], data = sea.level.mm) #hgt.changing
        #prob.of.flood[q] <- pevd.calculator(p=mcmc.chain[5e4+i,], param.names=param.names, temps[t], data.pt=hgt.changing[t]) 
        # new.vals.no.barrier[t] <- devd.calculator( p= mcmc.chain[(i + 5e4),], param.names = param.names, temps[t], data = hgt.changing)
        
        #prob.of.flood.barrier[q] <- pevd.calculator(p=mcmc.chain[5e4+i,], param.names=param.names, temps[t], data.pt=hgt.barrier.changing[t]) 
      #}
      
      #calculating changing damage costs as sea level rises
      damage.values <- rep(0, length(sea.level.m))
      damage.values.w.barrier <- rep(0, length(sea.level.m)) 
      for (j in 1:length(sea.level.m)){
        damage.values[j] <- damage.calc(sea.level.m[j], t)*1e6  
        damage.values.w.barrier[j] <- damage.calc.w.barrier(sea.level.m[j], t)*1e6
      }
     # print(damage.values)
      print(damage.values.w.barrier)
      
      expected.damage.no.barrier[t] <- sum(damage.values * new.vals.no.barrier) / sum(new.vals.no.barrier)
      expected.damage.w.barrier[t] <- sum(damage.values.w.barrier * new.vals.barrier) / sum(new.vals.barrier)
      total.costs.no.barrier[t] <- expected.damage.no.barrier[t]
      
      
      if (t == 1){
        total.costs.w.barrier[t] <- expected.damage.w.barrier[t] + 18e6
      }
      else{
        total.costs.w.barrier[t] <- expected.damage.w.barrier[t] + (75e3 / ((1.035)^(t-1)))
      }
        # expected.damage + build/maintenance costs
     }
    #expected damage costs
    print(expected.damage.no.barrier - expected.damage.w.barrier)
    #print(expected.damage.w.barrier)
    
    total.expected.no.barrier[i] <- mean(expected.damage.no.barrier)
    total.expected.w.barrier[i] <- mean(expected.damage.w.barrier)
    
    #expected costs including building & upkeeping barrier 
    expected.total.costs.no.barrier[i] <- mean(total.costs.no.barrier)
    expected.total.costs.w.barrier[i] <- mean(total.costs.w.barrier)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #flood.damage.no.barrier <- expected.total.costs.no.barrier * prob.of.flood #changed from total.expected.no.barrier
  #flood.damage.w.barrier <- expected.total.costs.w.barrier * prob.of.flood.barrier  #changed from total.expected.w.barrier
  
  #print(flood.damage.no.barrier)
  #print(flood.damage.w.barrier)
  #find the average values of damage costs for each iteration 
  
  #avg.damage.no.barrier <- mean(expected.total.costs.no.barrier) #apply(flood.damage.no.barrier, 1, mean) #mean(flood.damage.no.barrier)
  #avg.damage.w.barrier <- mean(expected.total.costs.w.barrier)#apply(flood.damage.w.barrier, 1, mean) #mean(flood.damage.w.barrier)
  
  #use which to find the number of places where we should and shouldnt build 
  should.build <- which(expected.total.costs.w.barrier < expected.total.costs.no.barrier)
  shouldnt.build <- which(expected.total.costs.no.barrier < expected.total.costs.w.barrier)
  #should.build <- which(avg.damage.w.barrier < avg.damage.no.barrier)
  #shouldnt.build <- which(avg.damage.no.barrier < avg.damage.w.barrier)
  
 # regret.vals <- rep(0, length(avg.damage.no.barrier))
  regret.vals <- rep(0, length(expected.total.costs.w.barrier))
  if (length(should.build) > length(shouldnt.build)){
    #we will be building - now need to calculate regret of the no build values 
    for (i in 1:length(shouldnt.build)){
      regret.vals[i] <- abs(expected.total.costs.no.barrier[shouldnt.build[i]] - expected.total.costs.w.barrier[shouldnt.build[i]])
     # regret.vals[i] <- abs(avg.damage.no.barrier[shouldnt.build[i]] - avg.damage.w.barrier[shouldnt.build[i]])
    }
  }else{
    for (i in 1:length(should.build)){
      regret.vals[i] <- abs(expected.total.costs.no.barrier[should.build[i]] - expected.total.costs.w.barrier[should.build[i]])
      #regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
    }
    #regret.vals[i] <- abs(avg.damage.no.barrier[should.build[i]] - avg.damage.w.barrier[should.build[i]])
  }
  
  regret.calc.out <- vector('list', 11)
  
  names(regret.calc.out) <- c('regret.array', 'should.build', 
                              'shouldnt.build', 
                              'new.vals.no.barrier', 'new.vals.barrier','expected.damage.no.barrier',
                              'expected.damage.w.barrier', 'total.w.barrier', 'total.no.barrier', 
                              'expected.total.costs.no.barrier', 'expected.total.costs.w.barrier')
  #taken out: 'no.barrier', 'w.barrier','flood.prob', 'flood.w.barrier.prob', 
  
 # regret.calc.out$no.barrier <- flood.damage.no.barrier
  
  #regret.calc.out$w.barrier <- flood.damage.w.barrier
  
  regret.calc.out$regret.array <- regret.vals
  
  regret.calc.out$should.build <- should.build
  
  regret.calc.out$shouldnt.build <- shouldnt.build 
  
 # regret.calc.out$avg.no.barrier <- avg.damage.no.barrier
  
 # regret.calc.out$avg.w.barrier <- avg.damage.w.barrier
  
#  regret.calc.out$flood.prob <- prob.of.flood
  
 # regret.calc.out$flood.w.barrier.prob <- prob.of.flood.barrier
  
  regret.calc.out$new.vals.no.barrier <- new.vals.no.barrier
  
  regret.calc.out$new.vals.barrier <- new.vals.barrier
  
  regret.calc.out$expected.damage.no.barrier <- expected.damage.no.barrier
  
  regret.calc.out$expected.damage.w.barrier <- expected.damage.w.barrier
  
  regret.calc.out$total.w.barrier <- total.expected.w.barrier
  
  regret.calc.out$total.no.barrier <- total.expected.no.barrier
  
  regret.calc.out$expected.total.costs.no.barrier <- expected.total.costs.no.barrier
  
  regret.calc.out$expected.total.costs.w.barrier <- expected.total.costs.w.barrier
  
  return(regret.calc.out)
  
}

setwd('~/codes/Klufas_NewLondon/')
#stuff to call for this :
temp.vals <- temperature_proj[166:186]
niter <- 10
years <- seq(from=2015, to = 2035, by =1)
effect.height.barrier <- 3280 
effect.height.no.barrier <- 2850
slr <- 2.57
#source('adaptive_mcmc_nonstat.R')
load('mcmc.test.stationary.parallel.RData') # <<<<<<<<<<<<<<<<<<<<<<<<< messing up 'damage.values' (and maybe some other stuff) FIX ME
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
