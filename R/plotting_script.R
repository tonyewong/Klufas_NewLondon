#plotting script 

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

#setwd('~/codes/Klufas_NewLondon/R/')
#source('gev_nonstationary_MCMC.R')

setwd('~/codes/Klufas_NewLondon/R/')
load('mcmc.test.stationary.parallel.RData')
load('DEOptim.allnonstat.RData')

#getting temperature for future , starts in 1850 
setwd('~/codes/Klufas_NewLondon/R/')
ncdata <- nc_open('../global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc')
temperature_proj <- ncvar_get(ncdata, 'tas')
time_proj <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

time_proj <- floor(time_proj/10000)

ind_norm <- which(time_proj==1901):which(time_proj==2000)

temperature_proj <- temperature_proj - mean(temperature_proj[ind_norm])

#Non-Stationary Ex Plot - Only Mu Non Stat--
par(mfrow = c(1,1))
x.lsl <- seq(from =0, to =3000, by = 1)
load('mu.nonstat.deoptim.RData')

#load('DEOptim.allnonstat2.RData') #why are these coming out so weirdly???? 
non.stat.2015 <- devd(x.lsl, loc = optim.like.temp.mu$optim$bestmem[1]  + optim.like.temp.mu.sigma.xi$optim$bestmem[2] * temperature_proj[165],
                    scale = optim.like.temp.mu$optim$bestmem[3], 
                    shape = optim.like.temp.mu$optim$bestmem[4],
                    type = c("GEV"))
non.stat.2065 <-  devd(x.lsl, loc = optim.like.temp.mu$optim$bestmem[1]  + optim.like.temp.mu$optim$bestmem[2] * temperature_proj[215],
                       scale = optim.like.temp.mu$optim$bestmem[3], 
                       shape = optim.like.temp.mu$optim$bestmem[4] ,
                       type = c("GEV"))
non.stat.2095 <-  devd(x.lsl, loc = optim.like.temp.mu$optim$bestmem[1]  + optim.like.temp.mu$optim$bestmem[2] * temperature_proj[245],
                       scale = optim.like.temp.mu$optim$bestmem[3], 
                       shape = optim.like.temp.mu$optim$bestmem[4],
                       type = c("GEV"))

plot(non.stat.2015, type='l', col='black', lwd = 2, xlab = 'Sea Level Height [mm]', ylab = 'Probability')
lines(non.stat.2065, col='red', lwd = 2)
lines(non.stat.2095, col='blue', lwd = 2)
title('Distribution of Maxmimum Sea Level Heights with Location Parameter Non-Stationary')
legend(2500,.0022, legend= c('2015', '2065', '2095'), lty=c(1,1,1), lwd=c(2.5, 2.5,2.5), col=c('black', 'red','blue'))
#-----------------------------------------------


#Non-Stationary Ex Plot - All Params Non Stat--
load('DEOptim.allnonstat2.RData')
all.non.stat.2015 <- devd(x.lsl, loc = optim.like.temp.mu.sigma.xi$optim$bestmem[1]  + optim.like.temp.mu.sigma.xi$optim$bestmem[2] * temperature_proj[165],
                      scale = optim.like.temp.mu.sigma.xi$optim$bestmem[3] + optim.like.temp.mu.sigma.xi$optim$bestmem[4]*temperature_proj[165], 
                      shape = optim.like.temp.mu.sigma.xi$optim$bestmem[5] +  optim.like.temp.mu.sigma.xi$optim$bestmem[6]*temperature_proj[165],
                      type = c("GEV"))
all.non.stat.2065 <-  devd(x.lsl, loc = optim.like.temp.mu.sigma.xi$optim$bestmem[1]  + optim.like.temp.mu.sigma.xi$optim$bestmem[2] * temperature_proj[215],
                           scale = optim.like.temp.mu.sigma.xi$optim$bestmem[3] + optim.like.temp.mu.sigma.xi$optim$bestmem[4]*temperature_proj[215], 
                           shape = optim.like.temp.mu.sigma.xi$optim$bestmem[5] +  optim.like.temp.mu.sigma.xi$optim$bestmem[6]*temperature_proj[215],
                           type = c("GEV"))
all.non.stat.2095 <-  devd(x.lsl, loc = optim.like.temp.mu.sigma.xi$optim$bestmem[1]  + optim.like.temp.mu.sigma.xi$optim$bestmem[2] * temperature_proj[245],
                           scale = optim.like.temp.mu.sigma.xi$optim$bestmem[3] + optim.like.temp.mu.sigma.xi$optim$bestmem[4]*temperature_proj[245], 
                           shape = optim.like.temp.mu.sigma.xi$optim$bestmem[5] +  optim.like.temp.mu.sigma.xi$optim$bestmem[6]*temperature_proj[245],
                           type = c("GEV"))

plot(all.non.stat.2095, type='l', col='blue', lwd = 2,  xlab = 'Sea Level Height [mm]', ylab = 'Probability')
lines(all.non.stat.2065, col='red', lwd = 2)
lines(all.non.stat.2015, col='black', lwd = 2)
title('Distribution of Maxmimum Sea Level Heights with All Parameters Non-Stationary')
legend(2500,.0022, legend= c('2015', '2065', '2095'), lty=c(1,1,1), lwd=c(2.5, 2.5,2.5), col=c('black', 'red','blue'))
#-----------------------------------------------

load('mcmc.parallel.longer.iter.RData')

#MCMC Plot, All Params Stationary --------------
par(mfrow = c(3,1))
plot(mcmc.stationary.parallel[[1]]$samples[,1], 
     type='l',xlab = 'Iteration', main = 'Location, Scale, and Shape of Stationary MCMC Chains', ylab = 'Location')
plot(mcmc.stationary.parallel[[1]]$samples[,2], 
     type='l', xlab = 'Iteration',ylab = 'Scale')
plot(mcmc.stationary.parallel[[1]]$samples[,3], 
     type='l', xlab = 'Iteration', ylab = 'Shape')
#-----------------------------------------------

#MCMC Plot, Only Mu Non Stationary--------------
par(mfrow = c(4,1))
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,1], 
     type='l',xlab = 'Iteration', main = 'Scale, and Shape of Stationary, Location Non Stat MCMC Chains', ylab = 'Location')
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,2], 
     type='l', xlab = 'Iteration',ylab = 'Location * Temperature')
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,3], 
     type='l', xlab = 'Iteration', ylab = 'Scale')
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,3], 
     type='l', xlab = 'Iteration', ylab = 'Shape')
#-----------------------------------------------

#Beginnings of Mega Plot of All Prior and Posterior Params----------
#load('mcmc.params.RData')



load('mcmc.allnonstat.for.real.RData')
dev.off()
vals <- seq(from = 0, to = 3000, by = 1)
vals2 <- seq(from = -10, to = 10, by = 1)
#for locaiton , sigma 
mu.sig <- plot(vals, dunif(vals,min = 0, max = 3000), type = 'l') 


#par(mfrow = c(7,1))
#dud.row<- c(43,43,43,43,43,43,43)
row1 <- c(1, 2, 3, 4, 5, 6)
row2 <- c(7, 8,9,10,11,12)
row3 <- c(13,14,15,16,17,18)
row4 <- c(19,20,21,22,23,24)
row5 <- c(25,26,27,28,29,30)
row6 <- c(31,32,33,34,35,36)
row7 <- c(37,38,39,40,41,42)

m <- rbind(row1,  row2, row3, row4, row5, row6, row7)
layout(m)

par(mar = c(0,0,0,0))


#stationary case 

plot(density(mcmc.stationary.parallel.post.burnin.1[,1]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.stationary.parallel.post.burnin.1[,2]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.stationary.parallel.post.burnin.1[,3]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot.new()

#non stationary mu 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot.new()

#non stationary sigma 
par(mar = c(0,0,0,0))
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot.new()

#non stationary xi 
par(mar = c(0,0,0,0))
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')

#non stationary mu and sigma 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,5]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')


#non stationary mu and xi 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,5]), col = 'blue')
lines(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')

#non stationary sigma and xi 
par(mar = c(0,0,0,0))
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot.new()
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')

#all non stationary 
par(mar = c(0,0,0,0))
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(vals, dunif(vals,min = 0, max = 3000), type = 'l', col='red')
plot(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
plot(vals2, dunif(vals2, min = -10, max = 10), type = 'l', col='red')
#-------------------------------------------------------------------

