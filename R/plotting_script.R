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
ind2<- best.params.mu.nonstat
#load('DEOptim.allnonstat2.RData') #why are these coming out so weirdly???? 
non.stat.2015 <- devd(x.lsl, loc = mcmc.mu.nonstat.parallel[[1]]$samples[ind,1]  + mcmc.mu.nonstat.parallel[[1]]$samples[ind,2] * temperature_proj[165],
                    scale = mcmc.mu.nonstat.parallel[[1]]$samples[ind,3], 
                    shape = mcmc.mu.nonstat.parallel[[1]]$samples[ind,4],
                    type = c("GEV"))
non.stat.2065 <-  devd(x.lsl, loc = mcmc.mu.nonstat.parallel[[1]]$samples[ind,1]  + mcmc.mu.nonstat.parallel[[1]]$samples[ind,2] * temperature_proj[215],
                       scale = mcmc.mu.nonstat.parallel[[1]]$samples[ind,3], 
                       shape = mcmc.mu.nonstat.parallel[[1]]$samples[ind,4],
                       type = c("GEV"))
non.stat.2095 <-  devd(x.lsl, loc = mcmc.mu.nonstat.parallel[[1]]$samples[ind,1]  + mcmc.mu.nonstat.parallel[[1]]$samples[ind,2] * temperature_proj[245],
                       scale = mcmc.mu.nonstat.parallel[[1]]$samples[ind,3], 
                       shape = mcmc.mu.nonstat.parallel[[1]]$samples[ind,4],
                       type = c("GEV"))

plot(non.stat.2015, type='l', col='black', lwd = 2, xlab = 'Sea Level Height [mm]', ylab = 'Probability')
lines(non.stat.2065, col='red', lwd = 2)
lines(non.stat.2095, col='blue', lwd = 2)
title('Distribution of Maxmimum Sea Level Heights with Location Parameter Non-Stationary')
legend(2500,.0017, legend= c('2015', '2065', '2095'), lty=c(1,1,1), lwd=c(2.5, 2.5,2.5), col=c('black', 'red','blue'))
#-----------------------------------------------


#Non-Stationary Ex Plot - All Params Non Stat--
#load('DEOptim.allnonstat2.RData')
dev.off()
x.lsl2 <- seq(from= 500, to =1500, by = 1)
par(mfrow = c(1,1))
ind <- best.params.mu.sigma.xi.nonstat 
all.non.stat.2015 <- devd(x.lsl2, loc = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 1]  + mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 2] * temperature_proj[165],
                      scale = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,3] + mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,4]*temperature_proj[165], 
                      shape = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,5] +  mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 6]*temperature_proj[165],
                      type = c("GEV"))
all.non.stat.2065 <-  devd(x.lsl2, loc = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 1]  + mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 2] * temperature_proj[215],
                           scale = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,3] + mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,4]*temperature_proj[215], 
                           shape = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,5] +  mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 6]*temperature_proj[215],
                           type = c("GEV"))
all.non.stat.2095 <-  devd(x.lsl2, loc = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 1]  + mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 2] * temperature_proj[245],
                           scale = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,3] + mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,4]*temperature_proj[245], 
                           shape = mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind,5] +  mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[ind, 6]*temperature_proj[245],
                           type = c("GEV"))

plot(all.non.stat.2095, type='l', col='blue', lwd = 2,  xlab = 'Sea Level Height [mm]', ylab = 'Probability')
lines(all.non.stat.2065, col='red', lwd = 2)
lines(all.non.stat.2015, col='black', lwd = 2)
title('Distribution of Maxmimum Sea Level Heights with All Parameters Non-Stationary')
legend(800,.15, legend= c('2015', '2065', '2095'), lty=c(1,1,1), lwd=c(2.5, 2.5,2.5), col=c('black', 'red','blue'))
#-----------------------------------------------

#load('mcmc.parallel.longer.iter.RData')

#MCMC Plot, All Params Stationary --------------
dev.off()
par(mfrow = c(3,1))
plot(mcmc.stationary.parallel[[1]]$samples[,1], 
     type='l',xlab = 'Iteration', main = 'Location, Scale, and Shape of Stationary MCMC Chains', ylab = 'Location')
lines(mcmc.stationary.parallel[[2]]$samples[,1], col = 'red')
lines(mcmc.stationary.parallel[[3]]$samples[,1], col = 'blue')
plot(mcmc.stationary.parallel[[1]]$samples[,2], 
     type='l', xlab = 'Iteration',ylab = 'Scale')
lines(mcmc.stationary.parallel[[2]]$samples[,2], col = 'red')
lines(mcmc.stationary.parallel[[3]]$samples[,2], col = 'blue')
plot(mcmc.stationary.parallel[[1]]$samples[,3], 
     type='l', xlab = 'Iteration', ylab = 'Shape')
lines(mcmc.stationary.parallel[[2]]$samples[,3], col = 'red')
lines(mcmc.stationary.parallel[[3]]$samples[,3], col = 'blue')
#-----------------------------------------------

#MCMC Plot, Only Mu Non Stationary--------------
#par(mfrow = c(4,1))
dev.off()
b <- rbind(c(1),c(2), c(3), c(4))
layout(b)
par(mar = c(0,1,0.5,1))
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,1], 
     type='l',xlab = 'Iteration', main = 'Scale, and Shape of Stationary, Location Non Stat MCMC Chains', ylab = 'Location')
lines(mcmc.mu.nonstat.parallel[[2]]$samples[,1], col = 'red')
lines(mcmc.mu.nonstat.parallel[[3]]$samples[,1], col = 'blue')
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,2], 
     type='l', xlab = 'Iteration',ylab = 'Location * Temperature')
lines(mcmc.mu.nonstat.parallel[[2]]$samples[,2], col = 'red')
lines(mcmc.mu.nonstat.parallel[[3]]$samples[,2], col = 'blue')
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,3], 
     type='l', xlab = 'Iteration', ylab = 'Scale')
lines(mcmc.mu.nonstat.parallel[[2]]$samples[,3], col = 'red')
lines(mcmc.mu.nonstat.parallel[[3]]$samples[,3], col = 'blue')
plot(mcmc.mu.nonstat.parallel[[1]]$samples[,4], 
     type='l', xlab = 'Iteration', ylab = 'Shape')
lines(mcmc.mu.nonstat.parallel[[2]]$samples[,4], col = 'red')
lines(mcmc.mu.nonstat.parallel[[3]]$samples[,4], col = 'blue')


#-----------------------------------------------

#Beginnings of Mega Plot of All Prior and Posterior Params----------
#load('mcmc.params.RData')
#load('mcmc.allnonstat.for.real.RData')

#load('mcmc.rerun2.RData')

#load('mcmc.rerun3.RData')
load('post.burnin.iter1e6.RData')
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

#load('post.burnin.iter1e5.RData')
#stationary case 
#load('mcmc.stat.post.burnin.RData')

plot(density(mcmc.stationary.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.stationary.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = 0, to = 400, by = 1),min = 0, max = 400), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.stationary.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot.new()

#non stationary mu 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = -500, to = 500, by = 1), min = -500, max = 500), type = 'l', col='blue', lty=2) #why the negs not showing up?
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = 0, to = 400, by = 1), min = 0, max = 400), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.mu.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot.new()
#plot(dunif(seq(from = -500, to = 500, by = 1), min = -500, max = 500), type = 'l', col='blue', lty=2)
#non stationary sigma 
par(mar = c(0,0,0,0))
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = -300, to = 300, by = 1), min = -300, max = 300), type = 'l', col='blue', lty=2)
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = -100, to = 200, by = 1), min = -100, max = 200), type = 'l', col='blue', lty=2)
plot(density(mcmc.sigma.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot.new()

#non stationary xi 
par(mar = c(0,0,0,0))
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = 0, to = 1000, by = 1), min = 0, max = 1000), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot(density(mcmc.xi.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)

#non stationary mu and sigma 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from =-500, to = 500, by = 1), min = -500, max = 500), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = 0, to = 100, by = 1), min = 0, max = 100), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -100, to = 100, by = 1), min = -100, max = 100), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.mu.sigma.nonstat.parallel.post.burnin.1[,5]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)


#non stationary mu and xi 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = -100, to = 100, by = 1), min = -100, max = 100), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = 0, to = 1000, by = 1), min = 0, max = 1000), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.xi.nonstat.parallel.post.burnin.1[,5]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)

#non stationary sigma and xi 
par(mar = c(0,0,0,0))
plot(density(mcmc.sigma.xi.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot.new()
plot(density(mcmc.sigma.xi.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = -100, to = 100, by = 1), min = -100, max = 100), type = 'l', col='blue', lty=2)
plot(density(mcmc.sigma.xi.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = -150, to = 150, by = 1), min = -150, max = 150), type = 'l', col='blue', lty=2)
plot(density(mcmc.sigma.xi.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot(density(mcmc.sigma.xi.nonstat.parallel.post.burnin.1[,5]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)

#load('post.burnin.iter1e5.RData')
#all non stationary 
par(mar = c(0,0,0,0))
plot(density(mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1[,1]), col = 'blue')
lines(dunif(seq(from = 0, to = 3000, by = 1), min = 0, max = 3000), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1[,2]), col = 'blue')
lines(dunif(seq(from = -500, to = 500, by = 1), min = -500, max = 500), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1[,3]), col = 'blue')
lines(dunif(seq(from = -500, to = 500, by = 1), min = -500, max = 500), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1[,4]), col = 'blue')
lines(dunif(seq(from = -500, to = 500, by = 1), min = -500, max = 500), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1[,5]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
plot(density(mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1[,6]), col = 'blue')
lines(dunif(seq(from = -10, to = 10, by = 1), min = -10, max = 10), type = 'l', col='blue', lty=2)
#-------------------------------------------------------------------

#year comparison plots w params---------------------------------------------------------------
#load('10.year.gaps.RData')
par(mfrow = c(3,2))
par(mar = c(0,4,2,1))
#mu0
#new london
nl.years <- c(10,20,30,40,50,60,76)
nl.mu0 <- c(new.london.10$mu.sigma.xi[1],
            new.london.20$mu.sigma.xi[1],
            new.london.30$mu.sigma.xi[1],
            new.london.40$mu.sigma.xi[1],
            new.london.50$mu.sigma.xi[1],
            new.london.60$mu.sigma.xi[1],
            new.london.all$mu.sigma.xi[1])
plot(nl.years, nl.mu0, type = 'l', xlim = c(0,110), ylim = c(950,2500), xaxt = 'n', ylab = 'mu0')
points(10, new.london.10$mu.sigma.xi[1], cex = 1.5)
points(20, new.london.20$mu.sigma.xi[1], cex = 1.5)
points(30, new.london.30$mu.sigma.xi[1], cex = 1.5)
points(40, new.london.40$mu.sigma.xi[1], cex = 1.5)
points(50, new.london.50$mu.sigma.xi[1], cex = 1.5)
points(60, new.london.60$mu.sigma.xi[1], cex = 1.5)
points(76, new.london.all$mu.sigma.xi[1], cex = 1.5)

#boston 
bos.years <- c(10,20,30,40,50,60,70,80,94)
bos.mu0 <- c(boston.10$mu.sigma.xi[1], 
             boston.20$mu.sigma.xi[1], 
             boston.30$mu.sigma.xi[1], 
             boston.40$mu.sigma.xi[1] ,
             boston.50$mu.sigma.xi[1], 
             boston.60$mu.sigma.xi[1], 
             boston.70$mu.sigma.xi[1], 
             boston.80$mu.sigma.xi[1], 
             boston.all$mu.sigma.xi[1])
lines(bos.years, bos.mu0)
points(10,boston.10$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(20,boston.20$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(30,boston.30$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(40,boston.40$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(50, boston.50$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(60, boston.60$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(70,boston.70$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(80,boston.80$mu.sigma.xi[1], pch = 0, cex = 1.5)
points(94, boston.all$mu.sigma.xi[1], pch = 0, cex = 1.5)

#atlantic city 
al.city.years <- c(10,20,30,40,50,60,70,80,90, 103)
al.city.mu0 <- c(atlantic.city.10$mu.sigma.xi[1], 
                 atlantic.city.20$mu.sigma.xi[1],
                 atlantic.city.30$mu.sigma.xi[1],
                 atlantic.city.40$mu.sigma.xi[1],
                 atlantic.city.50$mu.sigma.xi[1],
                 atlantic.city.60$mu.sigma.xi[1],
                 atlantic.city.70$mu.sigma.xi[1],
                 atlantic.city.80$mu.sigma.xi[1],
                 atlantic.city.90$mu.sigma.xi[1],
                 atlantic.city.all$mu.sigma.xi[1])
lines(al.city.years, al.city.mu0)
points(10,atlantic.city.10$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(20,atlantic.city.20$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(30,atlantic.city.30$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(40,atlantic.city.40$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(50,atlantic.city.50$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(60,atlantic.city.60$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(70, atlantic.city.70$mu.sigma.xi[1],pch = 2, cex = 1.5)
points(80, atlantic.city.80$mu.sigma.xi[1],pch = 2, cex = 1.5)
points(90,atlantic.city.90$mu.sigma.xi[1], pch = 2, cex = 1.5)
points(103,atlantic.city.all$mu.sigma.xi[1], pch = 2, cex = 1.5)

#portland 
portland.years <- c(10,20,30,40,50,60,70,80,90, 105)
portland.mu0 <- c(portland.10$mu.sigma.xi[1], 
                  portland.20$mu.sigma.xi[1],
                  portland.30$mu.sigma.xi[1], 
                  portland.40$mu.sigma.xi[1],
                  portland.50$mu.sigma.xi[1],
                  portland.60$mu.sigma.xi[1],
                  portland.70$mu.sigma.xi[1],
                  portland.80$mu.sigma.xi[1], 
                  portland.90$mu.sigma.xi[1], 
                  portland.all$mu.sigma.xi[1])
lines(portland.years, portland.mu0)
points(10,portland.10$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(20,portland.20$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(30,portland.30$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(40,portland.40$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(50, portland.50$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(60,portland.60$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(70,portland.70$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(80,portland.80$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(90,portland.90$mu.sigma.xi[1], pch = 3, cex = 1.5)
points(105,portland.all$mu.sigma.xi[1], pch = 3, cex = 1.5)

#legend(30, 2000, legend= c('30 year', '45 year', '60 year', 'all year'), 
#lty =c(1,1,1,1), col=c('black', 'red','blue', 'green'), cex = .5)
legend(60,2000, legend= c('New London', 'Boston', 'Atlantic City', 'Portland'), 
       pch =c(1,0,2,3), 
       col=c('black', 'black','black', 'black'), cex = .5)


#mu1
#new london
nl.mu1 <- c(new.london.10$mu.sigma.xi[2],
            new.london.20$mu.sigma.xi[2],
            new.london.30$mu.sigma.xi[2],
            new.london.40$mu.sigma.xi[2],
            new.london.50$mu.sigma.xi[2],
            new.london.60$mu.sigma.xi[2],
            new.london.all$mu.sigma.xi[2])
plot(nl.years, nl.mu1, type = 'l',xlim = c(0, 110), ylim = c(-150,150), xaxt = 'n', ylab = 'mu1')
points(10, new.london.10$mu.sigma.xi[2], cex = 1.5)
points(20, new.london.20$mu.sigma.xi[2], cex = 1.5)
points(30, new.london.30$mu.sigma.xi[2], cex = 1.5)
points(40, new.london.40$mu.sigma.xi[2], cex = 1.5)
points(50, new.london.50$mu.sigma.xi[2], cex = 1.5)
points(60, new.london.60$mu.sigma.xi[2], cex = 1.5)
points(76, new.london.all$mu.sigma.xi[2], cex = 1.5)

#boston 
bos.mu1 <- c(boston.10$mu.sigma.xi[2], 
             boston.20$mu.sigma.xi[2], 
             boston.30$mu.sigma.xi[2], 
             boston.40$mu.sigma.xi[2] ,
             boston.50$mu.sigma.xi[2], 
             boston.60$mu.sigma.xi[2], 
             boston.70$mu.sigma.xi[2], 
             boston.80$mu.sigma.xi[2], 
             boston.all$mu.sigma.xi[2])
lines(bos.years, bos.mu1)
points(10,boston.10$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(20,boston.20$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(30,boston.30$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(40,boston.40$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(50, boston.50$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(60, boston.60$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(70,boston.70$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(80,boston.80$mu.sigma.xi[2], pch = 0, cex = 1.5)
points(94, boston.all$mu.sigma.xi[2], pch = 0, cex = 1.5)

#atlantic city 
al.city.mu1 <- c(atlantic.city.10$mu.sigma.xi[2], 
                 atlantic.city.20$mu.sigma.xi[2],
                 atlantic.city.30$mu.sigma.xi[2],
                 atlantic.city.40$mu.sigma.xi[2],
                 atlantic.city.50$mu.sigma.xi[2],
                 atlantic.city.60$mu.sigma.xi[2],
                 atlantic.city.70$mu.sigma.xi[2],
                 atlantic.city.80$mu.sigma.xi[2],
                 atlantic.city.90$mu.sigma.xi[2],
                 atlantic.city.all$mu.sigma.xi[2])
lines(al.city.years, al.city.mu1)
points(10,atlantic.city.10$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(20,atlantic.city.20$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(30,atlantic.city.30$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(40,atlantic.city.40$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(50,atlantic.city.50$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(60,atlantic.city.60$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(70, atlantic.city.70$mu.sigma.xi[2],pch = 2, cex = 1.5)
points(80, atlantic.city.80$mu.sigma.xi[2],pch = 2, cex = 1.5)
points(90,atlantic.city.90$mu.sigma.xi[2], pch = 2, cex = 1.5)
points(103,atlantic.city.all$mu.sigma.xi[2], pch = 2, cex = 1.5)

#portland 
portland.mu1 <- c(portland.10$mu.sigma.xi[2], 
                  portland.20$mu.sigma.xi[2],
                  portland.30$mu.sigma.xi[2], 
                  portland.40$mu.sigma.xi[2],
                  portland.50$mu.sigma.xi[2],
                  portland.60$mu.sigma.xi[2],
                  portland.70$mu.sigma.xi[2],
                  portland.80$mu.sigma.xi[2], 
                  portland.90$mu.sigma.xi[2], 
                  portland.all$mu.sigma.xi[2])
lines(portland.years, portland.mu1)
points(10,portland.10$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(20,portland.20$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(30,portland.30$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(40,portland.40$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(50, portland.50$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(60,portland.60$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(70,portland.70$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(80,portland.80$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(90,portland.90$mu.sigma.xi[2], pch = 3, cex = 1.5)
points(105,portland.all$mu.sigma.xi[2], pch = 3, cex = 1.5)

#sigma0 
#new london
par(mar = c(0,4,0,1))
nl.sigma0 <- c(new.london.10$mu.sigma.xi[3],
               new.london.20$mu.sigma.xi[3],
               new.london.30$mu.sigma.xi[3],
               new.london.40$mu.sigma.xi[3],
               new.london.50$mu.sigma.xi[3],
               new.london.60$mu.sigma.xi[3],
               new.london.all$mu.sigma.xi[3])
plot(nl.years, nl.sigma0, xlim = c(0, 110), ylim = c(1.75,6), type = 'l', xaxt = 'n', ylab = 'sigma0')
points(10, new.london.10$mu.sigma.xi[3], cex = 1.5)
points(20, new.london.20$mu.sigma.xi[3], cex = 1.5)
points(30, new.london.30$mu.sigma.xi[3], cex = 1.5)
points(40, new.london.40$mu.sigma.xi[3], cex = 1.5)
points(50, new.london.50$mu.sigma.xi[3], cex = 1.5)
points(60, new.london.60$mu.sigma.xi[3], cex = 1.5)
points(76, new.london.all$mu.sigma.xi[3], cex = 1.5)

#boston 
bos.sigma0 <- c(boston.10$mu.sigma.xi[3], 
                boston.20$mu.sigma.xi[3], 
                boston.30$mu.sigma.xi[3], 
                boston.40$mu.sigma.xi[3] ,
                boston.50$mu.sigma.xi[3], 
                boston.60$mu.sigma.xi[3], 
                boston.70$mu.sigma.xi[3], 
                boston.80$mu.sigma.xi[3], 
                boston.all$mu.sigma.xi[3])
lines(bos.years, bos.sigma0)
points(10,boston.10$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(20,boston.20$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(30,boston.30$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(40,boston.40$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(50, boston.50$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(60, boston.60$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(70,boston.70$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(80,boston.80$mu.sigma.xi[3], pch = 0, cex = 1.5)
points(94, boston.all$mu.sigma.xi[3], pch = 0, cex = 1.5)

#atlantic city 
al.city.sigma0 <- c(atlantic.city.10$mu.sigma.xi[3], 
                    atlantic.city.20$mu.sigma.xi[3],
                    atlantic.city.30$mu.sigma.xi[3],
                    atlantic.city.40$mu.sigma.xi[3],
                    atlantic.city.50$mu.sigma.xi[3],
                    atlantic.city.60$mu.sigma.xi[3],
                    atlantic.city.70$mu.sigma.xi[3],
                    atlantic.city.80$mu.sigma.xi[3],
                    atlantic.city.90$mu.sigma.xi[3],
                    atlantic.city.all$mu.sigma.xi[3])
lines(al.city.years, al.city.sigma0)
points(10,atlantic.city.10$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(20,atlantic.city.20$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(30,atlantic.city.30$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(40,atlantic.city.40$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(50,atlantic.city.50$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(60,atlantic.city.60$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(70, atlantic.city.70$mu.sigma.xi[3],pch = 2, cex = 1.5)
points(80, atlantic.city.80$mu.sigma.xi[3],pch = 2, cex = 1.5)
points(90,atlantic.city.90$mu.sigma.xi[3], pch = 2, cex = 1.5)
points(103,atlantic.city.all$mu.sigma.xi[3], pch = 2, cex = 1.5)

#portland 
portland.sigma0 <- c(portland.10$mu.sigma.xi[3], 
                     portland.20$mu.sigma.xi[3],
                     portland.30$mu.sigma.xi[3], 
                     portland.40$mu.sigma.xi[3],
                     portland.50$mu.sigma.xi[3],
                     portland.60$mu.sigma.xi[3],
                     portland.70$mu.sigma.xi[3],
                     portland.80$mu.sigma.xi[3], 
                     portland.90$mu.sigma.xi[3], 
                     portland.all$mu.sigma.xi[3])
lines(portland.years, portland.sigma0)
points(10,portland.10$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(20,portland.20$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(30,portland.30$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(40,portland.40$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(50, portland.50$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(60,portland.60$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(70,portland.70$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(80,portland.80$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(90,portland.90$mu.sigma.xi[3], pch = 3, cex = 1.5)
points(105,portland.all$mu.sigma.xi[3], pch = 3, cex = 1.5)

#sigma1
#new london
nl.sigma1 <- c(new.london.10$mu.sigma.xi[4],
               new.london.20$mu.sigma.xi[4],
               new.london.30$mu.sigma.xi[4],
               new.london.40$mu.sigma.xi[4],
               new.london.50$mu.sigma.xi[4],
               new.london.60$mu.sigma.xi[4],
               new.london.all$mu.sigma.xi[4])
plot(nl.years, nl.sigma1,xlim = c(0, 110), ylim = c(-3,5), type = 'l' , xaxt = 'n', ylab = 'sigma1')
points(10, new.london.10$mu.sigma.xi[4], cex = 1.5)
points(20, new.london.20$mu.sigma.xi[4], cex = 1.5)
points(30, new.london.30$mu.sigma.xi[4], cex = 1.5)
points(40, new.london.40$mu.sigma.xi[4], cex = 1.5)
points(50, new.london.50$mu.sigma.xi[4], cex = 1.5)
points(60, new.london.60$mu.sigma.xi[4], cex = 1.5)
points(76, new.london.all$mu.sigma.xi[4], cex = 1.5)

#boston 
bos.sigma1 <- c(boston.10$mu.sigma.xi[4], 
                boston.20$mu.sigma.xi[4], 
                boston.30$mu.sigma.xi[4], 
                boston.40$mu.sigma.xi[4] ,
                boston.50$mu.sigma.xi[4], 
                boston.60$mu.sigma.xi[4], 
                boston.70$mu.sigma.xi[4], 
                boston.80$mu.sigma.xi[4], 
                boston.all$mu.sigma.xi[4])
lines(bos.years, bos.sigma1 )

points(10,boston.10$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(20,boston.20$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(30,boston.30$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(40,boston.40$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(50, boston.50$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(60, boston.60$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(70,boston.70$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(80,boston.80$mu.sigma.xi[4], pch = 0, cex = 1.5)
points(94, boston.all$mu.sigma.xi[4], pch = 0, cex = 1.5)

#atlantic city 
al.city.sigma1 <- c(atlantic.city.10$mu.sigma.xi[4], 
                    atlantic.city.20$mu.sigma.xi[4],
                    atlantic.city.30$mu.sigma.xi[4],
                    atlantic.city.40$mu.sigma.xi[4],
                    atlantic.city.50$mu.sigma.xi[4],
                    atlantic.city.60$mu.sigma.xi[4],
                    atlantic.city.70$mu.sigma.xi[4],
                    atlantic.city.80$mu.sigma.xi[4],
                    atlantic.city.90$mu.sigma.xi[4],
                    atlantic.city.all$mu.sigma.xi[4])
lines(al.city.years, al.city.sigma1)
points(10,atlantic.city.10$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(20,atlantic.city.20$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(30,atlantic.city.30$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(40,atlantic.city.40$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(50,atlantic.city.50$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(60,atlantic.city.60$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(70, atlantic.city.70$mu.sigma.xi[4],pch = 2, cex = 1.5)
points(80, atlantic.city.80$mu.sigma.xi[4],pch = 2, cex = 1.5)
points(90,atlantic.city.90$mu.sigma.xi[4], pch = 2, cex = 1.5)
points(103,atlantic.city.all$mu.sigma.xi[4], pch = 2, cex = 1.5)

#portland 
portland.sigma1 <- c(portland.10$mu.sigma.xi[4], 
                     portland.20$mu.sigma.xi[4],
                     portland.30$mu.sigma.xi[4], 
                     portland.40$mu.sigma.xi[4],
                     portland.50$mu.sigma.xi[4],
                     portland.60$mu.sigma.xi[4],
                     portland.70$mu.sigma.xi[4],
                     portland.80$mu.sigma.xi[4], 
                     portland.90$mu.sigma.xi[4], 
                     portland.all$mu.sigma.xi[4])
lines(portland.years, portland.sigma1)
points(10,portland.10$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(20,portland.20$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(30,portland.30$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(40,portland.40$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(50, portland.50$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(60,portland.60$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(70,portland.70$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(80,portland.80$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(90,portland.90$mu.sigma.xi[4], pch = 3, cex = 1.5)
points(105,portland.all$mu.sigma.xi[4], pch = 3, cex = 1.5)

#xi0
#new london
par(mar = c(4,4,0,1))
nl.xi0 <- c(new.london.10$mu.sigma.xi[5],
            new.london.20$mu.sigma.xi[5],
            new.london.30$mu.sigma.xi[5],
            new.london.40$mu.sigma.xi[5],
            new.london.50$mu.sigma.xi[5],
            new.london.60$mu.sigma.xi[5],
            new.london.all$mu.sigma.xi[5])
plot(nl.years, nl.xi0,xlim = c(0, 110), ylim = c(-1.5,1.5), type = 'l', ylab = 'xi0', xlab = 'Number of Years of Data')
points(10, new.london.10$mu.sigma.xi[5], cex = 1.5)
points(20, new.london.20$mu.sigma.xi[5], cex = 1.5)
points(30, new.london.30$mu.sigma.xi[5], cex = 1.5)
points(40, new.london.40$mu.sigma.xi[5], cex = 1.5)
points(50, new.london.50$mu.sigma.xi[5], cex = 1.5)
points(60, new.london.60$mu.sigma.xi[5], cex = 1.5)
points(76, new.london.all$mu.sigma.xi[5], cex = 1.5)

#boston 
bos.xi0 <- c(boston.10$mu.sigma.xi[5], 
             boston.20$mu.sigma.xi[5], 
             boston.30$mu.sigma.xi[5], 
             boston.40$mu.sigma.xi[5] ,
             boston.50$mu.sigma.xi[5], 
             boston.60$mu.sigma.xi[5], 
             boston.70$mu.sigma.xi[5], 
             boston.80$mu.sigma.xi[5], 
             boston.all$mu.sigma.xi[5])
lines(bos.years, bos.xi0)
points(10,boston.10$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(20,boston.20$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(30,boston.30$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(40,boston.40$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(50, boston.50$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(60, boston.60$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(70,boston.70$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(80,boston.80$mu.sigma.xi[5], pch = 0, cex = 1.5)
points(94, boston.all$mu.sigma.xi[5], pch = 0, cex = 1.5)

#atlantic city 
al.city.xi0 <- c(atlantic.city.10$mu.sigma.xi[5], 
                 atlantic.city.20$mu.sigma.xi[5],
                 atlantic.city.30$mu.sigma.xi[5],
                 atlantic.city.40$mu.sigma.xi[5],
                 atlantic.city.50$mu.sigma.xi[5],
                 atlantic.city.60$mu.sigma.xi[5],
                 atlantic.city.70$mu.sigma.xi[5],
                 atlantic.city.80$mu.sigma.xi[5],
                 atlantic.city.90$mu.sigma.xi[5],
                 atlantic.city.all$mu.sigma.xi[5])
lines(al.city.years, al.city.xi0)
points(10,atlantic.city.10$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(20,atlantic.city.20$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(30,atlantic.city.30$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(40,atlantic.city.40$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(50,atlantic.city.50$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(60,atlantic.city.60$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(70, atlantic.city.70$mu.sigma.xi[5],pch = 2, cex = 1.5)
points(80, atlantic.city.80$mu.sigma.xi[5],pch = 2, cex = 1.5)
points(90,atlantic.city.90$mu.sigma.xi[5], pch = 2, cex = 1.5)
points(103,atlantic.city.all$mu.sigma.xi[5], pch = 2, cex = 1.5)

#portland 
portland.xi0<- c(portland.10$mu.sigma.xi[5], 
                 portland.20$mu.sigma.xi[5],
                 portland.30$mu.sigma.xi[5], 
                 portland.40$mu.sigma.xi[5],
                 portland.50$mu.sigma.xi[5],
                 portland.60$mu.sigma.xi[5],
                 portland.70$mu.sigma.xi[5],
                 portland.80$mu.sigma.xi[5], 
                 portland.90$mu.sigma.xi[5], 
                 portland.all$mu.sigma.xi[5])
lines(portland.years, portland.xi0)
points(10,portland.10$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(20,portland.20$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(30,portland.30$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(40,portland.40$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(50, portland.50$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(60,portland.60$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(70,portland.70$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(80,portland.80$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(90,portland.90$mu.sigma.xi[5], pch = 3, cex = 1.5)
points(105,portland.all$mu.sigma.xi[5], pch = 3, cex = 1.5)

#xi1
#new london
nl.xi1 <- c(new.london.10$mu.sigma.xi[6],
            new.london.20$mu.sigma.xi[6],
            new.london.30$mu.sigma.xi[6],
            new.london.40$mu.sigma.xi[6],
            new.london.50$mu.sigma.xi[6],
            new.london.60$mu.sigma.xi[6],
            new.london.all$mu.sigma.xi[6])
plot(nl.years, nl.xi1, xlim = c(0, 110), ylim = c(-1.5,1.5), type = 'l', ylab = 'xi1', xlab = 'Number of Years of Data')
points(10, new.london.10$mu.sigma.xi[6], cex = 1.5)
points(20, new.london.20$mu.sigma.xi[6], cex = 1.5)
points(30, new.london.30$mu.sigma.xi[6], cex = 1.5)
points(40, new.london.40$mu.sigma.xi[6], cex = 1.5)
points(50, new.london.50$mu.sigma.xi[6], cex = 1.5)
points(60, new.london.60$mu.sigma.xi[6], cex = 1.5)
points(76, new.london.all$mu.sigma.xi[6], cex = 1.5)

#boston 
bos.xi1 <- c(boston.10$mu.sigma.xi[6], 
             boston.20$mu.sigma.xi[6], 
             boston.30$mu.sigma.xi[6], 
             boston.40$mu.sigma.xi[6] ,
             boston.50$mu.sigma.xi[6], 
             boston.60$mu.sigma.xi[6], 
             boston.70$mu.sigma.xi[6], 
             boston.80$mu.sigma.xi[6], 
             boston.all$mu.sigma.xi[6])
lines(bos.years, bos.xi1)
points(10,boston.10$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(20,boston.20$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(30,boston.30$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(40,boston.40$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(50, boston.50$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(60, boston.60$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(70,boston.70$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(80,boston.80$mu.sigma.xi[6], pch = 0, cex = 1.5)
points(94, boston.all$mu.sigma.xi[6], pch = 0, cex = 1.5)

#atlantic city 
al.city.xi1 <- c(atlantic.city.10$mu.sigma.xi[6], 
                 atlantic.city.20$mu.sigma.xi[6],
                 atlantic.city.30$mu.sigma.xi[6],
                 atlantic.city.40$mu.sigma.xi[6],
                 atlantic.city.50$mu.sigma.xi[6],
                 atlantic.city.60$mu.sigma.xi[6],
                 atlantic.city.70$mu.sigma.xi[6],
                 atlantic.city.80$mu.sigma.xi[6],
                 atlantic.city.90$mu.sigma.xi[6],
                 atlantic.city.all$mu.sigma.xi[6])
lines(al.city.years, al.city.xi1)
points(10,atlantic.city.10$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(20,atlantic.city.20$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(30,atlantic.city.30$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(40,atlantic.city.40$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(50,atlantic.city.50$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(60,atlantic.city.60$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(70, atlantic.city.70$mu.sigma.xi[6],pch = 2, cex = 1.5)
points(80, atlantic.city.80$mu.sigma.xi[6],pch = 2, cex = 1.5)
points(90,atlantic.city.90$mu.sigma.xi[6], pch = 2, cex = 1.5)
points(103,atlantic.city.all$mu.sigma.xi[6], pch = 2, cex = 1.5)

#portland 
portland.xi1 <- c(portland.10$mu.sigma.xi[6], 
                  portland.20$mu.sigma.xi[6],
                  portland.30$mu.sigma.xi[6], 
                  portland.40$mu.sigma.xi[6],
                  portland.50$mu.sigma.xi[6],
                  portland.60$mu.sigma.xi[6],
                  portland.70$mu.sigma.xi[6],
                  portland.80$mu.sigma.xi[6], 
                  portland.90$mu.sigma.xi[6], 
                  portland.all$mu.sigma.xi[6])
lines(portland.years, portland.xi1)
points(10,portland.10$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(20,portland.20$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(30,portland.30$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(40,portland.40$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(50, portland.50$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(60,portland.60$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(70,portland.70$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(80,portland.80$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(90,portland.90$mu.sigma.xi[6], pch = 3, cex = 1.5)
points(105,portland.all$mu.sigma.xi[6], pch = 3, cex = 1.5)
#-----------------------------------------------------------------------------