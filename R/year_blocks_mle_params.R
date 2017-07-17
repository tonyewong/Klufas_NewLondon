#working with different year blocks for estimating params 

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data(1939, 2014)

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
tide.data <- read.tide.data()

#setwd('~/codes/Klufas_NewLondon/R/')
#source('gev_nonstationary_MCMC.R')
setwd('~/codes/Klufas_NewLondon/R/')
source('optimization_sf.R')
#library(ismev)
#library(zoo)

lower.bound <- list(c(0,0,-5),
                    c(0,-100, 0, -5),
                    c(0, 0, 0, -200),
                    c(0, 0, -5, -5), 
                    c(0,-100, 0, 0, -5),
                    c(0, -100, 0, -1, -1),
                    c(0, -100, 0, -1, -10),
                    c(0,-100, 0, -100, -1, -1)
)
upper.bound <- list(c(3000,3000,5),
                    c(3000,1000, 1000, 5),
                    c(3000, 1000, 1000, 1000),
                    c(3000, 1000 , 100, 100),
                    c(3000,100, 1000, 10 , 5),
                    c(3000, 1000, 10 , 1, 1),
                    c(3000, 100, 1000, 10, 10),
                    c(3000,100, 1000, 10 , 1, 1)
)

get.prior.estimates <- function(max.data.clipped, city.temps){
  #for (i in 1:2){
  
  sublist <- vector('list', 8)
  names(sublist) <- other.p.names
  for (j in 1: length(all.p.names)){
    #print(j)
    optim.like <- DEoptim(neg.log.like.calc,
                          lower=lower.bound[[j]], 
                          upper=upper.bound[[j]], 
                          data = max.data.clipped, 
                          temps = city.temps, 
                          parnames = all.p.names[[j]])
    sublist[[other.p.names[j]]] <- optim.like$optim$bestmem
  }
  #optim.gev.fit$city.names[i] <- sublist
  return(sublist)
}
#New London Plots-----------------------------------------------------------------------------------
new.london.30 <- get.prior.estimates(tide.data$max[46:76], temps$values[46:76])
new.london.45 <- get.prior.estimates(tide.data$max[31:76], temps$values[31:76])
new.london.60 <- get.prior.estimates(tide.data$max[16:76], temps$values[16:76])
new.london.all <- get.prior.estimates(tide.data$max, temps$values)

par(mfrow = c(1,3))
plot(30, new.london.30$stat[1], xlim = c(25, 80))
points(45, new.london.45$stat[1], col = 'red')
points(60, new.london.60$stat[1], col = 'blue')
points(76, new.london.all$stat[1], col = 'green')
legend(30, 1400, legend= c('30 year', '45 year', '60 year', 'all year'), lty =c(1,1,1,1), col=c('black', 'red','blue', 'green'))

plot(30, new.london.30$stat[2], xlim = c(25, 80))
points(45, new.london.45$stat[2], col = 'red')
points(60, new.london.60$stat[2], col = 'blue')
points(76, new.london.all$stat[2], col = 'green')

plot(30, new.london.30$stat[3], ylim = c(0,.2), xlim = c(25, 80))
points(45, new.london.45$stat[3], col = 'red')
points(60, new.london.60$stat[3], col = 'blue')
points(76, new.london.all$stat[3], col = 'green')

par(mfrow = c(2,2))
plot(30, new.london.30$mu[1], xlim = c(25, 80))
points(45, new.london.45$mu[1], col = 'red')
points(60,new.london.60$mu[1], col = 'blue')
points(76,new.london.all$mu[1], col = 'green')
legend(25, 1000, legend= c('30 year', '45 year', '60 year', 'all year'), 
       lty =c(1,1,1,1), col=c('black', 'red','blue', 'green'), cex = .5)

plot(30,new.london.30$mu[2], ylim = c(-120, 0), xlim = c(25, 80))
points(45,new.london.45$mu[2], col = 'red')
points(60,new.london.60$mu[2], col = 'blue')
points(76,new.london.all$mu[2], col = 'green')

plot(30,new.london.30$mu[3], ylim = c(100, 180), xlim = c(25, 80))
points(45,new.london.45$mu[3], col = 'red')
points(60,new.london.60$mu[3], col = 'blue')
points(76,new.london.all$mu[3], col = 'green')

plot(30,new.london.30$mu[4], ylim = c(-.5,.5), xlim = c(25, 80))
points(45,new.london.45$mu[4], col = 'red')
points(60,new.london.60$mu[4], col = 'blue')
points(76,new.london.all$mu[4], col = 'green')

par(mfrow = c(2,3))

plot(30, new.london.30$mu.sigma.xi[1], xlim = c(25, 80))
points(45, new.london.30$mu.sigma.xi[1], col = 'red')
points(60,new.london.30$mu.sigma.xi[1], col = 'blue')
points(76,new.london.30$mu.sigma.xi[1], col = 'green')
legend(25, 1000, legend= c('30 year', '45 year', '60 year', 'all year'), 
       lty =c(1,1,1,1), col=c('black', 'red','blue', 'green'), cex = .75)

plot(30,new.london.30$mu.sigma.xi[2], ylim = c(-120, 0), xlim = c(25, 80))
points(45,new.london.30$mu.sigma.xi[2], col = 'red')
points(60,new.london.30$mu.sigma.xi[2], col = 'blue')
points(76,new.london.30$mu.sigma.xi[2], col = 'green')

plot(30,new.london.30$mu.sigma.xi[3], ylim = c(0, 10), xlim = c(25, 80))
points(45,new.london.30$mu.sigma.xi[3], col = 'red')
points(60,new.london.30$mu.sigma.xi[3], col = 'blue')
points(76,new.london.30$mu.sigma.xi[3], col = 'green')

plot(30,new.london.30$mu.sigma.xi[4], ylim = c(-5,10), xlim = c(25, 80))
points(45,new.london.30$mu.sigma.xi[4], col = 'red')
points(60,new.london.30$mu.sigma.xi[4], col = 'blue')
points(76,new.london.30$mu.sigma.xi[4], col = 'green')

plot(30,new.london.30$mu.sigma.xi[5], ylim = c(0,5), xlim = c(25, 80))
points(45,new.london.30$mu.sigma.xi[5], col = 'red')
points(60,new.london.30$mu.sigma.xi[5], col = 'blue')
points(76,new.london.30$mu.sigma.xi[5], col = 'green')

plot(30,new.london.30$mu.sigma.xi[6], ylim = c(-5,5), xlim = c(25, 80))
points(45,new.london.30$mu.sigma.xi[6], col = 'red')
points(60,new.london.30$mu.sigma.xi[6], col = 'blue')
points(76,new.london.30$mu.sigma.xi[6], col = 'green')
#-----------------------------------------------------------------------------------

#load('east.coast.RData')
new.london.10 <- get.prior.estimates(tide.data$max[66:76], temps$values[66:76])
new.london.20 <- get.prior.estimates(tide.data$max[56:76], temps$values[56:76])
new.london.30 <- get.prior.estimates(tide.data$max[46:76], temps$values[46:76])
new.london.40 <- get.prior.estimates(tide.data$max[36:76], temps$values[36:76])
new.london.50 <- get.prior.estimates(tide.data$max[26:76], temps$values[26:76])
new.london.60 <- get.prior.estimates(tide.data$max[16:76], temps$values[16:76])
new.london.all <- get.prior.estimates(tide.data$max, temps$values)

boston.10 <- get.prior.estimates(boston$max[84:94], boston.temps$values[84:94])
boston.20 <- get.prior.estimates(boston$max[74:94], boston.temps$values[74:94])
boston.30 <- get.prior.estimates(boston$max[64:94], boston.temps$values[64:94])
boston.40 <- get.prior.estimates(boston$max[54:94], boston.temps$values[54:94])
boston.50 <- get.prior.estimates(boston$max[44:94], boston.temps$values[44:94])
boston.60 <- get.prior.estimates(boston$max[34:94], boston.temps$values[34:94])
boston.70 <- get.prior.estimates(boston$max[24:94], boston.temps$values[24:94])
boston.80 <- get.prior.estimates(boston$max[14:94], boston.temps$values[14:94])
boston.all <- get.prior.estimates(boston$max, boston.temps$values)

atlantic.city.10 <-  get.prior.estimates(atlantic.city$max[93:103], atlantic.city.temps$values[93:103])
atlantic.city.20 <-  get.prior.estimates(atlantic.city$max[83:103], atlantic.city.temps$values[83:103])
atlantic.city.30 <-  get.prior.estimates(atlantic.city$max[73:103], atlantic.city.temps$values[73:103])
atlantic.city.40 <-  get.prior.estimates(atlantic.city$max[63:103], atlantic.city.temps$values[63:103])
atlantic.city.50 <-  get.prior.estimates(atlantic.city$max[53:103], atlantic.city.temps$values[53:103])
atlantic.city.60 <-  get.prior.estimates(atlantic.city$max[43:103], atlantic.city.temps$values[43:103])
atlantic.city.70 <-  get.prior.estimates(atlantic.city$max[33:103], atlantic.city.temps$values[33:103])
atlantic.city.80 <-  get.prior.estimates(atlantic.city$max[23:103], atlantic.city.temps$values[23:103])
atlantic.city.90 <-  get.prior.estimates(atlantic.city$max[13:103], atlantic.city.temps$values[13:103])
atlantic.city.all <-  get.prior.estimates(atlantic.city$max, atlantic.city.temps$values)

portland.10 <- get.prior.estimates(portland$max[95:105], portland.temps$values[95:105])
portland.20 <- get.prior.estimates(portland$max[85:105], portland.temps$values[85:105])
portland.30 <- get.prior.estimates(portland$max[75:105], portland.temps$values[75:105])
portland.40 <- get.prior.estimates(portland$max[65:105], portland.temps$values[65:105])
portland.50 <- get.prior.estimates(portland$max[55:105], portland.temps$values[55:105])
portland.60 <- get.prior.estimates(portland$max[45:105], portland.temps$values[45:105])
portland.70 <- get.prior.estimates(portland$max[35:105], portland.temps$values[35:105])
portland.80 <- get.prior.estimates(portland$max[25:105], portland.temps$values[25:105])
portland.90 <- get.prior.estimates(portland$max[15:105], portland.temps$values[15:105])
portland.all <- get.prior.estimates(portland$max, portland.temps$values)
save.image('10.year.gaps.RData')
#all stationary case comparison------------------------------------------------- 
par(mfrow = c(1,3))
#new london
plot(30, new.london.30$stat[1], xlim = c(25, 110), ylim = c(1000,2500))
points(45, new.london.45$stat[1], col = 'red')
points(60, new.london.60$stat[1], col = 'blue')
points(76, new.london.all$stat[1], col = 'green')

#boston 
points(30,boston.30$stat[1], pch = 0)
points(45, boston.45$stat[1], col = 'red', pch = 0)
points(60, boston.60$stat[1], col = 'blue', pch = 0)
points(94, boston.all$stat[1], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$stat[1], pch = 2)
points(45, atlantic.city.45$stat[1], col = 'red', pch = 2)
points(60, atlantic.city.60$stat[1], col = 'blue', pch = 2)
points(103,atlantic.city.all$stat[1], col = 'green', pch = 2)

#portland 
points(30,portland.30$stat[1], pch = 3)
points(45, portland.45$stat[1], col = 'red', pch = 3)
points(60, portland.60$stat[1], col = 'blue', pch = 3)
points(105,portland.all$stat[1], col = 'green', pch = 3)

legend(20, 1800, legend= c('30 year', '45 year', '60 year', 'all year'), lty =c(1,1,1,1), col=c('black', 'red','blue', 'green'))
legend(60, 1800, legend= c('New London', 'Boston', 'Atlantic City', 'Portland'), 
       pch =c(1,0,2,3), 
       col=c('black', 'black','black', 'black'))



#sigma 
#new london
plot(30, new.london.30$stat[2], xlim = c(25, 110), ylim = c(50,200))
points(45, new.london.45$stat[2], col = 'red')
points(60, new.london.60$stat[2], col = 'blue')
points(76, new.london.all$stat[2], col = 'green')

#boston 
points(30,boston.30$stat[2], pch = 0)
points(45, boston.45$stat[2], col = 'red', pch = 0)
points(60, boston.60$stat[2], col = 'blue', pch = 0)
points(94, boston.all$stat[2], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$stat[2], pch = 2)
points(45, atlantic.city.45$stat[2], col = 'red', pch = 2)
points(60, atlantic.city.60$stat[2], col = 'blue', pch = 2)
points(103,atlantic.city.all$stat[2], col = 'green', pch = 2)

#portland 
points(30,portland.30$stat[2], pch = 3)
points(45, portland.45$stat[2], col = 'red', pch = 3)
points(60, portland.60$stat[2], col = 'blue', pch = 3)
points(105,portland.all$stat[2], col = 'green', pch = 3)


#xi
#new london
plot(30, new.london.30$stat[3], xlim = c(25, 110), ylim = c(-.5,.5))
points(45, new.london.45$stat[3], col = 'red')
points(60, new.london.60$stat[3], col = 'blue')
points(76, new.london.all$stat[3], col = 'green')

#boston 
points(30,boston.30$stat[3], pch = 0)
points(45, boston.45$stat[3], col = 'red', pch = 0)
points(60, boston.60$stat[3], col = 'blue', pch = 0)
points(94, boston.all$stat[3], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$stat[3], pch = 2)
points(45, atlantic.city.45$stat[3], col = 'red', pch = 2)
points(60, atlantic.city.60$stat[3], col = 'blue', pch = 2)
points(103,atlantic.city.all$stat[3], col = 'green', pch = 2)

#portland 
points(30,portland.30$stat[3], pch = 3)
points(45, portland.45$stat[3], col = 'red', pch = 3)
points(60, portland.60$stat[3], col = 'blue', pch = 3)
points(105,portland.all$stat[3], col = 'green', pch = 3)
#-----------------------------------------------------------------------------------

#All non stationary comparison-----------------------------------------------------
par(mfrow = c(2,3))
#mu0
#new london
plot(30, new.london.30$mu.sigma.xi[1], xlim = c(25, 110), ylim = c(1000,2500))
points(45, new.london.45$mu.sigma.xi[1], col = 'red')
points(60, new.london.60$mu.sigma.xi[1], col = 'blue')
points(76, new.london.all$mu.sigma.xi[1], col = 'green')

#boston 
points(30,boston.30$mu.sigma.xi[1], pch = 0)
points(45, boston.45$mu.sigma.xi[1], col = 'red', pch = 0)
points(60, boston.60$mu.sigma.xi[1], col = 'blue', pch = 0)
points(94, boston.all$mu.sigma.xi[1], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$mu.sigma.xi[1], pch = 2)
points(45, atlantic.city.45$mu.sigma.xi[1], col = 'red', pch = 2)
points(60, atlantic.city.60$mu.sigma.xi[1], col = 'blue', pch = 2)
points(103,atlantic.city.all$mu.sigma.xi[1], col = 'green', pch = 2)

#portland 
points(30,portland.30$mu.sigma.xi[1], pch = 3)
points(45, portland.45$mu.sigma.xi[1], col = 'red', pch = 3)
points(60, portland.60$mu.sigma.xi[1], col = 'blue', pch = 3)
points(105,portland.all$mu.sigma.xi[1], col = 'green', pch = 3)

legend(30, 2000, legend= c('30 year', '45 year', '60 year', 'all year'), 
       lty =c(1,1,1,1), col=c('black', 'red','blue', 'green'), cex = .5)
legend(60,2000, legend= c('New London', 'Boston', 'Atlantic City', 'Portland'), 
       pch =c(1,0,2,3), 
       col=c('black', 'black','black', 'black'), cex = .5)


#mu1
#new london
plot(30, new.london.30$mu.sigma.xi[2], xlim = c(25, 110), ylim = c(-100,100))
points(45, new.london.45$mu.sigma.xi[2], col = 'red')
points(60, new.london.60$mu.sigma.xi[2], col = 'blue')
points(76, new.london.all$mu.sigma.xi[2], col = 'green')

#boston 
points(30,boston.30$mu.sigma.xi[2], pch = 0)
points(45, boston.45$mu.sigma.xi[2], col = 'red', pch = 0)
points(60, boston.60$mu.sigma.xi[2], col = 'blue', pch = 0)
points(94, boston.all$mu.sigma.xi[2], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$mu.sigma.xi[2], pch = 2)
points(45, atlantic.city.45$mu.sigma.xi[2], col = 'red', pch = 2)
points(60, atlantic.city.60$mu.sigma.xi[2], col = 'blue', pch = 2)
points(103,atlantic.city.all$mu.sigma.xi[2], col = 'green', pch = 2)

#portland 
points(30,portland.30$mu.sigma.xi[3], pch = 3)
points(45, portland.45$mu.sigma.xi[3], col = 'red', pch = 3)
points(60, portland.60$mu.sigma.xi[3], col = 'blue', pch = 3)
points(105,portland.all$mu.sigma.xi[3], col = 'green', pch = 3) 

#sigma0 
#new london
plot(30, new.london.30$mu.sigma.xi[3], xlim = c(25, 110), ylim = c(0,10))
points(45, new.london.45$mu.sigma.xi[3], col = 'red')
points(60, new.london.60$mu.sigma.xi[3], col = 'blue')
points(76, new.london.all$mu.sigma.xi[3], col = 'green')

#boston 
points(30,boston.30$mu.sigma.xi[3], pch = 0)
points(45, boston.45$mu.sigma.xi[3], col = 'red', pch = 0)
points(60, boston.60$mu.sigma.xi[3], col = 'blue', pch = 0)
points(94, boston.all$mu.sigma.xi[3], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$mu.sigma.xi[3], pch = 2)
points(45, atlantic.city.45$mu.sigma.xi[3], col = 'red', pch = 2)
points(60, atlantic.city.60$mu.sigma.xi[3], col = 'blue', pch = 2)
points(103,atlantic.city.all$mu.sigma.xi[3], col = 'green', pch = 2)

#portland 
points(30,portland.30$mu.sigma.xi[3], pch = 3)
points(45, portland.45$mu.sigma.xi[3], col = 'red', pch = 3)
points(60, portland.60$mu.sigma.xi[3], col = 'blue', pch = 3)
points(105,portland.all$mu.sigma.xi[3], col = 'green', pch = 3) 

#sigma1
#new london
plot(30, new.london.30$mu.sigma.xi[4], xlim = c(25, 110), ylim = c(-5,5))
points(45, new.london.45$mu.sigma.xi[4], col = 'red')
points(60, new.london.60$mu.sigma.xi[4], col = 'blue')
points(76, new.london.all$mu.sigma.xi[4], col = 'green')

#boston 
points(30,boston.30$mu.sigma.xi[4], pch = 0)
points(45, boston.45$mu.sigma.xi[4], col = 'red', pch = 0)
points(60, boston.60$mu.sigma.xi[4], col = 'blue', pch = 0)
points(94, boston.all$mu.sigma.xi[4], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$mu.sigma.xi[4], pch = 2)
points(45, atlantic.city.45$mu.sigma.xi[4], col = 'red', pch = 2)
points(60, atlantic.city.60$mu.sigma.xi[4], col = 'blue', pch = 2)
points(103,atlantic.city.all$mu.sigma.xi[4], col = 'green', pch = 2)

#portland 
points(30,portland.30$mu.sigma.xi[4], pch = 3)
points(45, portland.45$mu.sigma.xi[4], col = 'red', pch = 3)
points(60, portland.60$mu.sigma.xi[4], col = 'blue', pch = 3)
points(105,portland.all$mu.sigma.xi[4], col = 'green', pch = 3) 

#xi0
#new london
plot(30, new.london.30$mu.sigma.xi[5], xlim = c(25, 110), ylim = c(-1,2))
points(45, new.london.45$mu.sigma.xi[5], col = 'red')
points(60, new.london.60$mu.sigma.xi[5], col = 'blue')
points(76, new.london.all$mu.sigma.xi[5], col = 'green')

#boston 
points(30,boston.30$mu.sigma.xi[5], pch = 0)
points(45, boston.45$mu.sigma.xi[5], col = 'red', pch = 0)
points(60, boston.60$mu.sigma.xi[5], col = 'blue', pch = 0)
points(94, boston.all$mu.sigma.xi[5], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$mu.sigma.xi[5], pch = 2)
points(45, atlantic.city.45$mu.sigma.xi[5], col = 'red', pch = 2)
points(60, atlantic.city.60$mu.sigma.xi[5], col = 'blue', pch = 2)
points(103,atlantic.city.all$mu.sigma.xi[5], col = 'green', pch = 2)

#portland 
points(30,portland.30$mu.sigma.xi[5], pch = 3)
points(45, portland.45$mu.sigma.xi[5], col = 'red', pch = 3)
points(60, portland.60$mu.sigma.xi[5], col = 'blue', pch = 3)
points(105,portland.all$mu.sigma.xi[5], col = 'green', pch = 3) 

#xi1
#new london
plot(30, new.london.30$mu.sigma.xi[6], xlim = c(25, 110), ylim = c(-2,2))
points(45, new.london.45$mu.sigma.xi[6], col = 'red')
points(60, new.london.60$mu.sigma.xi[6], col = 'blue')
points(76, new.london.all$mu.sigma.xi[6], col = 'green')

#boston 
points(30,boston.30$mu.sigma.xi[6], pch = 0)
points(45, boston.45$mu.sigma.xi[6], col = 'red', pch = 0)
points(60, boston.60$mu.sigma.xi[6], col = 'blue', pch = 0)
points(94, boston.all$mu.sigma.xi[6], col = 'green', pch = 0)

#atlantic city 
points(30,atlantic.city.30$mu.sigma.xi[6], pch = 2)
points(45, atlantic.city.45$mu.sigma.xi[6], col = 'red', pch = 2)
points(60, atlantic.city.60$mu.sigma.xi[6], col = 'blue', pch = 2)
points(103,atlantic.city.all$mu.sigma.xi[6], col = 'green', pch = 2)

#portland 
points(30,portland.30$mu.sigma.xi[6], pch = 3)
points(45, portland.45$mu.sigma.xi[6], col = 'red', pch = 3)
points(60, portland.60$mu.sigma.xi[6], col = 'blue', pch = 3)
points(105,portland.all$mu.sigma.xi[6], col = 'green', pch = 3) 
#-----------------------------------------------------------------------------


#with 10 year plot segments---------------------------------------------------
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
par(mar = c(2,4,0,1))
nl.xi0 <- c(new.london.10$mu.sigma.xi[5],
               new.london.20$mu.sigma.xi[5],
               new.london.30$mu.sigma.xi[5],
               new.london.40$mu.sigma.xi[5],
               new.london.50$mu.sigma.xi[5],
               new.london.60$mu.sigma.xi[5],
               new.london.all$mu.sigma.xi[5])
plot(nl.years, nl.xi0,xlim = c(0, 110), ylim = c(-1.5,1.5), type = 'l', ylab = 'xi0')
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
plot(nl.years, nl.xi1, xlim = c(0, 110), ylim = c(-1.5,1.5), type = 'l', ylab = 'xi1')
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