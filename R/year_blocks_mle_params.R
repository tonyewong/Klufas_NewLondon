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

load('east.coast.RData')
boston.30 <- get.prior.estimates(boston$max[64:94], boston.temps$values[64:94])
boston.45 <- get.prior.estimates(boston$max[49:94], boston.temps$values[49:94])
boston.60 <- get.prior.estimates(boston$max[34:94], boston.temps$values[34:94])
boston.all <- get.prior.estimates(boston$max, boston.temps$values)

atlantic.city.30 <-  get.prior.estimates(atlantic.city$max[73:103], atlantic.city.temps$values[73:103])
atlantic.city.45 <-  get.prior.estimates(atlantic.city$max[58:103], atlantic.city.temps$values[58:103])
atlantic.city.60 <-  get.prior.estimates(atlantic.city$max[43:103], atlantic.city.temps$values[43:103])
atlantic.city.all <-  get.prior.estimates(atlantic.city$max, atlantic.city.temps$values)

portland.30 <- get.prior.estimates(portland$max[75:105], portland.temps$values[75:105])
portland.45 <- get.prior.estimates(portland$max[60:105], portland.temps$values[60:105])
portland.60 <- get.prior.estimates(portland$max[45:105], portland.temps$values[45:105])
portland.all <- get.prior.estimates(portland$max, portland.temps$values)

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
#-----------------------------------------------------------------------------------
