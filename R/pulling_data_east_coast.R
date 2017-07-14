#Script to pull all of the other tide data from places on the east coast

dat.dir <- './data/'
filetype <- 'csv'
septype <- ','

#naming code 
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.projout <- paste('../output_model/BRICK_project-lsl-surge_NewLondon_',today,'.nc',sep='')
filename.lslout  <- paste('../output_model/BRICK_project-lsl_NewLondon_',today,'.csv', sep="")

##==============================================================================
## Read tide gauge data
dat.dir <- './data2/'
filetype <- 'csv'
septype <- ','

#setwd('~/codes/Klufas_NewLondon/')

#files.tg <- list.files(path=dat.dir,pattern=filetype)

setwd('~/codes/Klufas_NewLondon/')

files.tg <- list.files(path=dat.dir,pattern=filetype)

AtlanticCity <- read.csv(paste(dat.dir,files.tg[1],sep=''), header=TRUE, sep=septype)
colnames(AtlanticCity) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Boston <-read.csv(paste(dat.dir,files.tg[2],sep=''), header=TRUE, sep=septype)
colnames(Boston) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Charleston <- read.csv(paste(dat.dir,files.tg[3],sep=''), header=TRUE, sep=septype)
colnames(Charleston) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
ChesapeakeBay <-read.csv(paste(dat.dir,files.tg[4],sep=''), header=TRUE, sep=septype)
colnames(ChesapeakeBay) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
DuckPier <- read.csv(paste(dat.dir,files.tg[5],sep=''), header=TRUE, sep=septype)
colnames(DuckPier) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Lewes <-read.csv(paste(dat.dir,files.tg[6],sep=''), header=TRUE, sep=septype)
colnames(Lewes) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Montauk <- read.csv(paste(dat.dir,files.tg[7],sep=''), header=TRUE, sep=septype)
colnames(Montauk) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Newport <- read.csv(paste(dat.dir,files.tg[8],sep=''), header=TRUE, sep=septype)
colnames(Newport) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
NewYorkCity <-read.csv(paste(dat.dir,files.tg[9],sep=''), header=TRUE, sep=septype)
colnames(NewYorkCity) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Portland <- read.csv(paste(dat.dir,files.tg[10],sep=''), header=TRUE, sep=septype)
colnames(Portland) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')

#city.list <- c(AtlanticCity, Boston,Charleston,ChesapeakeBay,DuckPier,Lewes,Montauk,Newport,NewYorkCity , Portland) 

read.data <- function (city){
  city.data <- vector('list', 3)
  names(city.data) <- c('max', 'mean','years')
  years         <- city$Year
  years.unique  <- unique(years)
  n.years       <- length(years.unique) 
  lsl.mean      <- rep(0,length(n.years))
  lsl.max       <- rep(0,length(n.years)) #max(data$lsl.norm[ind.thisyear])
  city$lsl.norm <- rep(NA,length(years))


##----------Finding Annual Block Maxima----------------------
#started at 2 so that the data the data stuff would start at 1939 not 1938 
  for (tt in 1:n.years) {
    ind.thisyear <- which(years==years.unique[tt])
    lsl.mean[tt] <- mean(city$Sea_Level[ind.thisyear])
    city$lsl.norm[ind.thisyear] <- city$Sea_Level[ind.thisyear] - lsl.mean[tt]
    lsl.max[tt] <- max(city$lsl.norm[ind.thisyear])
  }
  city.data$years <- years.unique 
  city.data$max <- lsl.max
  city.data$mean <- lsl.mean
  
  return(city.data)
}

atlantic.city <- read.data(AtlanticCity)
boston <- read.data(Boston)
chesapeake.bay <- read.data(ChesapeakeBay)
lewes<- read.data(Lewes)
new.york.city <- read.data(NewYorkCity)
portland <- read.data(Portland)
charleston <- read.data(Charleston)
duck.pier <- read.data(DuckPier)
montauk <- read.data(Montauk)
newport <- read.data(Newport)



library(ismev)
library(zoo)

gev.mle.atlantic.city <- fevd(coredata(atlantic.city$max), type='GEV') # extRemes
gev.mle2.atlantic.city <- gev.fit(coredata(atlantic.city$max), show = FALSE) 
print(gev.mle.atlantic.city$results$par)

gev.mle.boston <- fevd(coredata(boston$max), type='GEV') # extRemes
gev.mle2.boston <- gev.fit(coredata(boston$max), show = FALSE) 
print(gev.mle.boston$results$par)

gev.mle.chesapeake.bay <- fevd(coredata(chesapeake.bay$max), type='GEV') # extRemes
gev.mle2.chesapeake.bay <- gev.fit(coredata(chesapeake.bay$max), show = FALSE) 
print(gev.mle.chesapeake.bay$results$par)

gev.mle.lewes <- fevd(coredata(lewes$max), type='GEV') # extRemes
gev.mle2.lewes <- gev.fit(coredata(lewes$max), show = FALSE) 
print(gev.mle.lewes$results$par)

gev.mle.new.york.city <- fevd(coredata(new.york.city$max), type='GEV') # extRemes
gev.mle2.new.york.city <- gev.fit(coredata(new.york.city$max), show = FALSE) 
print(gev.mle.new.york.city$results$par)

gev.mle.portland <- fevd(coredata(portland$max), type='GEV') # extRemes
gev.mle2.portland <- gev.fit(coredata(portland$max), show = FALSE) 
print(gev.mle.portland$results$par)

gev.mle.charleston <- fevd(coredata(charleston$max), type='GEV') # extRemes
gev.mle2.charleston <- gev.fit(coredata(charleston$max), show = FALSE) 
print(gev.mle.charleston$results$par)

gev.mle.duck.pier <- fevd(coredata(duck.pier$max), type='GEV') # extRemes
gev.mle2.duck.pier <- gev.fit(coredata(duck.pier$max), show = FALSE) 
print(gev.mle.duck.pier$results$par)

gev.mle.montauk <- fevd(coredata(montauk$max), type='GEV') # extRemes
gev.mle2.montauk <- gev.fit(coredata(montauk$max), show = FALSE) 
print(gev.mle.montauk$results$par)

gev.mle.newport  <- fevd(coredata(newport$max), type='GEV') # extRemes
gev.mle2.montauk <- gev.fit(coredata(newport$max), show = FALSE) 
print(gev.mle.newport$results$par)

setwd('~/codes/Klufas_NewLondon/R/')
library('extRemes')
source('optimization_sf.R')

city.list.data <- list(atlantic.city, boston, chesapeake.bay, lewes, new.york.city, portland, charleston, duck.pier, montauk, newport )

all.p.names <- list(c('mu', 'sigma', 'xi'), 
                 c('mu0', 'mu1', 'sigma', 'xi'),
                 c('mu', 'sigma0', 'sigma1', 'xi'),
                 c('mu', 'sigma', 'xi0', 'xi1'),
                 c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi'),
                 c('mu', 'sigma0', 'sigma1' , 'xi0', 'xi1'),
                 c('mu0', 'mu1', 'sigma', 'xi0', 'xi1'),
                 c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi0', 'xi1'))

other.p.names <- c('stat', 'mu','sigma', 'xi', 'mu.sigma', 'sigma.xi', 'mu.xi', 'mu.sigma.xi')
city.names <- c('atlantic.city', 'boston', 'chesapeake.bay', 'lewes', 'new.york.city', 'portland', 'duck.pier', 'montauk', 'newport','charleston')
optim.gev.fit <- vector('list' , 10)
names(optim.gev.fit) <- city.names
lower.bound <- list(c(0,0,-5),
                 c(0,0, 0, -5),
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

get.prior.estimates <- function(city, city.temps){
  #for (i in 1:2){
  
    sublist <- vector('list', 8)
    names(sublist) <- other.p.names
     for (j in 1: length(all.p.names)){
      print(j)
        optim.like <- DEoptim(neg.log.like.calc,
                                  lower=lower.bound[[j]], 
                                  upper=upper.bound[[j]], 
                                  data = city$max, 
                                  temps = city.temps$values, 
                                  parnames = all.p.names[[j]])
        sublist[[other.p.names[j]]] <- optim.like$optim$bestmem
    }
    #optim.gev.fit$city.names[i] <- sublist
  return(sublist)
  }

setwd('~/codes/Klufas_NewLondon/R/')
source('read_temp_data.R')
#get the temps from all of the cities
boston.temps <- read.temp.data(1921, 2014)
atlantic.city.temps <- read.temp.data(1911,2014)
charleston.temps <- read.temp.data(1921, 2017)
montauk.temps <- read.temp.data(1959, 2014)
chesapeake.bay.temps <- read.temp.data(1975,2014)
lewes.temps <- read.temp.data(1957,2014)
new.york.city.temps <- read.temp.data(1920, 2014)
portland.temps <- read.temp.data(1910, 2014)
duck.pier.temps <- read.temp.data(1978, 2017)
newport.temps <- read.temp.data(1910, 2014)

boston.priors <- get.prior.estimates(boston, boston.temps)
atlantic.city.priors <- get.prior.estimates(atlantic.city, atlantic.city.temps)
charleston.priors <- get.prior.estimates(charleston, charleston.temps)
montauk.priors <- get.prior.estimates(montauk, montauk.temps)
chesapeake.bay.priors<- get.prior.estimates(chesapeake.bay, chesapeake.bay.temps)
lewes.priors<- get.prior.estimates(lewes,lewes.temps)
new.york.city.priors<- get.prior.estimates(new.york.city,new.york.city.temps)
portland.priors <- get.prior.estimates(portland,portland.temps)
duck.pier.priors<- get.prior.estimates(duck.pier,duck.pier.temps)
newport.priors <- get.prior.estimates(newport,newport.temps)

mega.priors.list <- list(boston.priors, 
                         atlantic.city.priors, 
                         charleston.priors,
                         montauk.priors, 
                         chesapeake.bay.priors, 
                         lewes.priors, 
                         new.york.city.priors, 
                         portland.priors, 
                         duck.pier.priors, 
                         newport.priors)

mu.all.stat <- c(boston.priors$stat[1], 
                 atlantic.city.priors$stat[1], 
                 charleston.priors$stat[1],
                 montauk.priors$stat[1], 
                 chesapeake.bay.priors$stat[1], 
                 lewes.priors$stat[1], 
                 new.york.city.priors$stat[1], 
                 portland.priors$stat[1], 
                 duck.pier.priors$stat[1], 
                 newport.priors$stat[1])

sigma.all.stat <- c(boston.priors$stat[2], 
                    atlantic.city.priors$stat[2], 
                    charleston.priors$stat[2],
                    montauk.priors$stat[2], 
                    chesapeake.bay.priors$stat[2], 
                    lewes.priors$stat[2], 
                    new.york.city.priors$stat[2], 
                    portland.priors$stat[2], 
                    duck.pier.priors$stat[2], 
                    newport.priors$stat[2])

xi.all.stat <-c(boston.priors$stat[3], 
                atlantic.city.priors$stat[3], 
                charleston.priors$stat[3],
                montauk.priors$stat[3], 
                chesapeake.bay.priors$stat[3], 
                lewes.priors$stat[3], 
                new.york.city.priors$stat[3], 
                portland.priors$stat[3], 
                duck.pier.priors$stat[3], 
                newport.priors$stat[3]) 

par(mfrow = c(1,3))
hist(mu.all.stat, freq = FALSE)
hist(sigma.all.stat, freq = FALSE)
hist(xi.all.stat, freq = FALSE)


mu0.mu.nonstat <- c(boston.priors$mu[1], 
                 atlantic.city.priors$mu[1], 
                 charleston.priors$mu[1],
                 montauk.priors$mu[1], 
                 chesapeake.bay.priors$mu[1], 
                 lewes.priors$mu[1], 
                 new.york.city.priors$mu[1], 
                 portland.priors$mu[1], 
                 duck.pier.priors$mu[1], 
                 newport.priors$mu[1])

mu1.mu.nonstat <- c(boston.priors$mu[2], 
                    atlantic.city.priors$mu[2], 
                    charleston.priors$mu[2],
                    montauk.priors$mu[2], 
                    chesapeake.bay.priors$mu[2], 
                    lewes.priors$mu[2], 
                    new.york.city.priors$mu[2], 
                    portland.priors$mu[2], 
                    duck.pier.priors$mu[2], 
                    newport.priors$mu[2])

sigma.mu.nonstat <-c(boston.priors$mu[3], 
                     atlantic.city.priors$mu[3], 
                     charleston.priors$mu[3],
                     montauk.priors$mu[3], 
                     chesapeake.bay.priors$mu[3], 
                     lewes.priors$mu[3], 
                     new.york.city.priors$mu[3], 
                     portland.priors$mu[3], 
                     duck.pier.priors$mu[3], 
                     newport.priors$mu[3])

xi.mu.nonstat <- c(boston.priors$mu[4], 
                   atlantic.city.priors$mu[4], 
                   charleston.priors$mu[4],
                   montauk.priors$mu[4], 
                   chesapeake.bay.priors$mu[4], 
                   lewes.priors$mu[4], 
                   new.york.city.priors$mu[4], 
                   portland.priors$mu[4], 
                   duck.pier.priors$mu[4], 
                   newport.priors$mu[4])

par(mfrow = c(2,2))
hist(mu0.mu.nonstat, freq= FALSE) 
hist(mu1.mu.nonstat, freq= FALSE)
hist(sigma.mu.nonstat, freq= FALSE)
hist(xi.mu.nonstat, freq= FALSE)


mu0.all.nonstat <- c(boston.priors$mu.sigma.xi[1], 
                    atlantic.city.priors$mu.sigma.xi[1], 
                    charleston.priors$mu.sigma.xi[1],
                    montauk.priors$mu.sigma.xi[1], 
                    chesapeake.bay.priors$mu.sigma.xi[1], 
                    lewes.priors$mu.sigma.xi[1], 
                    new.york.city.priors$mu.sigma.xi[1], 
                    portland.priors$mu.sigma.xi[1], 
                    duck.pier.priors$mu.sigma.xi[1], 
                    newport.priors$mu.sigma.xi[1])

mu1.all.nonstat <- c(boston.priors$mu.sigma.xi[2], 
                    atlantic.city.priors$mu.sigma.xi[2], 
                    charleston.priors$mu.sigma.xi[2],
                    montauk.priors$mu.sigma.xi[2], 
                    chesapeake.bay.priors$mu.sigma.xi[2], 
                    lewes.priors$mu.sigma.xi[2], 
                    new.york.city.priors$mu.sigma.xi[2], 
                    portland.priors$mu.sigma.xi[2], 
                    duck.pier.priors$mu.sigma.xi[2], 
                    newport.priors$mu.sigma.xi[2])

sigma0.all.nonstat <-c(boston.priors$mu.sigma.xi[3], 
                     atlantic.city.priors$mu.sigma.xi[3], 
                     charleston.priors$mu.sigma.xi[3],
                     montauk.priors$mu.sigma.xi[3], 
                     chesapeake.bay.priors$mu.sigma.xi[3], 
                     lewes.priors$mu.sigma.xi[3], 
                     new.york.city.priors$mu.sigma.xi[3], 
                     portland.priors$mu.sigma.xi[3], 
                     duck.pier.priors$mu.sigma.xi[3], 
                     newport.priors$mu.sigma.xi[3])

sigma1.all.nonstat <- c(boston.priors$mu.sigma.xi[4], 
                   atlantic.city.priors$mu.sigma.xi[4], 
                   charleston.priors$mu.sigma.xi[4],
                   montauk.priors$mu.sigma.xi[4], 
                   chesapeake.bay.priors$mu.sigma.xi[4], 
                   lewes.priors$mu.sigma.xi[4], 
                   new.york.city.priors$mu.sigma.xi[4], 
                   portland.priors$mu.sigma.xi[4], 
                   duck.pier.priors$mu.sigma.xi[4], 
                   newport.priors$mu.sigma.xi[4])

xi0.all.nonstat <- c(boston.priors$mu.sigma.xi[5], 
                     atlantic.city.priors$mu.sigma.xi[5], 
                     charleston.priors$mu.sigma.xi[5],
                     montauk.priors$mu.sigma.xi[5], 
                     chesapeake.bay.priors$mu.sigma.xi[5], 
                     lewes.priors$mu.sigma.xi[5], 
                     new.york.city.priors$mu.sigma.xi[5], 
                     portland.priors$mu.sigma.xi[5], 
                     duck.pier.priors$mu.sigma.xi[5], 
                     newport.priors$mu.sigma.xi[5])

xi1.all.nonstat <- c(boston.priors$mu.sigma.xi[6], 
                     atlantic.city.priors$mu.sigma.xi[6], 
                     charleston.priors$mu.sigma.xi[6],
                     montauk.priors$mu.sigma.xi[6], 
                     chesapeake.bay.priors$mu.sigma.xi[6], 
                     lewes.priors$mu.sigma.xi[6], 
                     new.york.city.priors$mu.sigma.xi[6], 
                     portland.priors$mu.sigma.xi[6], 
                     duck.pier.priors$mu.sigma.xi[6], 
                     newport.priors$mu.sigma.xi[6])
par(mfrow = c(2,3))
hist(mu0.all.nonstat, freq= FALSE) 
hist(mu1.all.nonstat, freq= FALSE)
plot(dnorm(mu1.all.nonstat), col = 'red', type = 'l')
hist(sigma0.all.nonstat, freq= FALSE)
hist(sigma1.all.nonstat, freq= FALSE)
hist(xi0.all.nonstat, freq= FALSE)
hist(xi1.all.nonstat, freq= FALSE)
#hist.prep <- vector('list', 5)
#names(hist.prep) <- c('all.stat.mu','all.stat.sigma','all.stat.xi', 'mu.nonstat.mu0', 'mu.nonstat.mu1')
#for( i in 1:length((mega.priors.list))){
  
#}
install.packages('fitdistrplus')
library(fitdistrplus)
dist <- fitdistr(mu1.all.nonstat, 'normal')
plot(dist)

save.image('east.coast.RData')
