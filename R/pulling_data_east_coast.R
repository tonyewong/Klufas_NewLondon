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


setwd('~/codes/Klufas_NewLondon/')
dat.dir <- './data/'
files.tg2 <- list.files(path=dat.dir,pattern=filetype)
NewLondon <- read.csv(paste(dat.dir, files.tg2, sep=''), header = TRUE, sep = septype)
colnames(NewLondon) <- c('Year','Month','Day','Hour','Sea_Level')
setwd('~/codes/Klufas_NewLondon/')
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

  
  #EXTRA ADDED STUFF ___________________________________
  #get all of the hours out of data 
  hours <- city$Hour
  #get the unique hours  expect there to be 24 (and there are!))
  hours.unique <- unique(hours)
  n.hours <- length(hours)
  
  #finding the difference in hours
  hour.diff <- diff(hours)
  #cheking to see if any places dont have numbers
  for (i in 1:n.hours){
    if (is.nan(hour.diff[i])){
      print("help")
    }
  }
  #conclusion - no nans 
  
  #for loop through all hours 
  
  #the hours that are weird
  day.counter <- 0
  month.counter <- 0
  #subtracting one to make it the same length as the differnece in hours 
  num <- 1
  hours.strange <- rep(0, length(n.hours)) 
   
  for (i in 1:(n.hours-1)){
    #print(hour.diff[i])
    #checking if difference between hours is not 1 or -23, will print 'index of data pt that has problem', else does nothing 
    if (hour.diff[i] != 1 & hour.diff[i] != -23){
      print(i+1)
      #the list of all of the weird hour gaps in the data 
      hours.strange[num] <- i+1
      num <- num + 1
    }
  }
  
  #now to check how big the gap is 
  days <- city$Day
  day.unique <- unique(days)
  n.days <- length(days)
  
  #also need to check to see if there is a month gap because that could also pose the same type of problem 
  
  months <- city$Month
  n.month <- length(months)
  
 
  #look at each index and the one after and check to see if days are within a day of each other 
  # if so - check hour difference 
  #else - check if in different months - cases -31, -30, -28, -29 (for feb)
  #go through the list of weird numbers we have 
  for (j in 1:length(hours.strange)){
    #need to check within the same day of each other - so that would be seeing if the difference between the days is 0 
    #index we will be using comes from the list that we just made (the list hours.strange)
    ind <- hours.strange[[j]]
    #print(j)
    # print(ind)
    #if the difference between the two data pts is in the same day #lets ignore it 
    # print(days[ind])
    if ((days[ind]-days[ind-1] ) == 0 & (months[ind] - months[ind-1]) == 0){
      print("in the same day")
    }else if(years[ind]- years[ind-1] > 1){
      print("big gap here**************************************************")
      print(ind)
    }else if ((months[ind] - months[ind-1]) >= 1 | (months[ind] - months[ind-1]) <= -1) { #if month gap greater than 1, make note of it
      if (abs(years[ind]- years[ind-1] > 1)){
        print("month gap greater than or equal to 1")
        print(abs(months[ind] - months[ind-1]))
        print(ind)
        month.counter <- month.counter + abs(months[ind] - months[ind-1])
      }}
    #when the days are not the same , need to check how far apart 
    else {
      print("no month gap, only day gap")
      print(days[ind]-days[ind-1])
      day.counter <- day.counter + days[ind]-days[ind-1]
      print(day.counter)
    }
  }
  
  #EXTRA ADDED STUFF ___________________________________
##----------Finding Annual Block Maxima----------------------
#started at 2 so that the data the data stuff would start at 1939 not 1938 
  for (tt in 1:n.years) {
    ind.thisyear <- which(years==years.unique[tt])
    if (length(ind.thisyear) >= 7884){
      lsl.mean[tt] <- mean(city$Sea_Level[ind.thisyear])
      city$lsl.norm[ind.thisyear] <- city$Sea_Level[ind.thisyear] - lsl.mean[tt]
      lsl.max[tt] <- max(city$lsl.norm[ind.thisyear])
    } #90% of number of hours in a year 
    else{
      lsl.mean[tt] <- 0 
      city$lsl.norm[ind.thisyear] <- 0 
      lsl.max[tt] <- 0
    }
  }
  
  #for(q in 1:length(years.unique)){
  ind <- which(lsl.max != 0)
  years.shift <- years.unique[ind]
  sl <- lsl.max[ind]
  sl.mean <- lsl.mean[ind]
    
  #}
  
  
  city.data$years <- years.shift 
  city.data$max <- sl
  city.data$mean <- sl.mean
  
  
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
new.london <- read.data(NewLondon)
save.image('all.loc.less.year.RData')

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
tmp <- c('stat2', 'mu2','sigma2', 'xi2', 'mu.sigma2', 'sigma.xi2', 'mu.xi2', 'mu.sigma.xi2')
optim.gev.fit <- vector('list' , 10)
names(optim.gev.fit) <- city.names
lower.bound <- list(c(0,0,-5),
                 c(0,-100, 0, -5),
                 c(0, 0, 0, -200),
                 c(0, 0, -5, -5), 
                 c(0,-100, -100, -100, -5),
                 c(0, -100, -100, -1, -1),
                 c(0, -100, 0, -1, -10),
                 c(0,-100, 0, -100, -1, -1)
                 )
upper.bound <- list(c(4000,3000,5),
                 c(4000,1000, 1000, 5),
                 c(4000, 1000, 1000, 1000),
                 c(4000, 1000 , 100, 100),
                 c(4000,100, 1000, 10 , 5),
                 c(4000, 1000, 10 , 1, 1),
                 c(4000, 100, 1000, 10, 10),
                 c(4000,100, 1000, 10 , 1, 1)
                 )

get.prior.estimates <- function(city, city.temps){
  #for (i in 1:2){
  
    sublist <- vector('list', 8)
    #sublist2 <- vector('list', 8)
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
        #sublist2[[tmp[j]]] <- optim.like$member$bestmemit
    }
    #optim.gev.fit$city.names[i] <- sublist
    #ouble.list <- vector('list', 2)
    #names(double.list) <- c('sublist', 'sublist2')
   # double.list$sublist <- sublist
    #double.list$sublist2 <- sublist2
  return(sublist)
  }

setwd('~/codes/Klufas_NewLondon/R/')
source('read_temp_data.R')
#get the temps from all of the cities
boston.temps <- read.temp.data(boston$years)
atlantic.city.temps <- read.temp.data(atlantic.city$years)
charleston.temps <- read.temp.data(charleston$years)
montauk.temps <- read.temp.data(montauk$years)
chesapeake.bay.temps <- read.temp.data(chesapeake.bay$years)
lewes.temps <- read.temp.data(lewes$years)
new.york.city.temps <- read.temp.data(new.york.city$years)
portland.temps <- read.temp.data(portland$years)
duck.pier.temps <- read.temp.data(duck.pier$years)
newport.temps <- read.temp.data(newport$years)
new.london.temps <- read.temp.data(new.london$years)
save.image('new.temps.RData')

boston.priors.2 <- get.prior.estimates(boston, boston.temps)
portland.priors.2 <- get.prior.estimates(portland,portland.temps)
save.image('priors.run2.RData')

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
#install.packages('MASS')
library(MASS)
#fitdistr(mu.all.stat, 'gamma')
par(mfrow = c(1,3), oma = c(4,2,4,0))

hist(mu.all.stat, freq = FALSE, main = '', xlab = expression(mu), ylab = '', cex.lab = 2)
mu.gamma <- fitdistr(mu.all.stat, 'gamma')
mu.seq <- seq(from = min (mu.all.stat) - 200, to = max(mu.all.stat) + 200, by = 1)
lines(mu.seq, dgamma(mu.seq, shape = mu.gamma$estimate[1], rate = mu.gamma$estimate[2]), col = 'red')

hist(sigma.all.stat, freq = FALSE, main='', xlab = expression(sigma), ylab = '', cex.lab = 2)
sigma.gamma <- fitdistr(sigma.all.stat, 'gamma')
seq.sigma <- seq(from = min(sigma.all.stat)-20, to = max(sigma.all.stat)+20, by =1)
lines( seq.sigma, dgamma(seq.sigma, shape = sigma.gamma$estimate[1], rate = sigma.gamma$estimate[2] ), col = 'red') 

hist(xi.all.stat, freq = FALSE, main ='', xlab = expression(xi), ylab = '' ,cex.lab = 2)
sd <- sd(xi.all.stat)
mean <- mean(xi.all.stat)
seq.xi <- seq(from =min(xi.all.stat) - .1, to =max(xi.all.stat) + .1, by = .01)
lines(seq.xi, dnorm(seq.xi, mean = mean, sd = sd), col = 'red')

par(oma = c(3,2,3,0))
mtext(text = 'Density', side = 2, outer = TRUE)
mtext(text = 'Parameter Values', side = 1, outer = TRUE)
title('Distribution of Stationary Parameters of Different Tide Stations', outer = TRUE)

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
par(oma = c(3,3,3,3))
hist(mu0.mu.nonstat, freq= FALSE, main = '', xlab = expression(mu[0]), ylab = '', cex.lab = 1.5) 
mu0.gamma <- fitdistr(mu0.mu.nonstat, 'gamma')
mu0.seq <- seq(from = min (mu0.mu.nonstat) - 500, to = max(mu0.mu.nonstat) + 300, by = 1)
lines(mu0.seq, dgamma(mu0.seq, shape = mu0.gamma$estimate[1], rate = mu0.gamma$estimate[2]), col = 'red')

hist(mu1.mu.nonstat, freq= FALSE, main = '', xlab = expression(mu[1]), ylab = '', cex.lab = 1.5)
sd <- sd(mu1.mu.nonstat)
mean <- mean(mu1.mu.nonstat)
mu1.seq <- seq(from = min (mu1.mu.nonstat) - 100, to = max(mu1.mu.nonstat) + 300, by = 1)
lines(mu1.seq, dnorm(mu1.seq, mean = mean, sd = sd), col = 'red')

hist(sigma.mu.nonstat, freq= FALSE, main = '', xlab = expression(sigma), ylab = '', cex.lab = 1.5)
sigma.gamma <- fitdistr(sigma.mu.nonstat, 'gamma')
sigma.seq <- seq(from = min (sigma.mu.nonstat) - 20, to = max(sigma.mu.nonstat) + 300, by = 1)
lines(sigma.seq, dgamma(sigma.seq, shape = sigma.gamma$estimate[1], rate = sigma.gamma$estimate[2]), col = 'red')

hist(xi.mu.nonstat, freq= FALSE, main = '', xlab = expression(xi), ylab = '', cex.lab = 1.5)
sd <- sd(xi.mu.nonstat)
mean <- mean(xi.mu.nonstat)
seq.xi <- seq(from =min(xi.mu.nonstat) - .1, to =max(xi.mu.nonstat) + .1, by = .01)
lines(seq.xi, dnorm(seq.xi, mean = mean, sd = sd), col = 'red')

par(oma = c(3,3,3,0))
title('Distribution of Parameters of Different Tide Stations With Mu Non-Stationary', outer = TRUE)
mtext(text = 'Parameter Values', side = 1, outer = TRUE)
mtext(text = 'Density', side = 2, outer = TRUE)

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

#save.image('east.coast.RData')

mu0.mu.sig.nonstat <- c(boston.priors$mu.sigma[1], 
                       atlantic.city.priors$mu.sigma[1], 
                       charleston.priors$mu.sigma[1],
                       montauk.priors$mu.sigma[1], 
                       chesapeake.bay.priors$mu.sigma[1], 
                       lewes.priors$mu.sigma[1], 
                       new.york.city.priors$mu.sigma[1], 
                       portland.priors$mu.sigma[1], 
                       duck.pier.priors$mu.sigma[1], 
                       newport.priors$mu.sigma[1])

mu1.mu.sig.nonstat <- c(boston.priors$mu.sigma[2], 
                        atlantic.city.priors$mu.sigma[2], 
                        charleston.priors$mu.sigma[2],
                        montauk.priors$mu.sigma[2], 
                        chesapeake.bay.priors$mu.sigma[2], 
                        lewes.priors$mu.sigma[2], 
                        new.york.city.priors$mu.sigma[2], 
                        portland.priors$mu.sigma[2], 
                        duck.pier.priors$mu.sigma[2], 
                        newport.priors$mu.sigma[2])

sig0.mu.sig.nonstat <- c(boston.priors$mu.sigma[3], 
                         atlantic.city.priors$mu.sigma[3], 
                         charleston.priors$mu.sigma[3],
                         montauk.priors$mu.sigma[3], 
                         chesapeake.bay.priors$mu.sigma[3], 
                         lewes.priors$mu.sigma[3], 
                         new.york.city.priors$mu.sigma[3], 
                         portland.priors$mu.sigma[3], 
                         duck.pier.priors$mu.sigma[3], 
                         newport.priors$mu.sigma[3])

sig1.mu.sig.nonstat <- c(boston.priors$mu.sigma[4], 
                         atlantic.city.priors$mu.sigma[4], 
                         charleston.priors$mu.sigma[4],
                         montauk.priors$mu.sigma[4], 
                         chesapeake.bay.priors$mu.sigma[4], 
                         lewes.priors$mu.sigma[4], 
                         new.york.city.priors$mu.sigma[4], 
                         portland.priors$mu.sigma[4], 
                         duck.pier.priors$mu.sigma[4], 
                         newport.priors$mu.sigma[4])

xi0.mu.sig.nonstat <- c(boston.priors$mu.sigma[5], 
                        atlantic.city.priors$mu.sigma[5], 
                        charleston.priors$mu.sigma[5],
                        montauk.priors$mu.sigma[5], 
                        chesapeake.bay.priors$mu.sigma[5], 
                        lewes.priors$mu.sigma[5], 
                        new.york.city.priors$mu.sigma[5], 
                        portland.priors$mu.sigma[5], 
                        duck.pier.priors$mu.sigma[5], 
                        newport.priors$mu.sigma[5])
#MEGA FIGURE----------------------------------------------
row1 <- c(1, 2, 3, 4, 5, 6)
row2 <- c(7, 8,9,10,11,12)
row3 <- c(13,14,15,16,17,18)
row4 <- c(19,20,21,22,23,24)

c<- rbind(row1, row2, row3, row4)
layout (c)
par(mar = c(.25,.25,.25,.25), oma = c(10,4,1,1))

#Stationary------------------------------------
par(mar = c(.25,1,2,1))
hist(mu.all.stat, freq = FALSE, xlab = '', ylab = '', cex.lab = 2, 
     yaxt = 'n',  cex.main = 2, cex.axis = 1.5, main = '')
#title(expresion(mu[0]))
mu.gamma <- fitdistr(mu.all.stat, 'gamma')
mu.seq <- seq(from = min (mu.all.stat) - 500, to = max(mu.all.stat) + 700, by = 1)
lines(mu.seq, dgamma(mu.seq, shape = mu.gamma$estimate[1], rate = mu.gamma$estimate[2]), col = 'red', lwd = 2)
mtext(text ='ST', side = 2, las = 1, line = 1, cex = 1.25)

plot.new()
#mtext(text = expression(mu[1]), line = 0, cex = 1.5)

hist(sigma.all.stat, freq = FALSE, xlab = expression(sigma), 
     ylab = '', cex.lab = 2, yaxt = 'n', cex.main = 2, cex.axis = 1.5, main = '')
sigma.gamma <- fitdistr(sigma.all.stat, 'gamma')
seq.sigma <- seq(from = min(sigma.all.stat)-200, to = max(sigma.all.stat)+700, by =1)
lines( seq.sigma, dgamma(seq.sigma, shape = sigma.gamma$estimate[1], rate = sigma.gamma$estimate[2] ), col = 'red', lwd = 2) 

plot.new()
#mtext(text = expression(sigma[1]), line = 0, cex = 1.5)

hist(xi.all.stat, freq = FALSE,  xlab = expression(xi), 
     ylab = '' ,cex.lab = 2, yaxt = 'n', cex.main = 2, cex.axis = 1.5, main = '')
sd <- sd(xi.all.stat)
mean <- mean(xi.all.stat)
seq.xi <- seq(from =min(xi.all.stat) - 5, to =max(xi.all.stat) + 2, by = .01)
lines(seq.xi, dnorm(seq.xi, mean = mean, sd = sd), col = 'red', lwd = 2)


plot.new()
#mtext(text = expression(xi[1]), line = 0, cex = 1.5)
#Mu NOn Stationary------------------------------------ 
hist(mu0.mu.nonstat, freq= FALSE, main = '', xlab = expression(mu[0]), ylab = '', cex.lab = 1.5, 
     yaxt = 'n', cex.axis = 1.5) 
mu0.gamma <- fitdistr(mu0.mu.nonstat, 'gamma')
mu0.seq <- seq(from = min (mu0.mu.nonstat) - 500, to = max(mu0.mu.nonstat) + 700, by = 1)
lines(mu0.seq, dgamma(mu0.seq, shape = mu0.gamma$estimate[1], rate = mu0.gamma$estimate[2]), col = 'red', lwd = 2)
mtext(text ='NS1', side = 2, las = 1, line = 1, cex = 1.25)

hist(mu1.mu.nonstat, freq= FALSE, main = '', xlab = expression(mu[1]), ylab = '', cex.lab = 1.5, 
     yaxt = 'n', cex.axis = 1.5)
sd <- sd(mu1.mu.nonstat)
mean <- mean(mu1.mu.nonstat)
mu1.seq <- seq(from = min (mu1.mu.nonstat) - 200, to = max(mu1.mu.nonstat) + 300, by = 1)
lines(mu1.seq, dnorm(mu1.seq, mean = mean, sd = sd), col = 'red', lwd = 2)

hist(sigma.mu.nonstat, freq= FALSE, main = '', xlab = expression(sigma), ylab = '', cex.lab = 1.5, 
     yaxt = 'n', cex.axis = 1.5)
sigma.gamma <- fitdistr(sigma.mu.nonstat, 'gamma')
sigma.seq <- seq(from = min (sigma.mu.nonstat) - 100, to = max(sigma.mu.nonstat) + 300, by = 1)
lines(sigma.seq, dgamma(sigma.seq, shape = sigma.gamma$estimate[1], rate = sigma.gamma$estimate[2]), col = 'red', lwd = 2)

plot.new()

hist(xi.mu.nonstat, freq= FALSE, main = '', xlab = expression(xi), ylab = '', cex.lab = 1.5, 
     yaxt = 'n', cex.axis = 1.5)
sd <- sd(xi.mu.nonstat)
mean <- mean(xi.mu.nonstat)
seq.xi <- seq(from =min(xi.mu.nonstat) - 1, to =max(xi.mu.nonstat) + 1, by = .01)
lines(seq.xi, dnorm(seq.xi, mean = mean, sd = sd), col = 'red', lwd = 2)

plot.new()

#------------------------------------ 
#mu, Sig, Non Stationary ------------------------ 
hist(mu0.mu.sig.nonstat, freq = FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mu0.gamma <- fitdistr(mu0.mu.sig.nonstat, 'gamma')
mu0.seq <- seq(from = min (mu0.mu.sig.nonstat) - 500, to = max(mu0.mu.sig.nonstat) + 300, by = 1)
lines(mu0.seq, dgamma(mu0.seq, shape = mu0.gamma$estimate[1], rate = mu0.gamma$estimate[2]), col = 'red', lwd = 2)
mtext(text = 'NS2', side = 2,  las = 1, line = 1, cex = 1.25)

hist(mu1.mu.sig.nonstat, freq = FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(mu1.mu.sig.nonstat)
sd <- sd(mu1.mu.sig.nonstat)
mu1.seq <- seq(from = min (mu1.mu.sig.nonstat)- 150, to = max(mu1.mu.sig.nonstat)+900, by = 1 )
lines(mu1.seq, dnorm(mu1.seq, mean = mean, sd = sd), col = 'red', lwd = 2)

hist(sig0.mu.sig.nonstat, freq = FALSE, main = '', yaxt = 'n', cex.axis = 1.5, ylim = c(0,.02))
sig0.gamma <- fitdistr(sig0.mu.sig.nonstat, 'gamma')
sig0.seq <- seq(from = min (sig0.mu.sig.nonstat) - 100, to = max(sig0.mu.sig.nonstat) + 300, by = 1)
lines(sig0.seq, dgamma(sig0.seq, shape = sig0.gamma$estimate[1], rate = sig0.gamma$estimate[2]), col = 'red', lwd = 2)

hist(sig1.mu.sig.nonstat, freq = FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(sig1.mu.sig.nonstat)
sd <- sd(sig1.mu.sig.nonstat)
sig1.seq <- seq(from = min(sig1.mu.sig.nonstat) - 10, to = max(sig1.mu.sig.nonstat) + 20, by = 1)
lines(sig1.seq, dnorm(sig1.seq, mean = mean, sd = sd), col = 'red', lwd = 2)

hist(xi0.mu.sig.nonstat, freq = FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(xi0.mu.sig.nonstat)
sd <- sd(xi0.mu.sig.nonstat)
xi0.seq <- seq(from = min (xi0.mu.sig.nonstat) - 2, to = max(xi0.mu.sig.nonstat) + .5, by = .05)
lines(xi0.seq, dnorm (xi0.seq , mean = mean, sd = sd), col = 'red', lwd = 2)

plot.new()

#------------------------------------
#All Non Stationary ------------------ 
hist(mu0.all.nonstat, freq= FALSE, main = '', yaxt = 'n', cex.axis = 1.5) 
mu0.gamma <- fitdistr(mu0.mu.nonstat, 'gamma')
mu0.seq <- seq(from = min (mu0.mu.nonstat) - 500, to = max(mu0.mu.nonstat) + 500, by = 1)
lines(mu0.seq, dgamma(mu0.seq, shape = mu0.gamma$estimate[1], rate = mu0.gamma$estimate[2]), col = 'red', lwd = 2)
mtext(text ='NS3', side = 2,  las = 1, line = 1, cex = 1.25)
mtext(text=expression(mu[0]), cex = 2, side = 1, line = 4)

hist(mu1.all.nonstat, freq= FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(mu1.all.nonstat)
sd <- sd(mu1.all.nonstat)
mu1.all.nonstat.seq <- seq(from = min(mu1.all.nonstat) -150, to = max(mu1.all.nonstat) + 900 , by = 1)
lines(mu1.all.nonstat.seq, dnorm(mu1.all.nonstat.seq, mean = mean, sd = sd), col = 'red', lwd = 2)
mtext(text=expression(mu[1]), cex = 2, side = 1, line = 4)

hist(sigma0.all.nonstat, freq= FALSE, main = '', yaxt = 'n', cex.axis = 1.5, xlim = c(3.5, 5.5), breaks =7)
sigma0.gamma <- fitdistr(sigma0.all.nonstat, 'gamma')
sigma0.seq <- seq(from = min (sigma0.all.nonstat) - 30, to = max(sigma0.all.nonstat) + 10, by = .05)
lines(sigma0.seq, dgamma(sigma0.seq, shape = sigma0.gamma$estimate[1], rate = sigma0.gamma$estimate[2]), col = 'red', lwd = 2)
mtext(text=expression(sigma[0]), cex = 2, side = 1, line = 4)

hist(sigma1.all.nonstat, freq= FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(sigma1.all.nonstat)
sd <- sd(sigma1.all.nonstat)
sigma1.seq <- seq(from = min (sigma1.all.nonstat) - 50, to = max(sigma1.all.nonstat) + 10, by = .05)
lines(sigma1.seq, dnorm(sigma1.seq, mean = mean, sd = sd), col = 'red', lwd = 2)
mtext(text=expression(sigma[1]), cex = 2, side = 1, line = 4)

hist(xi0.all.nonstat, freq= FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(xi0.all.nonstat)
sd <- sd(xi0.all.nonstat)
xi0.seq <- seq(from = min(xi0.all.nonstat) - 5, to = max(xi0.all.nonstat) + 3, by = .05)
lines(xi0.seq, dnorm(xi0.seq, mean = mean, sd = sd), col = 'red', lwd = 2)
mtext(text=expression(xi[0]), cex = 2, side = 1, line = 4)

hist(xi1.all.nonstat, freq= FALSE, main = '', yaxt = 'n', cex.axis = 1.5)
mean <- mean(xi1.all.nonstat)
sd <- sd(xi1.all.nonstat)
xi1.seq <- seq(from = min (xi1.all.nonstat) - .2, to = max(xi1.all.nonstat) + .2, by = .05)
lines(xi1.seq, dnorm(xi1.seq, mean = mean, sd = sd), col = 'red', lwd = 2)
mtext(text=expression(xi[1]), cex = 2,side = 1, line = 4)

#mtext(text = 'MLE Distributions of Different Tide Stations', outer = TRUE, cex = 1.5, line = 1)
mtext(text = 'Parameter Values', side = 1, outer = TRUE, line = 6, cex = 1.25)
#------------------------------------