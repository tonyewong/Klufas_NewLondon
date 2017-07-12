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
ChesapeakeBay <-read.csv(paste(dat.dir,files.tg[3],sep=''), header=TRUE, sep=septype)
colnames(ChesapeakeBay) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Lewes <-read.csv(paste(dat.dir,files.tg[4],sep=''), header=TRUE, sep=septype)
colnames(Lewes) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
NewYorkCity <-read.csv(paste(dat.dir,files.tg[5],sep=''), header=TRUE, sep=septype)
colnames(NewYorkCity) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')
Portland <- read.csv(paste(dat.dir,files.tg[6],sep=''), header=TRUE, sep=septype)
colnames(Portland) <- c('Year', 'Month', 'Day', 'Hour', 'Sea_Level')

city.list <- c(AtlanticCity,Boston,ChesapeakeBay,Lewes,NewYorkCity , Portland) 

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