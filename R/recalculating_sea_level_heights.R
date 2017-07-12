
#calling csv file
dat.dir <- './data/'
filetype <- 'csv'
septype <- ','

#naming code 
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.projout <- paste('../output_model/BRICK_project-lsl-surge_NewLondon_',today,'.nc',sep='')
filename.lslout  <- paste('../output_model/BRICK_project-lsl_NewLondon_',today,'.csv', sep="")

##==============================================================================
## Read tide gauge data

setwd('~/codes/Klufas_NewLondon/')

files.tg <- list.files(path=dat.dir,pattern=filetype)

data <- read.csv(paste(dat.dir,files.tg[1],sep=''), header=TRUE, sep=septype)
if(length(files.tg) > 1) {
  for (ff in 2:length(files.tg)) {
    data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
  }
}

#===============================================================================
# Storm surge
#===============================================================================

# Get tide gauge data and prepare to analyze.

years         <- data$Year
years.unique  <- unique(years)
n.years       <- length(years.unique) 
lsl.mean      <- rep(0,length(n.years))
lsl.max       <- rep(0,length(n.years)) #max(data$lsl.norm[ind.thisyear])
data$lsl.norm <- rep(NA,length(years))
lsl.max2      <- rep(0, length(n.years))
fit <- lm(data$Sea_Level ~ c(1:length(data$Sea_Level)))

##----------Finding Annual Block Maxima----------------------
#started at 2 so that the data the data stuff would start at 1939 not 1938 
for (tt in 1:n.years) {
  ind.thisyear <- which(years==years.unique[tt])
  lsl.mean[tt] <- mean(data$Sea_Level[ind.thisyear])
  data$lsl.norm[ind.thisyear] <- data$Sea_Level[ind.thisyear] - lsl.mean[tt]
  lsl.max[tt] <- max(data$lsl.norm[ind.thisyear])
  lsl.max2[tt] <- max(fit$residuals[ind.thisyear]) 
}

#water level anomaly = verified water level measurement - [sea level trend(date - date 
#zero SLR)] - predicted water level measurement 

#SLTi <- just the sea level trend, when they did the paper it was 2.13mm / year 
# I will use current which is 2.57 mm / year 

SLTi <- 2.57    #mm/day 
#need to figure out hours 

WLA <- rep(0, length(n.years))

# i is not important here - normally would refer to the station being used, but here
#only wokring w one station 
MSL <- -.092 

#WLV <- 
  
Toi <- 44 #currently storing as index to pull later but don't really know how to use this

for (i in 1:n.years){
  ind.thisyear <- which(years==years.unique[i])
  lsl.mean[i] <- mean(data$Sea_Level[ind.thisyear])
  data$lsl.norm[ind.thisyear] <- data$Sea_Level[ind.thisyear] - lsl.mean[i]
  lsl.max[i] <- max(data$lsl.norm[ind.thisyear])
  lsl.max2[i] <- max(fit$residuals[ind.thisyear]) 
  WLA[i] <- lsl.max[i] - SLTi*(i-44) - 400
}

#fit <- lm(data$Sea_Level ~ c(1:length(data$Sea_Level)))
#plot(years, data$SeaLevel)

months <- data$Month
months.per.yr <- rep(0,n.years)
months.unique <- unique(months)
n.months <- length(months.unique)
#for (j in 1:n.years){
 # for (k in 1:n.months){
#  months.per.yr[j] <- 
  #}
 # months.per.yr[j] <- length(unique(months) && n.years[j])
#}

lsl.max.2  <- rep(0,length(n.years))
#going to try and find monthly block max to see if that helps 
for (i in 1:n.years){
  #ind.month <- which(years==years.unique[i] && months==months.unique[i]) #something about only taking the months that are in this specific year 
  #n.months.per.yr <- which(months==)
  ind.thisyear <- which(years==years.unique[i])
  month.max <- rep(0, length(months.unique))
  lsl.mean <- rep(0,n.months)
  for (q in 1:n.months){
    ind.thismonth <- which(months==months.unique[q] & years==years.unique[i])
    lsl.mean[q] <- mean(data$Sea_Level[ind.thismonth])
    data$lsl.norm[ind.thismonth] <- data$Sea_Level[ind.thismonth] - lsl.mean[q]
    month.max[q] <- data$lsl.norm[ind.thismonth]
  }
  lsl.max.2[i] <- max(month.max)
}

lsl.max <- lsl.max[2: n.years]

years.unique <- years.unique[2: n.years]


fit <- lm(lsl.max ~ years.unique) 
