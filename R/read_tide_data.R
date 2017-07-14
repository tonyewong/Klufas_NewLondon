# Alex Klufas
# read_tide_data.R
# Written on June 14, 2017 
# Modified on June 14, 2017
#
# Script in order to read New London tide guage data 

read.tide.data <- function(){
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
  #library(fExtremes)
  #library(extRemes)
  #library(lubridate)
  #library(zoo)
  
  # Get tide gauge data and prepare to analyze.
  
  
  
  years         <- data$Year
  years.unique  <- unique(years)
  n.years       <- length(years.unique) 
  lsl.mean      <- rep(0,length(n.years))
  lsl.max       <- rep(0,length(n.years))
  data$lsl.norm <- rep(NA,length(years))
  
  
  ##----------Finding Annual Block Maxima----------------------
  #started at 2 so that the data the data stuff would start at 1939 not 1938 
  for (tt in 1:n.years) {
    ind.thisyear <- which(years==years.unique[tt])
    lsl.mean[tt] <- mean(data$Sea_Level[ind.thisyear])
    data$lsl.norm[ind.thisyear] <- data$Sea_Level[ind.thisyear] - lsl.mean[tt]
    lsl.max[tt] <- max(data$lsl.norm[ind.thisyear])
  }
  
  lsl.max <- lsl.max[2: n.years]
  
  years.unique <- years.unique[2: n.years]
  
  lsl <- vector('list', 3)
  names(lsl) <- c('max', 'mean','years')
  lsl$max <- lsl.max 
  lsl$mean <- lsl.mean
  lsl$years <- years.unique
  
  return(lsl)
}
