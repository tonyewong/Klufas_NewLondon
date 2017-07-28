# Alex Klufas
# read_temp_data.R
# Written on June 14, 2017
# Modified on June 14, 2017
#
# Script in order to read NOAA temperature data
# Before using this, must have read the file in to take out years that have less than 10% of the available
# data 


read.temp.data <- function(city.years){
  #sets working directory for script
  setwd('~/codes/Klufas_NewLondon/')

  temp.data <- read.table('noaa-temp-1880-2017.csv', header = TRUE, sep=',')
  temp.years <- rep(0,length(city.years))
  temp.values <- rep(0,length(city.years))     #limit temp data to just the years I am working with

  for(i in 1:length(city.years)){
    ind <- which(city.years[i] == temp.data$Year)
    temp.years[i] <- temp.data$Year[[ind]]
    temp.values[i] <- temp.data$Value[[ind]]
  }

  temps <- vector('list', 2)

  names(temps) <- c('years','values')

  temps$years <- temp.years
  temps$values <- temp.values

  return(temps)
}
