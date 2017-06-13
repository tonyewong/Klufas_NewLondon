# Alex Klufas
# read_temp_data.R
# Written on June 14, 2017 
# Modified on June 14, 2017
#
# Script in order to read NOAA temperature data 

#sets working directory for script
setwd('~/codes/Klufas_NewLondon/')

#read temp data 
temp.data <- read.table('noaa-temp-1880-2017.csv', header = TRUE, sep=',')
temp.years <- temp.data$Year
temp.values <- temp.data$Value[59:134]     #limit temp data to just the years I am working with 
