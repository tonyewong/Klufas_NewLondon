# Alex Klufas
# tide_guage_driver.R
# Written on June 16, 2017 
# Modified on June 16, 2017

install.packages('fExtremes')
install.packages('extRemes')
install.packages('lubridate')
install.packages('zoo')

#-------libraries needed to run code-----------------------
library(fExtremes)
library(extRemes)
library(lubridate)
library(zoo)
library(lhs)
library(raster)
library(DEoptim)
library(lhs)
library(raster)


setwd('~/codes/Klufas_NewLondon/R/')

#-------getting tide guage data----------------------------
source('read_tide_data.R')
data <- read.tide.data()

#-------getting temperature data---------------------------
source('read_temp_data.R')
temps <- read.temp.data(1939,2014)

#------creating survival function--------------------------
source('optimization_sf.R')
