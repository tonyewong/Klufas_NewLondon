# Alex Klufas
# tide_guage_driver.R
# Written on June 16, 2017 
# Modified on June 16, 2017

#-------getting tide guage data----------------------------
source('read_tide_data.R')
data <- lsl.max 

#-------getting temperature data---------------------------
source('read_temp_data.R')
temps <- temp.values

#------creating survival function--------------------------