# Alex Klufas
# read_tide_data.R
# Written on June 14, 2017 
# Modified on July 27, 2017
#
# Script in order to read tide data, takes out years that are missing more than 10% of data
# Need to already have column names as part of the data file (Year, Month, Day, Hour, Sea_Level)


read.tide.data <- function (city){
  
  #set up vector for pulling the data 
  city.data <- vector('list', 3)
  names(city.data) <- c('max', 'mean','years')
  years         <- city$Year
  years.unique  <- unique(years)
  n.years       <- length(years.unique) 
  lsl.mean      <- rep(0,length(n.years))
  lsl.max       <- rep(0,length(n.years)) 
  city$lsl.norm <- rep(NA,length(years))
  
  
  
  #get all of the hours out of data 
  hours <- city$Hour
  hours.unique <- unique(hours)
  n.hours <- length(hours)
  
  #finding the difference in hours
  hour.diff <- diff(hours)
  #cheking to see if any places dont have numbers
  for (i in 1:n.hours){
    if (is.nan(hour.diff[i])){
      print("Hour Value Does Not Exist, is Nan at Index:")
      print(i)
    }
  }
  
  #the hours that are weird
  day.counter <- 0
  month.counter <- 0
  #subtracting one to make it the same length as the differnece in hours 
  num <- 1
  hours.strange <- rep(0, length(n.hours)) 
  
  for (i in 1:(n.hours-1)){
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
  months <- city$Month
  n.month <- length(months)
  
  
  #look at each index and the one after and check to see if days are within a day of each other 
  #if so - check hour difference 
  #else - check if in different months - cases -31, -30, -28, -29 (for feb)
  #go through the list of weird numbers we have 
  for (j in 1:length(hours.strange)){
    #need to check within the same day of each other - so that would be seeing if the difference between the days is 0 
    #index we will be using comes from the list that we just made (the list hours.strange)
    ind <- hours.strange[[j]]
    
    #if the difference between the two data pts is in the same day - ignore it 
    if ((days[ind]-days[ind-1] ) == 0 & (months[ind] - months[ind-1]) == 0){
      print("Gap in the same day at index:")
      print(ind)
    }else if(years[ind]- years[ind-1] > 1){
      print("Year gap or greater at index:")
      print(ind)
    }else if ((months[ind] - months[ind-1]) >= 1 | (months[ind] - months[ind-1]) <= -1) { #if month gap greater than 1, make note of it
      if (abs(years[ind]- years[ind-1] > 1)){
        print("More than one month gap at index: ")
        print(ind)
        print("Length of gap: ")
        print(abs(months[ind] - months[ind-1]))
        month.counter <- month.counter + abs(months[ind] - months[ind-1])
      }}
    #when the days are not the same , need to check how far apart 
    else {
      print("Only month gap at index: ")
      print(ind)
      print("Length of gap: ")
      print(days[ind]-days[ind-1])
      day.counter <- day.counter + days[ind]-days[ind-1]
      print("Total Days Missing: ")
      print(day.counter)
    }
  }
  
  ##----------Finding Annual Block Maxima----------------------
  for (tt in 1:n.years) {
    ind.thisyear <- which(years==years.unique[tt])
    #checking if more than 90% of the year's data is in the year 
    if (length(ind.thisyear) >= 7884){ #90% of number of hours in a year 
      lsl.mean[tt] <- mean(city$Sea_Level[ind.thisyear])
      city$lsl.norm[ind.thisyear] <- city$Sea_Level[ind.thisyear] - lsl.mean[tt]
      lsl.max[tt] <- max(city$lsl.norm[ind.thisyear])
    }
    else{
      lsl.mean[tt] <- 0 
      city$lsl.norm[ind.thisyear] <- 0 
      lsl.max[tt] <- 0
    }
  }
  
  #removing the zeros 
  ind <- which(lsl.max != 0)
  years.shift <- years.unique[ind]
  sl <- lsl.max[ind]
  sl.mean <- lsl.mean[ind]
  
  city.data$years <- years.shift 
  city.data$max <- sl
  city.data$mean <- sl.mean
  
  return(city.data)
}