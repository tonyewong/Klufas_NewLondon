# Alex Klufas 
#
# Created on 6/1/17
# Modified on 6/1/17

#starting with basics 

#importing all of the packages

install.packages('adaptMCMC')
install.packages('compiler')
install.packages('DEoptim')
install.packages('doParallel')
install.packages('fExtremes')
install.packages('fields')
install.packages('fMultivar')
install.packages('foreach')
install.packages('gplots')
install.packages('graphics')
install.packages('lhs')
install.packages('maps')
install.packages('methods')
install.packages('ncdf4')
install.packages('plotrix')
install.packages('pscl')
install.packages('RColorBrewer')
install.packages('sensitivity')
install.packages('sn')
install.packages('stats')

#also adding extRemes b/c not here 

install.packages('extRemes')

#and also ismev 
install.packages('ismev')
install.packages('lubridate')
install.packages('zoo')
install.packages('ncdf4')
install.packages('Bolstad')

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
#files.tg <- list.files(path=dat.dir,pattern=filetype)

setwd('~/codes/Klufas_NewLondon/')

files.tg <- list.files(path=dat.dir,pattern=filetype)

data <- read.csv(paste(dat.dir,files.tg[1],sep=''), header=TRUE, sep=septype)
if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
        data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
    }
}
##==============================================================================

#===============================================================================
# Storm surge
#===============================================================================
  
  library(extRemes)
  library(fExtremes)
  library(ismev)
  library(lubridate)
  library(zoo)
  
  # Get tide gauge data and prepare to analyze.
  
  years         <- data$Year
  years.unique  <- unique(years)
  n.years       <- length(years.unique)
  lsl.mean      <- rep(0,length(n.years))
  lsl.max       <- rep(0,length(n.years))
  data$lsl.norm <- rep(NA,length(years))
  
  #get all of the hours out of data 
  hours <- data$Hour
  #get the unique hours - expect there to be 24 (and there are!))
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
  
  #subtracting one to make it the same length as the differnece in hours 
  num <- 1
 # hours.strange <- list()
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
  days <- data$Day
  day.unique <- unique(days)
  n.days <- length(days)
  
  #also need to check to see if there is a month gap because that could also pose the same type of problem 
  
  months <- data$Month
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
  }
    else if(years[ind]- years[ind-1] > 1){
      print("big gap here**************************************************")
      print(ind)
    }
  #if month gap greater than 1, make note of it 
    else if ((months[ind] - months[ind-1]) >= 1 | (months[ind] - months[ind-1]) <= -1){
     print("month gap greater than or equal to 1")
      print(months[ind] - months[ind-1])
      #print(ind)
    }
    #when the days are not the same , need to check how far apart 
    else {
     print("no month gap, only day gap")
      print(days[ind]-days[ind-1])
    }
  }
  
  #now I gotta figure out how to organize it all or something 
  #ideas ideas...
  
  # 1. get unique years of data
  # 2. subtract annual block means
  # 3. get annual block maxima <-- these are what we fit GEV distribution to
  
  #started at 2 so that the data the data stuff would start at 1939 not 1938 
  for (tt in 2:n.years) {
    #print(years.unique[tt])
    ind.thisyear <- which(years==years.unique[tt])
    lsl.mean[tt] <- mean(data$Sea_Level[ind.thisyear])
    data$lsl.norm[ind.thisyear] <- data$Sea_Level[ind.thisyear] - lsl.mean[tt]
    lsl.max[tt] <- max(data$lsl.norm[ind.thisyear])
  }
  

fit <- lm(lsl.mean ~ years.unique) #linear model

library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)

# fit a preliminary maximum likelihood estimate

gev.mle <- fevd(coredata(lsl.max), type='GEV') # extRemes
gev.mle2 <- gev.fit(coredata(lsl.max), show = FALSE) 

print(gev.mle$results$par) #printing location, scale, shape for hist

#x.lsl <- seq(from=0, to=3000, by=50)

x.lsl <- seq(from=0, to=3000, by=1)

#creating histogram 
hist(lsl.max, freq=FALSE)

#curve that will be placed on top of hist 
curve <- devd(x.lsl, loc=1056.2090549, scale=164.5221833, shape=0.1833216, threshold=0, log=FALSE, type=c("GEV"))

lines(curve, col = 'red')

##=====
# likelihood function stuff

#creating steps of the liklihood function
loc.move <- seq(from=0, to=3000, by=1)

#how many steps there are total 
n.loc.move <- length(loc.move)

#creating array to be filled in a momnet with changing location values 
curve.loc <- rep(0, n.loc.move)

#devd(lsl.max, loc=0, scale=164.5221833, shape=0.1833216, threshold=0, log=TRUE, type=c("GEV"))

#creates all the likelihood functions for location movement 
for (i in 1:n.loc.move){
   curve.loc[i] <- sum(devd(lsl.max, loc=loc.move[i], scale=164.5221833, shape=0.1833216, log=FALSE, type=c("GEV")))
   
}
plot(curve.loc)

#closest I have gotten 
#for (i in 1:n.loc.move){
  #curve.loc[i] <- sum(devd(lsl.max, loc=loc.move[i], scale=164.5221833, shape=log10(0.1833216), log=FALSE, type=c("GEV")))
#  
#}
#plot(curve.loc)

#creates all the likelihood functions for scale movement 
curve.scale <- rep(0, n.loc.move)
for (i in 1:n.loc.move){
  curve.scale[i] <- sum(devd(lsl.max, loc=1056.2090549, scale=loc.move[i], shape=0.1833216, log=FALSE, type=c("GEV")))
  
}
plot(curve.scale)

#creates all the likelihood functions for shape movement 
curve.shape <- rep(0, n.loc.move)
for (i in 1:n.loc.move){
  curve.shape[i] <- sum(devd(lsl.max, loc=1056.2090549, scale=164.5221833, shape=loc.move[i], log=FALSE, type=c("GEV")))
  
}
plot(curve.shape)



#example(slice3D)
#some.surface <- isosurf3D(lsl.max, curve.loc, curve.shape )


neg.lik <- function(gev.parameters, lsl.data){
  # fill this in!
  return(neg.lik)
}


xsquared <- function(x){
  output <- x*x
  return(output)
}

param.init <- 2
optim.out <- optim(param.init, fn=xsquared)


#creating latin hyper cubes 
library(lhs)

lhs.loc.scale <- randomLHS(n=3000, k=2)
#location (column 1 ) and then scale (column 2)
lhs.loc.scale.update <- cbind(lhs.loc.scale[,1]*(3000), lhs.loc.scale[,2]*(3000))

#creating empty array to store data 
lhs.loc.scale.like <- rep(0, 3000)

#similar for loop as above, summing probability values with two changing parameters 
for (j in 1:3000){
    lhs.loc.scale.like[j] <- sum(devd(lsl.max, loc=lhs.loc.scale.update[j,1], scale=lhs.loc.scale.update[j,2], shape=0.1833216, log=FALSE, type=c("GEV")))
  }

#plotting probability of each of the different changes in parameters 
plot(lhs.loc.scale.like)

#location v. liklihood with changing parameters
plot(lhs.loc.scale.update[,1], lhs.loc.scale.like)

#scale v. liklihood with changing parameters
plot(lhs.loc.scale.update[,2], lhs.loc.scale.like) #THESE ACTUALLY LOOK RIGHT :) 

#creating array of liklihoods with three changing parameters 

#initalize latin hyper cube, 3000 steps, 3 columns 
lhs.cube <- randomLHS(n=3000,k=3)

#shift by 3000 and create new variable 

lhs.cube.update <- cbind(lhs.cube[,1]*3000, lhs.cube[,2]*3000, lhs.cube[,3]*10)

#create even newer array to store all of this stuff / data 
lhs.all.params <- rep(0,3000)

for (q in 1:3000){
  lhs.all.params[q] <- sum(devd(lsl.max, loc=lhs.cube.update[q,1], scale=lhs.cube.update[q,2], shape=lhs.cube.update[q,3], log=FALSE, type=c("GEV")))
}

#plotting all changing parameters
plot(lhs.all.params)

#plotting location v. lillihood all changing parameters
plot(lhs.cube.update[,1],lhs.all.params)

#plotting scale v. lillihood all changing parameters
plot(lhs.cube.update[,2],lhs.all.params)


#plotting shape  v. lillihood all changing parameters
plot(lhs.cube.update[,3],lhs.all.params)


hist(lhs.all.params)

#trying to create surface plot 
install.packages('plot3D')
library(plot3D)

#x<-lhs.loc.scale.update[,1]
#y<-lhs.loc.scale.update[,2]
#z<-lhs.loc.scale.like
x <- seq(0,10, by =1)
y <- seq(0,10, by=1)
z <- outer(x,y)

#vals <- xyz.coords(x, y, 1)

lhs.loc.like <- cbind(lhs.loc.scale.like, lhs.loc.scale[,1])


contour(x, y, z)

#this isn't doing what I want it to do , so I will comment it out
#image.maker <- function(coords, value){
 # N <- length(coords[,1])
 # image.out <- matrix(0, nrow = N, ncol = N)
 # A <- as.numeric(factor(coords[,1])) #coords[,1] <- as.numeric(factor(coords[,1]))
#  B <- as.numeric(factor(coords[,2])) #coords[,2] <- as.numeric(factor(coords[,2]))
 # for (i in 1:nrow(coords)){
 #   image.out[A[i], B[i]] <- value[i]
#  }
# print(dim(image.out))
 # return(image.out)
#}


#my.image <- image.maker(lhs.loc.scale.update, lhs.loc.scale.like)
#contour(my.image)


install.packages('ggplot2')
library(ggplot2)

likelihood <- data.frame(like <- lhs.loc.scale.like)
vals <- ggplot(likelihood$like....lhs.loc.scale.like, aes(x,y))
vals + geom_contour() 

qplot(lhs.loc.scale.like, aes(x, y, z))
vals + geom_contour() 
#persp(z=lhs.loc.like, theta=-30)


filled.contour(x = 10*1:nrow(volcano),y = 10*1:ncol(volcano),
               z = volcano, color.palette = terrain.colors,
               plot.title = title(main = "The Topography of Maunga Whau",
                                  xlab = "Meters North",ylab = "Meters West"),
               plot.axes = {axis(1, seq(100, 800, by = 100))
                 axis(2, seq(100, 600, by = 100))},
               key.title = title(main="Height\n(meters)"),
               key.axes = axis(4, seq(90, 190, by = 10)))


install.packages('raster')
s100 <- matrix(c(267573.9, 2633781, 213.29545, 262224.4, 2633781, 69.78261, 263742.7, 2633781, 51.21951, 259328.4, 2633781, 301.98413, 264109.8, 2633781, 141.72414, 255094.8, 2633781, 88.90244),  ncol=3,  byrow=TRUE)
colnames(s100) <- c('X', 'Y', 'Z')

library(raster)
# set up an 'empty' raster, here via an extent object derived from your data
e <- extent(s100[,1:2])
e <- e + 1000 # add this as all y's are the same

r <- raster(e, ncol=10, nrow=2)
# or r <- raster(xmn=, xmx=,  ...

# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(s100[, 1:2], r, s100[,3], fun=mean)
plot(x)
#my idea, first make sequenced grid, I think I need to interpolate the values tbh 

#-----------------contour type plot of mean (x axis), sigma (y axis), likelihood (z-axis) ,scale was kept constant ------------------------
#making matrix w three columns, mean, sigma, likelihood lhs.all.params

#using scale here, not shape, fix this post lunch_________________________________________________________________________
lhs.new.stuff <- cbind(lhs.loc.scale.update[,1], lhs.loc.scale.update[,2], lhs.loc.scale.like)
#labelling the matrix data 
colnames(lhs.new.stuff) <- c('loc', 'scale', 'likelihood')

library(raster)

#creating empty raster <- still not really sure what this does 
empty <- extent(lhs.new.stuff[,1:2])

next.step <- raster(empty, ncol=10, nrow=10)
last.step <- rasterize(lhs.new.stuff[,1:2], next.step, lhs.new.stuff[,3], fun=mean)
plot(last.step, xlab='location', ylab='scale')

## fun and nice things
t0 <- proc.time()
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for (i in 1:niter){
  # do stuff
  Sys.sleep(1)
  setTxtProgressBar(pb, i)
}
close(pb)
t1 <- proc.time()
time_that_took <- t1-t0
## end fun and nice things

#contoured plot of the above section!!!!!!!!!!!
#I dont undesrstand why it is not plotting below 164 ish (the place where shape supposed to end up at)
contour(last.step, ylab ='shape', xlab='location', nlevels=50, col='red', ylim = c(0,500), xlim=c(700,1500)) #only works with low resolution???  

# idea / thought - when I make it high resolution, it centers on the place where there is a high probability, so maybe I should think about zooming in?? (or
# something like that 

# lol okay what is my next step
# maybe I can get a surface out of this

#persp3D(x=1:10, y=1:10, z=last.step)

##-----------------contour type plot of mean (x axis), sigma (y axis), likelihood (z-axis) , shape was also changing 
lhs.new.stuff2 <- cbind(lhs.cube.update[,1], lhs.cube.update[,2], lhs.all.params)
#labelling the matrix data 
colnames(lhs.new.stuff2 ) <- c('loc', 'scale', 'likelihood')

library(raster)

#creating empty raster <- still not really sure what this does 
empty <- extent(lhs.new.stuff2[,1:2])

next.step2 <- raster(empty, ncol=10, nrow=10)
last.step2 <- rasterize(lhs.new.stuff2[,1:2], next.step2, lhs.new.stuff2[,3], fun=mean)
plot(last.step2, ylab='scale', xlab='location')

#contoured plot of the above section!!!!!!!!!!!
contour(last.step2, ylab ='scale', xlab='location',nlevels=50) #only works with low resolution??? 

#-----------------------contouring location, shpae, likelihood, scale also changing 

lhs.new.stuff3 <- cbind(lhs.cube.update[,1], lhs.cube.update[,3], lhs.all.params)
#labelling the matrix data 
colnames(lhs.new.stuff3 ) <- c('loc', 'shape', 'likelihood')

library(raster)

#creating empty raster <- still not really sure what this does 
empty <- extent(lhs.new.stuff3[,1:2])

next.step3 <- raster(empty, ncol=10, nrow=10)
last.step3 <- rasterize(lhs.new.stuff3[,1:2], next.step3, lhs.new.stuff3[,3], fun=mean)
plot(last.step3, ylab ='shape', xlab='location')

#contoured plot of the above section!!!!!!!!!!!
contour(last.step3, ylab ='shape', xlab='location', nlevels=50) #only works with low resolution??? 

#-----------------------contouring shape, scale, likelihood, location changing 

lhs.new.stuff4 <- cbind(lhs.cube.update[,2], lhs.cube.update[,3], lhs.all.params)
#labelling the matrix data 
colnames(lhs.new.stuff4 ) <- c('scale', 'shape', 'likelihood')

library(raster)

#creating empty raster <- still not really sure what this does 
empty <- extent(lhs.new.stuff4[,1:2])

next.step4 <- raster(empty, ncol=100, nrow=100)
last.step4 <- rasterize(lhs.new.stuff4[,1:2], next.step4, lhs.new.stuff4[,3], fun=mean)
plot(last.step4, xlab='scale', ylab='shape')

#contoured plot of the above section!!!!!!!!!!!
contour(last.step4, ylab ='shape', xlab='scale', nlevels=50) #only works with low resolution??? 


lhs.cube2 <- randomLHS(n=1000000,k=3)

#shift by 3000 and create new variable , shape shifted by 5

lhs.cube.update2 <- cbind(lhs.cube2[,1]*3000, lhs.cube2[,2]*3000, lhs.cube2[,3]*(10) - 5)

#create even newer array to store all of this stuff / data 
lhs.all.params2 <- rep(0,1000000)

niter <-1000000
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for (q in 1:niter){
  lhs.all.params2[q] <- sum(devd(lsl.max, loc=lhs.cube.update2[q,1], scale=lhs.cube.update2[q,2], shape=lhs.cube.update2[q,3], log=FALSE, type=c("GEV")))
  setTxtProgressBar(pb, q)
}
close(pb)

#--------------------LOCATION v SCALE - 1000000 LHS 
lhs.many.values <- cbind(lhs.cube.update2[,1], lhs.cube.update2[,2], lhs.all.params2)
#labelling the matrix data 
colnames(lhs.many.values) <- c('location', 'scale', 'likelihood')

library(raster)

#creating empty raster <- still not really sure what this does 
empty <- extent(lhs.many.values[,1:2])

next.step5 <- raster(empty, ncol=100, nrow=100)
last.step5 <- rasterize(lhs.many.values[,1:2], next.step5, lhs.many.values[,3], fun=mean)
plot(last.step5, xlab='location', ylab='scale')
title("location v. scale, z axis = likelihood")

#contoured plot of the above section!!!!!!!!!!!

#------------------- SCALE v. SHAPE - 1000000 LHS 

lhs.many.values2 <- cbind(lhs.cube.update2[,2], lhs.cube.update2[,3], lhs.all.params2)
colnames(lhs.many.values2) <- c('scale', 'shape', 'likelihood')
empty2 <- extent(lhs.many.values2[,1:2])

next.step6 <- raster(empty2, ncol=1000, nrow=10)
last.step6 <- rasterize(lhs.many.values2[,1:2], next.step6, lhs.many.values2[,3], fun=mean)
plot(last.step6, xlab='scale', ylab='shape', ylim=c(0,50), xlim=c(1000,2000) )
title("scale v. shape, z axis = likelihood")

#contoured plot of the above section!!!!!!!!!!!


#-----------------LOCATION V. SHAPE - 1000000 LHS
lhs.many.values3 <- cbind(lhs.cube.update2[,1], lhs.cube.update2[,3], lhs.all.params2)
colnames(lhs.many.values3) <- c('scale', 'shape', 'likelihood')
empty3 <- extent(lhs.many.values3[,1:2])

next.step7 <- raster(empty3, ncol=1000, nrow=10)
last.step7 <- rasterize(lhs.many.values3[,1:2], next.step7, lhs.many.values3[,3], fun=mean)
plot(last.step7, xlab='location', ylab='shape', ylim=c(0,10), xlim=c(1400,1450))
title("location v. shape, z axis = likelihood")

#contoured plot of the above section!!!!!!!!!!!

#x = seq(0,3000,10)
#y = seq(0,3000,10)
#x <- seq(0,54,1)
#y <- seq(0,54,1)
##z <- lhs.loc.scale.like

#order the likelihood in order to find the most likely values of certain parameters 
lhs.sorted <- rev(order(lhs.all.params2))

#taking top 10 percent of all of this, gives us all of the top , these are just the indicies though
lhs.top1p <- round(length(lhs.sorted)*.01)

#this is list of all of the indicies that we need 
ind.to.use <- lhs.sorted[1:lhs.top1p]

#plotting location and scale of top 1 percent - no colors yet 
plot(lhs.cube.update2[ind.to.use,1], lhs.cube.update2[ind.to.use,2])

colnames(lhs.cube.update2)
#need to figure out how to get indicies to take from other lists - ie from locaiton and scale 

#making list of colors - top 1 percent that work 
lhs.color.rep <- rep(0,lhs.top1p)

#assigning colors depending on the value of the shape parameter 
lhs.color.rep[lhs.cube.update2[ind.to.use,3] > 0& lhs.cube.update2[ind.to.use,3] <=3]= 'blue'
lhs.color.rep[lhs.cube.update2[ind.to.use,3] > -3& lhs.cube.update2[ind.to.use,3] <= 0]= 'black'
lhs.color.rep[lhs.cube.update2[ind.to.use,3] <= -3] = 'red'
lhs.color.rep[lhs.cube.update2[ind.to.use,3] >3] ='pink'

#plotting location v. scale with color representing shape value colors for top 1 percent 
plot(lhs.cube.update2[ind.to.use,1], lhs.cube.update2[ind.to.use,2], col=lhs.color.rep, xlab ='location', ylab ='scale')
title('location v. scale, top 1 percent')
legend(0,3000,c('>3 & <6', '<=3', '>=6'))
lty=c(1,1)
lwd=c(2.5,2.5, 2.5)
col=c('black','red','blue')
#want to plot location v. scale from this 

#trying to see what happens when we use the top .1% - looking to see if any pattern exists 
lhs.top01p <- round(length(lhs.sorted)*.001)
ind.to.use2 <- lhs.sorted[1:lhs.top01p]

lhs.color.rep2 <- rep(0,lhs.top01p)

#same color idea as above but for smaller percent of the data 
lhs.color.rep2[lhs.cube.update2[ind.to.use2,3] > 0& lhs.cube.update2[ind.to.use2,3] <=3]= 'blue'
lhs.color.rep2[lhs.cube.update2[ind.to.use2,3] > -3& lhs.cube.update2[ind.to.use2,3] <= 0]= 'black'
lhs.color.rep2[lhs.cube.update2[ind.to.use2,3] <= -3] = 'red'
lhs.color.rep2[lhs.cube.update2[ind.to.use2,3] >3] ='pink'

plot(lhs.cube.update2[ind.to.use2,1], lhs.cube.update2[ind.to.use2,2], col=lhs.color.rep2, xlab ='location', ylab ='scale')
title('location v. scale, top .1 percent')

#plot(lhs.many.values[,1], lhs.many.values[,2])
lhs.top10p <- round(length(lhs.sorted)*.1)
ind.to.use3 <- lhs.sorted[1:lhs.top10p]

lhs.color.rep3 <- rep(0,lhs.top10p)

#same color idea as above but for larger percent of the data
lhs.color.rep3[lhs.cube.update2[ind.to.use3,3] > 0& lhs.cube.update2[ind.to.use3,3] <=3]= 'blue'
lhs.color.rep3[lhs.cube.update2[ind.to.use3,3] > -3& lhs.cube.update2[ind.to.use3,3] <=0]= 'black'
lhs.color.rep3[lhs.cube.update2[ind.to.use3,3] <= -3] = 'red'
lhs.color.rep3[lhs.cube.update2[ind.to.use3,3] >3] ='pink'

plot(lhs.cube.update2[ind.to.use3,1], lhs.cube.update2[ind.to.use3,2], col=lhs.color.rep3, xlab ='location', ylab ='scale')
title('location v. scale, top 10 percent')


#-------trying to create survival function plot

#first lets find the best fitting parameters of the entire data set - aka the highest likelihood

ind.high.like <- lhs.sorted[1]

colnames(lhs.cube.update2) <- c('location', 'scale', 'shape')

sf.loc1 <- lhs.cube.update2[ind.high.like,1]
sf.scale1 <- lhs.cube.update2[ind.high.like,2]
sf.shape1 <- lhs.cube.update2[ind.high.like,3]

#I think I use pevd 

#I have lsl.max - how would I use it though 

x.hgt <- seq (0,3000, by=10)
#list of indicies of max sea level sorted from lowest to highest 
lsl.sort <- order(lsl.max)

lsl.sorted.vals <- lsl.max[lsl.sort]

#with the location, scale and shape that was given to us with the GEV curve parameters way above (around line 200-210)
sf.hgt <- 1-pevd(x.hgt, loc=1056.2090549, scale=164.5221833, shape = 0.1833216, type=c("GEV"), lower.tail=TRUE)

sf.hgt.data <- 1-pevd(lsl.max, loc=1056.2090549, scale=164.5221833, shape = 0.1833216, type=c("GEV"), lower.tail=TRUE)
plot(x.hgt, sf.hgt, log='y', type='l')
plot(x.hgt, log10(sf.hgt), type='l')

sf.hgt2 <- 1-pevd(x.hgt, loc=sf.loc1, scale=sf.scale1, shape = sf.shape1, type=c("GEV"), lower.tail=TRUE)

plot(x.hgt, sf.hgt2, log='y', type='l')
plot(x.hgt, log10(sf.hgt2), type='l')

#points(lsl.max)

esf.vals <- seq(from=(length(lsl.sorted.vals)), to=1, by=(-1)) / (length(lsl.sorted.vals)+1)
plot(x.hgt, log10(sf.hgt), type='l')
points(lsl.sorted.vals, log10(esf.vals)) 

#-------TAKING OUT ALL OF THE WEIRD SHAPE PARAMETERS 

#need to find max of the sea levels expereinced 

max.sea.hgt <- lsl.max[rev(order(lsl.max))[1]]

#max.sea.hgt is the standard - if the upper bound of some combination of location, scale and shape are below this - then we have to pull it out of the latin hyper cube
lhs.cube.refined <- lhs.cube.update2

# (location -  scale) / shape 

niter <-10000
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for (z in 1:niter){
  if(((lhs.cube.refined[z,1] - lhs.cube.refined[z,2]) / lhs.cube.refined[z,3]) < max.sea.hgt){
    lhs.cube.refined[-z,]
  }
  setTxtProgressBar(pb, z)
}
close(pb)


##==============================================================================
## End
##==============================================================================
