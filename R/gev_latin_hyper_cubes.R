# Alex Klufas
# gev_latin_hyper_cubes.R
# Written on June 14, 2017 
# Modified on June 14, 2017
#
# Script in order to create Latin Hyper Cubes with data from New London Tide Guage (read_tide_data.R)

#creating latin hyper cube

#function to create lhs with certain size, number params will always be 3 (b/c GEV), upper and lower bounds of params in vector form , returns LHS of size size with 
#liklihood function calculated 

create.lhs.value <- function(size, upper.bound, lower.bound){

  lhs.cube <- randomLHS(n=size,k=3)
  
  #shift by 3000 and create new variable , shape shifted by 5
  lhs.cube.update <- cbind(((lhs.cube[,1]*(upper.bound - lower.bound)) -lower.bound ), 
                           (lhs.cube[,2]*(upper.bound - lower.bound) -lower.bound), 
                           (lhs.cube[,3](upper.bound - lower.bound) -lower.bound))

  #create even newer array to store all of this stuff / data 
  lhs.all.params <- rep(0,size)

  niter <-size
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for (q in 1:niter){
    lhs.all.params[q] <- sum(devd(lsl.max, loc=lhs.cube.update[q,1], scale=lhs.cube.update[q,2], shape=lhs.cube.update[q,3], log=TRUE, type=c("GEV")))
    setTxtProgressBar(pb, q)
  }
  close(pb)
  
  lhs.out <- vector('list', 2) #lhs, likelihood of lhs sample 
  names(lhs.out) <- c('lhs', 'likelihood') 
  lhs.out$lhs <- lhs.cube.update
  lhs.out$likelihood <- lhs.all.params
  
  return(lhs.out)
}


#creates a function that returns the plots of location v. scale, scale v. shape, location v. shape , returns the plots? 
create.lhs.plots <- function(lhs, lhs.likelihood.vals){
  
  #--------------------LOCATION v SCALE - 1000000 LHS --------------------------------
  lhs.many.values <- cbind(lhs[,1], lhs[,2], lhs.likelihood.vals)
  #labelling the matrix data 
  colnames(lhs.many.values) <- c('location', 'scale', 'likelihood')

  #creating empty raster <- still not really sure what this does 
  empty <- extent(lhs.many.values[,1:2])

  next.step5 <- raster(empty, ncol=100, nrow=100)
  last.step5 <- rasterize(lhs.many.values[,1:2], next.step5, lhs.many.values[,3], fun=mean)
  plot(last.step5, xlab='location', ylab='scale')
  title("location v. scale, z axis = likelihood")

  #log liklelihood
  plot(lhs[seq(from =1, to =1e6, by =1000),1], lhs[seq(from =1, to =1e6, by =1000)], ylim=c(-1000,0))

  plot(lhs[seq(from =1, to =1e6, by =1000),1], exp(lhs.all.params2[seq(from =1, to =1e6, by =1000)]))

  #------------------- SCALE v. SHAPE - 1000000 LHS --------------------------------


  lhs.many.values2 <- cbind(lhs[,2], lhs[,3], lhs.likelihood.vals)
  colnames(lhs.many.values2) <- c('scale', 'shape', 'likelihood')
  empty2 <- extent(lhs.many.values2[,1:2])

  next.step6 <- raster(empty2, ncol=1000, nrow=10)
  last.step6 <- rasterize(lhs.many.values2[,1:2], next.step6, lhs.many.values2[,3], fun=mean)
  plot(last.step6, xlab='scale', ylab='shape', ylim=c(0,50), xlim=c(1000,2000) )
  title("scale v. shape, z axis = likelihood")

  #-----------------LOCATION V. SHAPE - 1000000 LHS--------------------------------
  lhs.many.values3 <- cbind(lhs[,1], lhs[,3], lhs.likelihood.vals)
  colnames(lhs.many.values3) <- c('scale', 'shape', 'likelihood')
  empty3 <- extent(lhs.many.values3[,1:2])

  next.step7 <- raster(empty3, ncol=1000, nrow=10)
  last.step7 <- rasterize(lhs.many.values3[,1:2], next.step7, lhs.many.values3[,3], fun=mean)
  plot(last.step7, xlab='location', ylab='shape', ylim=c(0,10), xlim=c(1400,1450))
  title("location v. shape, z axis = likelihood")

  #order the likelihood in order to find the most likely values of certain parameters 
  lhs.sorted <- rev(order(lhs.likelihood.vals))

  #taking top 10 percent of all of this, gives us all of the top , these are just the indicies though
  lhs.top1p <- round(length(lhs.sorted)*.01)

  #this is list of all of the indicies that we need 
  ind.to.use <- lhs.sorted[1:lhs.top1p]

  #plotting location and scale of top 1 percent 
  plot(lhs[ind.to.use,1], lhs[ind.to.use,2], xlab = 'location', ylab = 'scale')
  title('location v. scale')
  plot(lhs[ind.to.use,2], lhs[ind.to.use,3], xlab = 'scale', ylab ='shape')
  title('scale v. shape')
  plot(lhs[ind.to.use,1], lhs[ind.to.use,3], xlab = 'location', ylab = 'scale')
  title('location v. shape')
  
  output <- vector('list', 4) #return indcies to use, sorted list, 
  names(output) <- c('ind.top1p', 'sorted')
  
  output$ind.top1p <- ind.to.use
  output$sorted <- lhs.sorted
  
  return(output)
}

#sf plots using lhs, idk what this would return though 
sf.using.lhs <- function(sorted.vals, lhs, data){
  #-------trying to create survival function plot with Latin Hyper Cubes ----------------------------------------------------------------

  #first lets find the best fitting parameters of the entire data set - aka the highest likelihood

  ind.high.like <- sorted.vals[1]

  colnames(lhs) <- c('location', 'scale', 'shape')
  
  sf.loc1 <- lhs[ind.high.like,1]
  sf.scale1 <- lhs[ind.high.like,2]
  sf.shape1 <- lhs[ind.high.like,3]

  x.hgt <- seq (0,3000, by=10)
  #list of indicies of max sea level sorted from lowest to highest 
  lsl.sort <- order(data)

  lsl.sorted.vals <- data[lsl.sort]

  #with the location, scale and shape that was given to us with the GEV curve parameters way above (around line 200-210)
  sf.hgt <- 1-pevd(x.hgt, loc=gev.mle$results$par[1], scale=gev.mle$results$par[2], shape = gev.mle$results$par[3], type=c("GEV"), lower.tail=TRUE)

  plot(x.hgt, log10(sf.hgt), type='l')

  sf.hgt2 <- 1-pevd(x.hgt, loc=sf.loc1, scale=sf.scale1, shape = sf.shape1, type=c("GEV"), lower.tail=TRUE)

  plot(x.hgt, log10(sf.hgt2), type='l')

  esf.vals <- seq(from=(length(lsl.sorted.vals)), to=1, by=(-1)) / (length(lsl.sorted.vals)+1)
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT')

  points(lsl.sorted.vals, log10(esf.vals)) 
  lines(x.hgt, log10(sf.hgt2), col='red')

  sf.top10 <- lhs.sorted[1:10]

  sf.location.top10 <- lhs.cube.update2[sf.top10, 1]
  sf.scale.top10 <- lhs.cube.update2[sf.top10, 2]
  sf.shapes.top10 <-lhs.cube.update2[sf.top10,3]

  #------Top 10 Highest Likelihood Overall---------------------------------------------------------------- 
  n.lhs.sorted <- length(lhs.sorted)
  for (s in 1:10){
   sf.hgt.test <- 1-pevd(x.hgt, loc=sf.location.top10[s], scale=sf.scale.top10[s], shape=sf.shapes.top10[s], type=c("GEV"), lower.tail=TRUE)
   lines(x.hgt, log10(sf.hgt.test), col='blue')
  }
}
