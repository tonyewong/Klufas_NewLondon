  # Alex Klufas 
  #
  # Created on 6/1/17
  # Modified on 6/1/17
  
  #starting with basics 
  
  #importing all of the packages
  
  install.packages('DEoptim')
  install.packages('fExtremes')
  install.packages('lhs')
  install.packages('ncdf4')
  install.packages('extRemes')
  install.packages('ismev')
  install.packages('lubridate')
  install.packages('zoo')
  install.packages('ncdf4')
  install.packages('raster')
  install.packages('DEoptim')
  
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
    library(fExtremes)
    library(extRemes)
    library(lubridate)
    library(zoo)
    
    # Get tide gauge data and prepare to analyze.
    
    years         <- data$Year
    years.unique  <- unique(years)
    n.years       <- length(years.unique) 
    lsl.mean      <- rep(0,length(n.years))
    lsl.max       <- rep(max(data$lsl.norm[ind.thisyear]),length(n.years))
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
    
  
  fit <- lm(lsl.mean ~ years.unique) #linear model
  
  # fit a preliminary maximum likelihood estimate
  library(ismev)
  
  gev.mle <- fevd(coredata(lsl.max), type='GEV') # extRemes
  gev.mle2 <- gev.fit(coredata(lsl.max), show = FALSE) 
  
  print(gev.mle$results$par) #printing location, scale, shape for hist
  
  x.lsl <- seq(from=0, to=3000, by=1)
  
  #creating histogram 
  hist(lsl.max, freq=FALSE)
  
  #curve that will be placed on top of hist 
  curve <- devd(x.lsl, loc=gev.mle$results$par[1], scale=gev.mle$results$par[2], shape=gev.mle$results$par[3], threshold=0, log=FALSE, type=c("GEV"))
  lines(curve, col = 'red')
  
  #creating steps of the liklihood function
  loc.move <- seq(from=0, to=3000, by=1)
  
  #how many steps there are total 
  n.loc.move <- length(loc.move)
  
  #creating array to be filled in a momnet with changing location values 
  curve.loc <- rep(0, n.loc.move)
  
  #creates all the likelihood functions for location movement 
  for (i in 1:n.loc.move){
     curve.loc[i] <- sum(devd(lsl.max, loc=loc.move[i], scale=gev.mle$results$par[2], shape=gev.mle$results$par[3], log=FALSE, type=c("GEV")))
     
  }
  plot(curve.loc)
  
  #creates all the likelihood functions for scale movement 
  curve.scale <- rep(0, n.loc.move)
  for (i in 1:n.loc.move){
    curve.scale[i] <- sum(devd(lsl.max, loc=gev.mle$results$par[1], scale=loc.move[i], shape=gev.mle$results$par[3], log=FALSE, type=c("GEV")))
    
  }
  plot(curve.scale)
  
  shape.move <- seq(from=-2, to=2, by=.02)
  
  n.shape.move <- length(shape.move)
  #creates all the likelihood functions for shape movement 
  curve.shape <- rep(0, n.shape.move)
  for (i in 1:n.shape.move){
    curve.shape[i] <- prod(devd(lsl.max, loc=gev.mle$results$par[1], scale=gev.mle$results$par[2], shape=shape.move[i], log=FALSE, type=c("GEV")))
    
  }
  plot(shape.move, log(curve.shape))
  
  #creating latin hyper cubes 
  library(lhs)
  library(raster)
  
  lhs.cube2 <- randomLHS(n=1000000,k=3)
  
  #shift by 3000 and create new variable , shape shifted by 5
  lhs.cube.update2 <- cbind(lhs.cube2[,1]*3000, lhs.cube2[,2]*3000, lhs.cube2[,3]*(10) - 5)
  
  #create even newer array to store all of this stuff / data 
  lhs.all.params2 <- rep(0,1000000)
  
  niter <-1000000
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for (q in 1:niter){
    lhs.all.params2[q] <- sum(devd(lsl.max, loc=lhs.cube.update2[q,1], scale=lhs.cube.update2[q,2], shape=lhs.cube.update2[q,3], log=TRUE, type=c("GEV")))
    setTxtProgressBar(pb, q)
  }
  close(pb)
  
  #--------------------LOCATION v SCALE - 1000000 LHS --------------------------------
  lhs.many.values <- cbind(lhs.cube.update2[,1], lhs.cube.update2[,2], lhs.all.params2)
  #labelling the matrix data 
  colnames(lhs.many.values) <- c('location', 'scale', 'likelihood')
  
  #creating empty raster <- still not really sure what this does 
  empty <- extent(lhs.many.values[,1:2])
  
  next.step5 <- raster(empty, ncol=100, nrow=100)
  last.step5 <- rasterize(lhs.many.values[,1:2], next.step5, lhs.many.values[,3], fun=mean)
  plot(last.step5, xlab='location', ylab='scale')
  title("location v. scale, z axis = likelihood")
  
  #log liklelihood
  plot(lhs.cube.update2[seq(from =1, to =1e6, by =1000),1], lhs.all.params2[seq(from =1, to =1e6, by =1000)], ylim=c(-1000,0))
  
  plot(lhs.cube.update2[seq(from =1, to =1e6, by =1000),1], exp(lhs.all.params2[seq(from =1, to =1e6, by =1000)]))
  
  #------------------- SCALE v. SHAPE - 1000000 LHS --------------------------------
  
  
  lhs.many.values2 <- cbind(lhs.cube.update2[,2], lhs.cube.update2[,3], lhs.all.params2)
  colnames(lhs.many.values2) <- c('scale', 'shape', 'likelihood')
  empty2 <- extent(lhs.many.values2[,1:2])
  
  next.step6 <- raster(empty2, ncol=1000, nrow=10)
  last.step6 <- rasterize(lhs.many.values2[,1:2], next.step6, lhs.many.values2[,3], fun=mean)
  plot(last.step6, xlab='scale', ylab='shape', ylim=c(0,50), xlim=c(1000,2000) )
  title("scale v. shape, z axis = likelihood")
  
  #-----------------LOCATION V. SHAPE - 1000000 LHS--------------------------------
  lhs.many.values3 <- cbind(lhs.cube.update2[,1], lhs.cube.update2[,3], lhs.all.params2)
  colnames(lhs.many.values3) <- c('scale', 'shape', 'likelihood')
  empty3 <- extent(lhs.many.values3[,1:2])
  
  next.step7 <- raster(empty3, ncol=1000, nrow=10)
  last.step7 <- rasterize(lhs.many.values3[,1:2], next.step7, lhs.many.values3[,3], fun=mean)
  plot(last.step7, xlab='location', ylab='shape', ylim=c(0,10), xlim=c(1400,1450))
  title("location v. shape, z axis = likelihood")
  
  #order the likelihood in order to find the most likely values of certain parameters 
  lhs.sorted <- rev(order(lhs.all.params2))
  
  #taking top 10 percent of all of this, gives us all of the top , these are just the indicies though
  lhs.top1p <- round(length(lhs.sorted)*.01)
  
  #this is list of all of the indicies that we need 
  ind.to.use <- lhs.sorted[1:lhs.top1p]
  
  #plotting location and scale of top 1 percent - no colors yet 
  plot(lhs.cube.update2[ind.to.use,1], lhs.cube.update2[ind.to.use,2], xlab = 'location', ylab = 'scale')
  title('location v. scale')
  plot(lhs.cube.update2[ind.to.use,2], lhs.cube.update2[ind.to.use,3], xlab = 'scale', ylab ='shape')
  title('scale v. shape')
  plot(lhs.cube.update2[ind.to.use,1], lhs.cube.update2[ind.to.use,3], xlab = 'location', ylab = 'scale')
  title('location v. shape')
  
  #-------trying to create survival function plot with Latin Hyper Cubes ----------------------------------------------------------------
  
  #first lets find the best fitting parameters of the entire data set - aka the highest likelihood
  
  ind.high.like <- lhs.sorted[1]
  
  colnames(lhs.cube.update2) <- c('location', 'scale', 'shape')
  
  sf.loc1 <- lhs.cube.update2[ind.high.like,1]
  sf.scale1 <- lhs.cube.update2[ind.high.like,2]
  sf.shape1 <- lhs.cube.update2[ind.high.like,3]
  
  x.hgt <- seq (0,3000, by=10)
  #list of indicies of max sea level sorted from lowest to highest 
  lsl.sort <- order(lsl.max)
  
  lsl.sorted.vals <- lsl.max[lsl.sort]
  
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
  
  #---------OPTIMIZATION---------------------------------------------------------------------------------------- 
  library(DEoptim)
  
  log.like.calc <- function(data, location1, scale1, shape1){
    ll<-sum(devd(data, loc=location1, scale=scale1, shape=shape1, type=c('GEV'), log=TRUE))
    return(ll)
  }
  
  neg.log.like.calc <- function(p, data){
    mu <- p[1]
    sigma <- p[2]
    xi <- p[3]
    
    nll <- -1*sum(devd(data, loc=mu, scale=sigma, shape=xi, type=c('GEV'), log=TRUE))
    
    return(nll)
  }
  
  
  upper.bound <- c(3000,3000,5)
  lower.bound <- c(0,0,-5)
  
  p.names<- c('mu', 'sigma', 'xi')
  
  de.optim.val <- DEoptim.control(VTR = -Inf, strategy = 2, bs = FALSE, NP = 200,
                  itermax = 200, CR = 0.5, F = 0.8, trace = TRUE, initialpop = NULL,
                  storepopfrom = 101, storepopfreq = 1, p = 0.2, c = 0,
               parallelType = 0, cluster = NULL, packages = c(), parVar = c(),
                  foreachArgs = list())
  
  optim.like <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data= lsl.max)
  
  value.to.use <- pevd(x.hgt, loc = optim.like$optim$bestmem[1], scale=optim.like$optim$bestmem[2], shape=optim.like$optim$bestmem[3], type=c('GEV'), lower.tail = FALSE)
  
  #plotting sf's (both LHS and DEoptim)
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT')
  
  points(lsl.sorted.vals, log10(esf.vals)) 
  lines(x.hgt, log10(sf.hgt2), col='red')
  
  for (s in 1:10){
    sf.hgt.test <- 1-pevd(x.hgt, loc=sf.location.top10[s], scale=sf.scale.top10[s], shape=sf.shapes.top10[s], type=c("GEV"), lower.tail=TRUE)
    lines(x.hgt, log10(sf.hgt.test), col='blue')
  }
  
  lines(x.hgt, log10(value.to.use), col='green')
  
  #-------------non-stationary work 
  
  #read temp data 
  temp.data <- read.table('noaa-temp-1880-2017.csv', header = TRUE, sep=',')
  temp.years <- temp.data$Year
  temp.values <- temp.data$Value[59:134]     #limit temp data to just the years I am working with 
  
  #log likelihood
  log.like.calc <- function(data, location1, scale1, shape1){
    ll<-sum(devd(data, loc=location1, scale=scale1, shape=shape1, type=c('GEV'), log=TRUE))
    return(ll)
  }
  
  #negative log likelihood calculator 
  neg.log.like.calc <- function(p, parnames, data, temps){
   n.parnames <- length(parnames)
   #checks length of parameters, which will be used to determine which parameters are non stationary 
   
   #if all parameters are stationary 
   if (n.parnames == 3){
      mu <- p[1]
      sigma <- p[2]
      xi <- p[3]
   }
   
   #one stationary parameter
   else if(n.parnames ==4){
      if (parnames[1] == 'mu0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi <- p[4]
  
      mu <- mu0 + mu1*temps
      sigma <- exp(sigma)
      }
     else if (parnames[2] == 'sigma0'){
       mu <- p[1]
       sigma0 <- p[2]
       sigma1 <- p[3]
       xi <- p[4]

       sigma <- exp(sigma0+ sigma1*temps)
     }
     else{
       mu <- p[1]
       sigma <- p[2]
       xi0 <- p[3]
       xi1 <- p[4]
       
       sigma <- exp(sigma)
       xi <- xi0 + xi1*temps
     }
   }
   
   #two stationary parameters
     else if (n.parnames == 5){
       if (parnames[1] == 'mu0' & parnames[3]=='sigma0'){
        mu0 <- p[1]
        mu1 <- p[2]
        sigma0 <- p[3]
        sigma1 <- p[4]
        xi <- p[5]
        
        mu <- mu0 + mu1*temps
        sigma <- exp(sigma0 + sigma1*temps)
       }
       else if(parnames[2] == 'sigma0' & parnames[4] == 'xi0'){
         mu <- p[1]
         sigma0 <- p[2]
         sigma1 <- p[3]
         xi0 <- p[4]
         xi1 <- p[5]
         
         sigma <- exp(sigma0 + sigma1*temps)
         xi <- xi0 + xi1*temps
         
       }
       else{
         mu0 <- p[1]
         mu1 <- p[2]
         sigma <- p[3]
         xi0 <- p[4]
         xi1 <- p[5]
         
         mu <- mu0 + mu1*temps
         sigma <- exp(sigma)
         xi <- xi0 + xi1*temps
       }
     }
   #all parameters non stationary 
    else if (n.parnames ==6){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma0 <- p[3]
      sigma1 <- p[4]
      xi0 <- p[5]
      xi1 <- p[6]
      
      mu <- mu0 + mu1*temps
      sigma <- exp(sigma0 + sigma1*temps)
      xi <- xi0 + xi1*temps
      
    }
    nll <- -1*sum(devd(data, loc=mu, scale=sigma, shape=xi, type=c('GEV'), log=TRUE))
  
    return(nll)
  }
  
  #---------Non Stationary Plots Using DE Optim-------------------------------------------------------------------
  
  #----non stationary mu only---------------------------------works!!!!!----------------------
  p.names <- c('mu0', 'mu1', 'sigma', 'xi')
  
  upper.bound <- c(3000,1000, 1000, 5)
  lower.bound <- c(0,0, -100, -5)
  
  optim.like.temp.mu <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)
  
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT, Non Stationary Mu')
  points(lsl.sorted.vals, log10(esf.vals))
  
  fit.vals.optim.mu <- rep(0, length(lsl.max))
  
  for (j in 1:length(temp.values)){
    sf.optim.mu <- 1-pevd(lsl.max[j], loc = (optim.like.temp.mu$optim$bestmem[1] + optim.like.temp.mu$optim$bestmem[2] * temp.values[j]), scale = exp(optim.like.temp.mu$optim$bestmem[3]) , shape = optim.like.temp.mu$optim$bestmem[4],type=c('GEV'))
    fit.vals.optim.mu[j] <- sf.optim.mu
    points(lsl.max[j], log10(sf.optim.mu), pch = 2, col='red')
  }
  #-----non sationary mu and sigma ----------- works ---------------------------------------- --------------------------------------------------- --------------------------------------------------- 
  
  p.names <- c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi')
  
  upper.bound <- c(3000,100, 1000, 10 , 5)
  lower.bound <- c(0,-100, 0, 0, -5)
  
  optim.like.temp.mu.sigma <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values,parnames = p.names)
  
  fit.vals.optim.mu.sigma <- rep(0, length(lsl.max))
  
  #-------run for plot of non sationary mu and sigma 
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT, non stationary mu and sigma')
  points(lsl.sorted.vals, log10(esf.vals))
  for (j in 1:length(temp.values)){
      sf.optim.mu.sigma <- 1 - pevd(lsl.max[j], loc = (optim.like.temp.mu.sigma$optim$bestmem[1] + optim.like.temp.mu.sigma$optim$bestmem[2] * temp.values[j]), scale = exp(optim.like.temp.mu.sigma$optim$bestmem[3] + optim.like.temp.mu.sigma$optim$bestmem[4]*temp.values[j]), shape = optim.like.temp.mu.sigma$optim$bestmem[5],type=c('GEV'))
      fit.vals.optim.mu.sigma[j] <- sf.optim.mu.sigma
      points(lsl.max[j], log10(sf.optim.mu.sigma), pch = 2, col='green')
    
  }
  #------all variables non sationary --------------------------------------this one does work ------------- --------------------------------------------------- --------------------------------------------------- 
  p.names <- c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi0', 'xi1')
  
  upper.bound <- c(3000,100, 1000, 10 , 1, 1)
  lower.bound <- c(0,-100, 0, 0, -1, -1)
  
  optim.like.temp.mu.sigma.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)
  
  fit.vals.optim.mu.sigma.xi <- rep(0, length(lsl.max))
  
  #---------run for plot of non sationary mu and sigma and xi
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT, mu, sigma, xi all nonstationary')
  points(lsl.sorted.vals, log10(esf.vals))
  for (j in 1:length(temp.values)){
      sf.optim.mu.sigma.xi <-1-pevd(lsl.max[j], loc = optim.like.temp.mu.sigma.xi$optim$bestmem[1] + optim.like.temp.mu.sigma.xi$optim$bestmem[2] * temp.values[j], scale = exp(optim.like.temp.mu.sigma.xi$optim$bestmem[3] + optim.like.temp.mu.sigma.xi$optim$bestmem[4]*temp.values[j]), shape = optim.like.temp.mu.sigma.xi$optim$bestmem[5] + optim.like.temp3$optim$bestmem[6]*temp.values[j],type=c('GEV'))
      fit.vals.optim.mu.sigma.xi[j] <- sf.optim.mu.sigma.xi
      points(lsl.max[j], log10(sf.optim.mu.sigma.xi), pch = 2, col='red')
    }
  
  
  
  #------sigma and xi non stationary----------------------works now! ----------------------- --------------------------------------------------- --------------------------------------------------- 
  p.names <- c('mu', 'sigma0', 'sigma1' , 'xi0', 'xi1')
  
  upper.bound <- c(3000, 1000, 10 , 1, 1)
  lower.bound <- c(0, -100, 0, -1, -1)
  
  optim.like.temp.sigma.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)
  
  fit.vals.optim.sigma.xi <- rep(0, length(lsl.max))
  
  
  #---------run for plot of non sationary sigma and xi
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT, non stationary sigma and xi ')
  points(lsl.sorted.vals, log10(esf.vals))
  for (j in 1:length(temp.values)){
    sf.optim5.sigma.xi <-1-pevd(lsl.max[j], loc = optim.like.temp.sigma.xi$optim$bestmem[1], scale = exp(optim.like.temp.sigma.xi$optim$bestmem[2] + optim.like.temp.sigma.xi$optim$bestmem[3]*temp.values[j]), shape = optim.like.temp.sigma.xi$optim$bestmem[4] + optim.like.temp.sigma.xi$optim$bestmem[5]*temp.values[j],type=c('GEV'))
    fit.vals.optim.sigma.xi[j] <- sf.optim5.sigma.xi
    points(lsl.max[j], log10(sf.optim5.sigma.xi), pch = 2, col='red')
  }
  
  
  
  #------non stationary mu and xi----------this one works----------------------------------------- --------------------------------------------------- --------------------------------------------------- 
  p.names <- c('mu0', 'mu1', 'sigma', 'xi0', 'xi1')
  
  upper.bound <- c(3000, 100, 1000, 10, 10)
  lower.bound <- c(0, -100, 0, -1, -10)
  
  optim.like.temp.mu.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)
  
  fit.vals.optim.mu.xi <- rep(0, length(lsl.max))
  #---------run for plot of non sationary mu and xi
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT, nonstationary mu and xi')
  points(lsl.sorted.vals, log10(esf.vals))
  for (j in 1:length(temp.values)){
    sf.optim.mu.xi<-1-pevd(lsl.max[j], loc = optim.like.temp.mu.xi5$optim$bestmem[1] +optim.like.temp.mu.xi$optim$bestmem[2]*temp.values[j], scale = exp(optim.like.temp.mu.xi$optim$bestmem[3]), shape = optim.like.temp.mu.xi$optim$bestmem[4] + optim.like.temp.mu.xi$optim$bestmem[5]*temp.values[j],type=c('GEV'))
    fit.vals.optim.mu.xi[j] <- sf.optim.mu.xi
    points(lsl.max[j], log10(sf.optim.mu.xi), pch = 2, col='red')
  }
  
  #------only sigma non stationary -------------------------this works!!-------------------------- --------------------------------------------------- --------------------------------------------------- 
  
  p.names <- c('mu', 'sigma0', 'sigma1', 'xi')
  
  upper.bound <- c(3000, 1000, 10, 1)
  lower.bound <- c(0, 0, -10, -1)
  
  optim.like.temp.sigma <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)
  
  fit.vals.optim.sigma <- rep(0, length(lsl.max))
  #---------run for plot of non sationary sigma 
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT, sigma non stationary')
  points(lsl.sorted.vals, log10(esf.vals))
  for (j in 1:length(temp.values)){
    sf.optim.sigma <-1-pevd(lsl.max[j], loc = optim.like.temp.sigma$optim$bestmem[1], scale = exp(optim.like.temp.sigma$optim$bestmem[2] + optim.like.temp.sigma$optim$bestmem[3]*temp.values[j]), shape = optim.like.temp.sigma$optim$bestmem[4],type=c('GEV'))
    fit.vals.optim.sigma[j] <- sf.optim.sigma
    points(lsl.max[j], log10(sf.optim.sigma), pch = 2, col='red')
  }
  
  #------ only xi non stationary---------------------------------this works! ------------------ --------------------------------------------------- ---------------------------------------------------  
  p.names <- c('mu', 'sigma', 'xi0', 'xi1')
  
  upper.bound <- c(3000, 500 , 1, 1)
  lower.bound <- c(0, 0, -1, -1)
  
  optim.like.temp.xi <- DEoptim(neg.log.like.calc, lower=lower.bound, upper=upper.bound, control=de.optim.val, data = lsl.max, temps = temp.values, parnames = p.names)
  
  fit.vals.optim.xi <- rep(0, length(lsl.max))
  #---------run for plot of non sationary xi
  plot(x.hgt, log10(sf.hgt), type='l', ylab = 'log probability', xlab = 'height [mm] of sea')
  title('Survival Function of height of sea [mm], New London CT')
  points(lsl.sorted.vals, log10(esf.vals))
  for (j in 1:length(temp.values)){
    sf.optim.xi <-1-pevd(lsl.max[j], loc = optim.like.temp.xi$optim$bestmem[1], scale = exp(optim.like.temp.xi$optim$bestmem[2]), shape = optim.like.temp.xi$optim$bestmem[3] + optim.like.temp.xi$optim$bestmem[4]*temp.values[j],type=c('GEV'))
    fit.vals.optim.xi[j] <- sf.optim.xi
    points(lsl.max[j], log10(sf.optim.xi), pch = 2, col='red')
  }
  
  
  
  #--------------------------RMSE For Optim Fits------------------------------------------------------------------------------------------------
  
  #only mu non stationary 
  rmse.mu <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.mu)^2))
  
  #only sigma non stationary 
  rmse.sigma<- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.sigma)^2))
  
  #xi non stationary 
  rmse.xi <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.xi))^2)
  
  #mu and sigma non stationary 
  rmse.mu.sigma <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.mu.sigma)^2))
  
  #mu and xi non stationary 
  rmse.mu.xi<- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.mu.xi))^2)
  
  #sigma and xi non stationary 
  rmse.sigma.xi <- sqrt(1/(length(lsl.max))*(sum(esf.vals - fit.vals.optim.sigma.xi))^2)
  
  #mu, sigma, xi non stationary 
  rmse.mu.sigma.xi <- sqrt(1/(length(lsl.max))*(sum(esf.vals - sf.optim.mu.sigma.xi))^2)
  
  #----------------AIC For Optim Fits------------------------------------------------------------------------------------------------------------------------
  #2k - 2 ln (L)
  
  aic.mu <- 2*4 - 2*log(optim.like.temp.mu$optim$bestval)
  
  aic.sigma <- 2*4 - 2*log(optim.like.temp.sigma$optim$bestval)
  
  aic.xi <- 2*4 - 2*log(optim.like.temp.xi$optim$bestval)
  
  aic.mu.sigma <- 2*5 - 2*log(optim.like.temp.mu.sigma$optim$bestval)
  
  aic.mu.xi<- 2*5 - 2*log(optim.like.temp.mu.xi$optim$bestval)
  
  aic.sigma.xi<- 2*5 - 2*log(optim.like.temp.sigma.xi$optim$bestval)
  
  aic.mu.sigma.xi <- 2*6 - 2*log(optim.like.temp.mu.sigma.xi$optim$bestval)
  
  #------------------BIC For Optim Fits------------------------------------------------------------------------------------------------------------------------ 
  #kln(n) - 2 ln(L)
  
  bic.mu<- 4*log(length(lsl.max)) - 2*log(optim.like.temp.mu$optim$bestval)
  
  bic.sigma <- 4*log(length(lsl.max))- 2*log(optim.like.temp.sigma$optim$bestval)
  
  bic.xi <- 4*log(length(lsl.max))- 2*log(optim.like.temp.xi$optim$bestval)
  
  bic.mu.sigma<- 5*log(length(lsl.max))- 2*log(optim.like.temp.mu.sigma$optim$bestval)
  
  bic.mu.xi <- 5*log(length(lsl.max))- 2*log(optim.like.temp.mu.xi$optim$bestval)
  
  bic.sigma.xi <- 5*log(length(lsl.max))- 2*log(optim.like.temp.sigma.xi$optim$bestval)
  
  bic.mu.sigma.xi <- 6*log(length(lsl.max))- 2*log(optim.like.temp.mu.sigma.xi$optim$bestval)
  
  
  #
  #==============================================================================
  ## End
  ##==============================================================================
