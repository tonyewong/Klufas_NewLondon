
#make a fake GEV 
library(mvtnorm)
library(extRemes)
library(coda)

#set params - true values 
mu <- 400

sigma  <- 150

xi <- .02 

x.gev <- seq(from=0,to=1000,by=5)

#the true value curve
stationary.curve <- devd(x.gev, loc = mu, scale = sigma, shape = xi, log=TRUE, type=c('GEV'))
plot(stationary.curve, type='l')

#time to make some data

epsilon.gev <- revd(length(x.gev), loc=mu, scale = sigma, shape = xi, type=c("GEV"))

epsilon.gev.devd <- devd(epsilon.gev, loc = mu, scale = sigma, shape = xi, log=TRUE, type=c('GEV'))

x.gev.star <-  epsilon.gev.devd 

plot(x.gev.star)
lines(stationary.curve)

plot(epsilon.gev, epsilon.gev.devd)
lines(x.gev, stationary.curve)

log.like2 <- function(p, data, loc, sigma.guess, xi){ # p = theta.1, data = y.func 
  val <- sum(devd(data, loc = p[1], scale = p[2], shape = p[3], log=TRUE, type=c("GEV")))#make mu = 0 here 
  return(val)
}

log.like.pri2 <- function(p, loc, sd, min, max){
  #should be for loop for all of the parameters in the function, will just start w two - alpha and beta 
  mu <- p[1]
  sigma <- p[2]
  xi <- p[3]
  
  #let alpha be normally distrtibuted 
  p.mu <- sum(dunif(x=mu, min= 0, max=3000, log = TRUE))
  #let beta be uniformly distributed 
  p.sigma <- sum(dunif(x=sigma, min= 0, max=3000, log = TRUE))
  
  p.xi <- sum(dunif(x=xi,min= -10, max=10, log = TRUE))
  
  #add together b/c logs
  p.mu.sigma.xi <- p.mu  + p.sigma + p.xi
  return(p.mu.sigma.xi)
}



                                   


log.post.final2 <- function(p, data, loc, sigma.guess, sd, xi, min, max){
  
  log.like.pri <- log.like.pri2(p, loc, sd, min, max)
  
  #log.like.pri <- log.like.pri2(p, loc, sd, xi, min, max)
  
  if (is.finite(log.like.pri)){
    log.like <- log.like2(p, data, loc, sigma.guess, xi)
    
    log.like.final <- (log.like + log.like.pri)
  }
  else{
    log.like.final <-log.like.pri  
  }
  #log.like.final <- (log.like + log.like.pri)
  
  return(log.like.final)
}

#list prior / fixed params here 
priors <- vector('list', 5)

names(priors) <- c('mu', 'sigma', 'xi', 'min', 'max')

priors$mu <- 400

priors$sigma <- 150 

priors$xi <- .02

priors$min <- 0 

priors$max <- 600

mcmc.function.gev <- function(niter, initial.mu, initial.sigma, initial.xi, data, initial.mu.step , initial.sigma.step, initial.xi.step){
  #theta.vals <- cbind(rep(NA, niter), rep(NA, niter), rep(NA, niter)) 
  
  #creating matrix that is niter by 3 large 
  theta.vals <- matrix(nrow=niter, ncol=3)
  #theta.vals <- mat.or.vec( niter, 3)
  
  #making array to store the probabilities of theta params 
  theta.prob.vals <- rep(0, niter)
  
  #start with theta.1
  theta.1 <- c(initial.mu, initial.sigma, initial.xi)
  
  #calculate probability of this given theta set I chose 
  theta.1.prob <- log.post.final2(p=theta.1, 
                                  data=data, 
                                  loc=priors$mu,
                                  sigma.guess = priors$sigma , 
                                  sd =priors$sigma , 
                                  xi = priors$xi, 
                                  min=priors$min, 
                                  max=priors$max) #use y.func as parameter when using this function 
  
  #place theta values and probability values into vectors to old
  theta.vals[1,] <- theta.1
  theta.prob.vals[1] <- theta.1.prob
  
  param.changes <- rep(0,niter)
  #set stepsizes for alpha and beta 
  
  step.size.mu <- initial.mu.step
  
  step.size.sigma <-  initial.sigma.step
  
  step.size.xi <- initial.xi.step 
  
  #covar.vals <- list(rep(0,niter))
  
  sigma.sd <- ((2.4)^2 /3) *matrix(c(step.size.mu, 0, 0, 0, step.size.sigma , 0, 0, 0, step.size.beta), ncol=3)
  
  #covar.vals[1] <- sigma.sd
  
  #for loop 
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for (i in 2:niter){
    #going to try and use mvn to calculate alpha and beta
    theta.prev <- theta.vals[i-1,]
    
    if (i > 3e3){ #cov of the alpha and beta vectors 
     # print("if statement")
      sigma.sd2 <- ((2.4)^2 /3) * cov(theta.vals[1:(i-1),]) #changed from theta.vals[1:i-1,]
      
      #cov <- matrix(c(step.size.alpha, 0, 0, step.size.beta), ncol=2)
      #covar.vals[i] <- sigma.sd2
      
      mu <- c(theta.prev[1],theta.prev[2], theta.prev[3])
      
      #mu.sigma.xi.dependent <- revd(n=1, loc = theta.prev[1], scale = theta.prev[2], shape =theta.prev[3] , type=c("GEV"))
      
      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = mu, sigma = sigma.sd2) #sigma = covariance matrix
      
      mu.star <- mu.sigma.xi.dependent[1]
      sigma.star <- mu.sigma.xi.dependent[2]
      xi.star <- mu.sigma.xi.dependent[3]
      
    }
    else{
    
     # covar.vals[i] <- sigma.sd
      mu <- c(theta.prev[1],theta.prev[2], theta.prev[3])
      
      #mu.sigma.xi.dependent <- revd(n=1, loc = theta.prev[1], scale = theta.prev[2], shape =theta.prev[3] , type=c("GEV"))
      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = mu, sigma = sigma.sd) 
      #print(mu.sigma.xi.dependent)
      mu.star <- mu.sigma.xi.dependent[1]
      sigma.star <- mu.sigma.xi.dependent[2]
      xi.star <- mu.sigma.xi.dependent[3]
      
      
    }
    
    current.theta <- c(mu.star, sigma.star, xi.star)
    
    #print(current.theta)
    
    
    #calculate probability of new params
    theta.star.prob <- log.post.final2(p=current.theta, 
                                       data = data, 
                                       loc = priors$mu, 
                                       sigma.guess = priors$sigma , 
                                       sd = priors$sigma , 
                                       xi = priors$xi, 
                                       min = priors$min, 
                                       max = priors$max )
    
   # print(theta.star.prob)
    
    
   # print(theta.prob.prev)
    
    #pull probability value of prev params 
    theta.prob.prev <- theta.prob.vals[i-1]
    
    #calculate alpha
    alpha.prob <- min(0, theta.star.prob - theta.prob.prev) 
    
    #decide which parameter to continue with 
    u <- runif(n=1, min=0, max=1)
    
    #if less than alpha.prob, new theta value is assigned 
    if (log(u) < alpha.prob) {
      current.theta.prob <- theta.star.prob
      param.changes[i] <- 1
    }else{
      #if greater than or equal to alpha.prob, prev theta value assigned, theta params from before also replace current theta values 
      
      current.theta.prob <- theta.prob.prev
      current.theta <- theta.vals[i-1, ]
    }
   
    #theta vals are the values for alpha and beta that were used 
    theta.vals[i,] <- current.theta
    theta.prob.vals[i] <- current.theta.prob
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #need to return list of a lot of values from this function 
  mcmc.out <- vector('list', 3) #chains, acceptances, lpost, covar
  names(mcmc.out) <- c('chains', 'acceptances', 'lpost')
  
  mcmc.out$chains <- theta.vals
  mcmc.out$acceptances <- param.changes
  mcmc.out$lpost <- theta.prob.vals
  return(mcmc.out)
}

gev.run1 <- mcmc.function.gev(niter = 1e5, initial.mu = 50, initial.sigma = 50, initial.xi = 1, data =epsilon.gev, initial.mu.step = .2 , initial.sigma.step = 1, initial.xi.step = .05)                        
gev.run2 <- mcmc.function.gev(niter = 1e5, initial.mu = 100, initial.sigma = 100, initial.xi = .2, data =epsilon.gev, initial.mu.step = .5 , initial.sigma.step = .05, initial.xi.step = .1)                        
gev.run3 <- mcmc.function.gev(niter = 1e5, initial.mu = 200, initial.sigma = 75, initial.xi = .3, data =epsilon.gev, initial.mu.step = .02 , initial.sigma.step = .5, initial.xi.step = .5)                        


mcmc1 <- as.mcmc(gev.run1$chains)
mcmc2 <- as.mcmc(gev.run2$chains)
mcmc3 <- as.mcmc(gev.run3$chains)

#compare the different values between the chains of alpha and beta 
mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
plot(mcmc.chain.list)

gelman.plot(mcmc.chain.list)


#--------------i_burnin-----------------------------
niter <- 1e5
my.seq <- seq(from = 5e3, to = niter, by = 1e3)
list.of.gelman.diags <- rep(0, length(my.seq))

for (i in 1:length(my.seq)){
  mcmc1 <- as.mcmc(gev.run1$chains[1:my.seq[i],])
  mcmc2 <- as.mcmc(gev.run2$chains[1:my.seq[i],])
  mcmc3 <- as.mcmc(gev.run3$chains[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags[i] <- gelman.diag(mcmc.chain.list)[2]
}

#-----------histograms and plots of location----------------------
plot(gev.run1$chains[,1], type= 'l')
plot(gev.run2$chains[,1], type = 'l')
plot(gev.run3$chains[,1], type = 'l')


#location histograms 
hist(gev.run1$chains[1000:1e5,1])
abline(v=mean(gev.run1$chains[1000:1e5,1]), col='red')

hist(gev.run2$chains[1000:1e5,1])
abline(v=mean(gev.run2$chains[1000:1e5,1]), col='red')

hist(gev.run3$chains[1000:1e5,1])
abline(v=mean(gev.run3$chains[1000:1e5,1]), col='red') #5% and 95% quantile lines also , also true value 

hist(gev.run1$chains[6000:1e5,1])
hist(gev.run2$chains[6000:1e5,1])
hist(gev.run3$chains[6000:1e5,1])

#-----------histograms and plots of scale----------------------

plot(gev.run1$chains[,2], type= 'l')
plot(gev.run2$chains[,2], type = 'l')
plot(gev.run3$chains[,2], type = 'l')


#location histograms 
hist(gev.run1$chains[1000:1e5,2])
abline(v=mean(gev.run1$chains[1000:1e5,2]), col='red')

hist(gev.run2$chains[1000:1e5,2])
abline(v=mean(gev.run2$chains[1000:1e5,2]), col='red')

hist(gev.run3$chains[1000:1e5,2])
abline(v=mean(gev.run3$chains[1000:1e5,2]), col='red') #5% and 95% quantile lines also , also true value 

hist(gev.run1$chains[6000:1e5,2])
hist(gev.run2$chains[6000:1e5,2])
hist(gev.run3$chains[6000:1e5,2])

#-----------histograms and plots of shape----------------------

plot(gev.run1$chains[,3], type= 'l')
plot(gev.run2$chains[,3], type = 'l')
plot(gev.run3$chains[,3], type = 'l')


#location histograms 
hist(gev.run1$chains[1000:1e5,3])
abline(v=mean(gev.run1$chains[1000:1e5,3]), col='red')

hist(gev.run2$chains[1000:1e5,3])
abline(v=mean(gev.run2$chains[1000:1e5,3]), col='red')

hist(gev.run3$chains[1000:1e5,3])
abline(v=mean(gev.run3$chains[1000:1e5,3]), col='red') #5% and 95% quantile lines also , also true value 

hist(gev.run1$chains[6000:1e5,3])
hist(gev.run2$chains[6000:1e5,3])
hist(gev.run3$chains[6000:1e5,3])

#------------------trying to incorperate new data wow!------------------------------------------------------
setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
vals <- read.tide.data()


gev.run.sl.1 <- mcmc.function.gev(niter = 1e5, 
                              initial.mu = 50, 
                              initial.sigma = 50, 
                              initial.xi = 1, 
                              data = vals$max, 
                              initial.mu.step = .2 , 
                              initial.sigma.step = 1, 
                              initial.xi.step = .05)   


gev.run2 <- mcmc.function.gev(niter = 1e5, 
                              initial.mu = 100, 
                              initial.sigma = 100, 
                              initial.xi = .2, 
                              data =vals$max, 
                              initial.mu.step = .5 , 
                              initial.sigma.step = .05, 
                              initial.xi.step = .1) 


gev.run3 <- mcmc.function.gev(niter = 1e5, 
                              initial.mu = 200, 
                              initial.sigma = 75, 
                              initial.xi = .3, 
                              data =vals$max, 
                              initial.mu.step = .02 ,
                              initial.sigma.step = .5, 
                              initial.xi.step = .5)   


mcmc1 <- as.mcmc(gev.run.sl.1$chains)
mcmc2 <- as.mcmc(gev.run2$chains)
mcmc3 <- as.mcmc(gev.run3$chains)

#compare the different values between the chains of alpha and beta 
mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
plot(mcmc.chain.list)

gelman.plot(mcmc.chain.list)

#--------------i_burnin-----------------------------
niter <- 1e5
my.seq <- seq(from = 5e3, to = niter, by = 1e3)
list.of.gelman.diags <- rep(0, length(my.seq))

for (i in 1:length(my.seq)){
  mcmc1 <- as.mcmc(gev.run.sl.1$chains[1:my.seq[i],])
  mcmc2 <- as.mcmc(gev.run2$chains[1:my.seq[i],])
  mcmc3 <- as.mcmc(gev.run3$chains[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags[i] <- gelman.diag(mcmc.chain.list)[2]
}

#20th elt of this list is where it converges --> 23000 place 

#--------------location ------------------------------------------

plot(gev.run.sl.1$chains[(2.3e4:1e5),1], type='l')
plot(gev.run2$chains[(2.3e4:1e5),1], type='l')
plot(gev.run3$chains[(2.3e4:1e5),1], type='l')


#location histograms 
hist(gev.run.sl.1$chains[(2.3e4:1e5),1])
percentiles <- quantile(gev.run.sl.1$chains[(2.3e4:1e5),1], c(.05,.95))
abline(v=mean(gev.run.sl.1$chains[(2.3e4:1e5),1]), col='red')
abline(v=percentiles, col='blue')

hist(gev.run2$chains[(2.3e4:1e5),1])
percentiles <- quantile(gev.run2$chains[(2.3e4:1e5),1], c(.05,.95))
abline(v=mean(gev.run2$chains[(2.3e4:1e5),1]), col='red')
abline(v=percentiles, col='blue')

hist(gev.run3$chains[(2.3e4:1e5),1])
percentiles <- quantile(gev.run.3$chains[(2.3e4:1e5),1], c(.05,.95))
abline(v=mean(gev.run3$chains[(2.3e4:1e5),1]), col='red') #5% and 95% quantile lines also , also true value 
abline(v=percentiles, col='blue')


#-----------scale-------------------------------------------------


plot(gev.run.sl.1$chains[(2.3e4:1e5),2], type='l')
plot(gev.run2$chains[(2.3e4:1e5),2], type='l')
plot(gev.run3$chains[(2.3e4:1e5),2], type='l')


#location histograms 
hist(gev.run.sl.1$chains[(2.3e4:1e5),2])
percentiles <- quantile(gev.run.sl.1$chains[(2.3e4:1e5),2], c(.05,.95))
abline(v=mean(gev.run.sl.1$chains[(2.3e4:1e5),2]), col='red')
abline(v=percentiles, col='blue')

hist(gev.run2$chains[(2.3e4:1e5),2])
percentiles <- quantile(gev.run2$chains[(2.3e4:1e5),2], c(.05,.95))
abline(v=mean(gev.run2$chains[(2.3e4:1e5),2]), col='red')
abline(v=percentiles, col='blue')

hist(gev.run3$chains[(2.3e4:1e5),2])
percentiles <- quantile(gev.run.3$chains[(2.3e4:1e5),2], c(.05,.95))
abline(v=mean(gev.run3$chains[(2.3e4:1e5),2]), col='red') #5% and 95% quantile lines also , also true value 
abline(v=percentiles, col='blue')


#-----------shape-------------------------------------------------

plot(gev.run.sl.1$chains[(2.3e4:1e5),3], type='l')
plot(gev.run2$chains[(2.3e4:1e5),3], type='l')
plot(gev.run3$chains[(2.3e4:1e5),3], type='l')


#location histograms 
hist(gev.run.sl.1$chains[(2.3e4:1e5),3])
percentiles <- quantile(gev.run.sl.1$chains[(2.3e4:1e5),3], c(.05,.95))
abline(v=mean(gev.run.sl.1$chains[(2.3e4:1e5),3]), col='red')
abline(v=percentiles, col='blue')

hist(gev.run2$chains[(2.3e4:1e5),3])
percentiles <- quantile(gev.run2$chains[(2.3e4:1e5),3], c(.05,.95))
abline(v=mean(gev.run2$chains[(2.3e4:1e5),3]), col='red')
abline(v=percentiles, col='blue')

hist(gev.run3$chains[(2.3e4:1e5),3])
percentiles <- quantile(gev.run.3$chains[(2.3e4:1e5),3], c(.05,.95))
abline(v=mean(gev.run3$chains[(2.3e4:1e5),3]), col='red') #5% and 95% quantile lines also , also true value 
abline(v=percentiles, col='blue') 

