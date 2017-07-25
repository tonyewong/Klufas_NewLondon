

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data(1939,2014)

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
tide.data <- read.tide.data()


#setwd('~/codes/Klufas_NewLondon/R/')
#source('gev_fitting_MCMC.R')

#doing GEV MCMC fitting for non stationary models

#leggo

log.like.calc.nonstat <- function(p, parnames, data){ #non stationary
  n.parnames <- length(parnames)
  #checks length of parameters, which will be used to determine which parameters are non stationary

  #if all parameters are stationary
  if (n.parnames == 3){
    mu <- p[1]
    sigma <- p[2]
    xi <- p[3]
    #vals <- temps$values*2  #because I was getting error that "temps" wasn't being used for some reason
  }

  #one stationary parameter
  else if(n.parnames == 4){
    if (parnames[1] == 'mu0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi <- p[4]

      mu <- mu0 + mu1*temps$values
     # print(mu)
    }
    else if (parnames[2] == 'sigma0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi <- p[4]

      sigma <- exp(sigma0+ sigma1*temps$values)
    }
    else{
      mu <- p[1]
      sigma <- p[2]
      xi0 <- p[3]
      xi1 <- p[4]

      #sigma <- exp(sigma)
      xi <- xi0 + xi1*temps$values
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

      mu <- mu0 + mu1*temps$values
      sigma <- exp(sigma0 + sigma1*temps$values)
    }
    else if(parnames[2] == 'sigma0' & parnames[4] == 'xi0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]

      sigma <- exp(sigma0 + sigma1*temps$values)
      xi <- xi0 + xi1*temps$values

    }
    else{
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]

      mu <- mu0 + mu1*temps$values
      #sigma <- exp(sigma)
      xi <- xi0 + xi1*temps$values
    }
  }
  #all parameters non stationary
  else if (n.parnames == 6){
    mu0 <- p[1]
    mu1 <- p[2]
    sigma0 <- p[3]
    sigma1 <- p[4]
    xi0 <- p[5]
    xi1 <- p[6]

    mu <- mu0 + mu1*temps$values #mu0 + mu1*temps$values
    sigma <- exp(sigma0 + sigma1*temps$values)
    xi <- xi0 + xi1*temps$values

  }

  nll <- sum(devd(data, loc=mu, scale=sigma, shape=xi, type=c('GEV'), log=TRUE))
  return(nll)
}

log.like.pri.nonstat<- function(p, parnames, min, max, temps){
  #should be for loop for all of the parameters in the function, will just start w two - alpha and beta
  n.parnames <- length(parnames)
  #checks length of parameters, which will be used to determine which parameters are non stationary

  #if all parameters are stationary
  if (n.parnames == 3){
    mu <- p[1]
    sigma <- p[2]
    xi <- p[3]
    
    p.mu <- dunif(x = mu, min = 0, max = 3000, log = TRUE)
    p.sigma <- dunif(x=sigma, min = 0, max = 400, log = TRUE)
    p.xi <- dunif(x=xi, min = -10, max = 10, log = TRUE)
    
    lpri <- p.mu + p.sigma + p.xi
  }

  #one stationary parameter
  else if(n.parnames ==4){
    if (parnames[1] == 'mu0'){
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi <- p[4]

      p.mu0    <-  dunif(x= mu0, min = 0, max = 3000, log = TRUE)
      p.mu1    <-  dunif(x= mu1, min = -500, max = 500, log = TRUE)
      p.sigma  <-  dunif(x = sigma, min = 0, max = 400, log = TRUE)
      p.xi     <-  dunif(x = xi, min = -10, max = 10, log = TRUE) 
      
      lpri <- p.mu0 + p.mu1 + p.sigma + p.xi

# 
#       mu <-mu0 + mu1*temps
#       #print(mu)

    }
    else if (parnames[2] == 'sigma0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi <- p[4]
      
      p.mu <-  dunif( x= mu, min = 0, max = 3000, log = TRUE)
      p.sigma0 <- dunif(x =sigma0, min = -300, max = 300, log = TRUE)
      p.sigma1 <- dunif(x=sigma1, min = -100, max = 200, log = TRUE)
      p.xi <- dunif(x = xi, min = -10, max = 10, log = TRUE) 

      lpri <- p.mu + p.sigma0 + p.sigma1 + p.xi


      sigma <- exp(sigma0+ sigma1*temps)

    }
    else{
      mu <- p[1]
      sigma <- p[2]
      xi0 <- p[3]
      xi1 <- p[4]

      
      p.mu <-  dunif( x= mu, min = 0, max = 3000, log = TRUE)
      p.sigma <- dunif(x = sigma, min = 0, max = 1000, log = TRUE)
      p.xi0 <- dunif(x = xi0, min = -10, max = 10, log = TRUE)
      p.xi1 <- dunif(x=xi1, min = -10, max = 10, log = TRUE) 
      
      lpri <- p.mu + p.sigma + p.xi0 + p.xi1


      # xi <- xi0 + xi1*temps

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
      
      p.mu0 <-  dunif( x= mu0, min = 0, max = 3000, log = TRUE)
      p.mu1 <- dunif(x=mu1, min = -500, max = 1500, log = TRUE)
      p.sigma0 <- dunif(x = sigma0, min = 0, max = 300, log = TRUE)
      p.sigma1 <- dunif(x= sigma1, min = -100, max = 100, log = TRUE)
      p.xi <- dunif(x = xi, min = -10, max = 10, log = TRUE)
      
      lpri <- p.mu0 + p.mu1 + p.sigma0 + sigma1 + p.xi
      

      # print('here')
      # mu <- mu0 + mu1*temps
      # sigma <- exp(sigma0 + sigma1*temps)

    }
    else if(parnames[2] == 'sigma0' & parnames[4] == 'xi0'){
      mu <- p[1]
      sigma0 <- p[2]
      sigma1 <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]

      
      p.mu <-  dunif( x= mu, min = 0, max = 3000, log = TRUE)
      p.sigma0 <- dunif(x = sigma0 , min = -100, max = 100, log = TRUE)
      p.sigma1 <- dunif(x = sigma1, min = -150 , max = 150, log = TRUE)
      p.xi0 <- dunif(x = xi0, min= -10, max = 10, log = TRUE)
      p.xi1 <- dunif( x = xi1, min = -10, max = 10, log = TRUE)
      
      lpri <- p.mu + p.sigma0 +p.sigma1 + p.xi0 +p.xi1
      #print('here')


    }
    else{
      mu0 <- p[1]
      mu1 <- p[2]
      sigma <- p[3]
      xi0 <- p[4]
      xi1 <- p[5]
      #print('here')
      
      p.mu0 <-  dunif( x= mu0 , min = 0, max = 3000, log = TRUE)
      p.mu1 <-  dunif( x= mu1, min = -100, max = 100, log = TRUE)
      p.sigma <- dunif(x = sigma, min = 0, max = 1000, log = TRUE)
      p.xi0 <- dunif(x = xi0 , min = -10, max = 10, log = TRUE)
      p.xi1 <- dunif(x = xi1, min = -10, max = 10, log = TRUE)
      
      lpri <- p.mu0 +p.mu1 + p.sigma + p.xi0 +p.xi1
     # mu <- mu0 + mu1*temps
     # xi <- xi0 + xi1*temps
    }
  }
  #all parameters non stationary
  else if (n.parnames ==6){
    #print('hererrr')
    mu0 <- p[1]
    mu1 <- p[2]
    sigma0 <- p[3]
    sigma1 <- p[4]
    xi0 <- p[5]
    xi1 <- p[6]

    
    p.mu0 <-  dunif( x= mu0, min = 0, max = 3000, log = TRUE)
    p.mu1 <-  dunif( x=mu1, min = -500, max = 1500, log = TRUE)
    p.sigma0 <- dunif(x = sigma0 , min = -500, max = 500, log = TRUE)
    p.sigma1 <- dunif(x = sigma1, min = -500, max = 500, log = TRUE)
    p.xi0 <- dunif(x = xi0, min = -10, max = 10, log = TRUE)
    p.xi1 <- dunif(x = xi1, min = -10, max = 10, log = TRUE)
    
    lpri <- p.mu0 + p.mu1 + p.sigma0 + p.sigma1 + p.xi0 + p.xi1

  }
  


# TW -- ah-ha! So this is a spot where toruble might be brewing. The prior distributions should
# be defined for the actual parameters, mu0, mu1, sigma0, sigma1, etc... and not only for
# mu/sigma/xi. So you'll have this calculation specifically for each of the cases above in
# the "if" statement checking for the different models. And this will help with sigma, because
# the stationary model has a reasonable prior for sigma on [0, 3000] (as you have below) but
# the nonstationary parameters sigma0 and sigma1 will be much smaller. sigma0 I think can be
# distributed uniformly on [0 something] (play around with the upper bound, and see from the
# maximum likelihood estimates what a reasonable value for that is), and sigma1 can be positive or negative, so
# use maybe a uniform prior centered at 0, something like [-something, +something] (again using
# the MLE values you found for sigma1 as a guide for what a reasonable upper/lower bound is).

  #p.mu <- sum(dunif(x=mu, min= 0, max=3000, log = TRUE))

  #print('here')
 # p.sigma <- sum(dunif(x=sigma, min= 0, max=3000, log = TRUE))

  #p.xi <- sum(dunif(x=xi, min= -10, max=10, log = TRUE))


  #add together b/c logs
  #p.mu.sigma.xi <- p.mu  + p.sigma + p.xi
  return(lpri)
}

log.post.final2 <- function(p, parnames, data, min, max, temps){

  log.like.pri.nonstat <- log.like.pri.nonstat(p, parnames, min, max, temps)

  #log.like.pri <- log.like.pri2(p, loc, sd, xi, min, max)

  if (is.finite(log.like.pri.nonstat)){
    log.like <- log.like.calc.nonstat(p, parnames, data)
    #print('here2')
    log.like.final <- (log.like + log.like.pri.nonstat)
  }
  else{

    log.like.final <-log.like.pri.nonstat 
    #print('here3')

    log.like.final <-log.like.pri.nonstat
    print('here3')

  }

  return(log.like.final)
}

matrix.maker <- function(size, step.sizes){#size = dimension x dimension of the matrix, step size = vector of the values down the diagonal
  matrix <- matrix(rep(0, size^2), ncol=size)
  for (i in 1:size){
    matrix[i,i] <- step.sizes[i]
  }
  return(matrix)
}
#list prior / fixed params here
priors <- vector('list', 5)

names(priors) <- c('mu', 'sigma', 'xi', 'min', 'max')

priors$mu <- 400

priors$sigma <- 150

priors$xi <- .02

priors$min <- 0

priors$max <- 3000

mcmc.function.gev <- function(niter, params, parnames, data, initial.steps){
  #theta.vals <- cbind(rep(NA, niter), rep(NA, niter), rep(NA, niter))
  n.parnames <- length(parnames)

  #creating matrix that is niter by a certain number of parameters large
  theta.vals <- matrix(nrow=niter, ncol=n.parnames)

  #theta.vals <- mat.or.vec( niter, 3)

  #making array to store the probabilities of theta params
  theta.prob.vals <- rep(0,niter)

  #start with theta.1
  for (i in 1:n.parnames){
    theta.vals[1,i] <- params[i]
  }
  #theta.1 <- c(initial.mu, initial.sigma, initial.xi)
  theta.1 <- theta.vals[1,]
  #calculate probability of this given theta set I chose
  theta.1.prob <- log.post.final2(p=theta.1,
                                  parnames = parnames,
                                  data = data,
                                  temps = temps$values,
                                  min=priors$min,
                                  max=priors$max) #use y.func as parameter when using this function

  #place theta values and probability values into vectors to old
  theta.prob.vals[1] <- theta.1.prob

  param.changes <- rep(0,niter)

  sigma.sd <- ((2.4)^2 /n.parnames) *matrix.maker(n.parnames , initial.steps)


  #for loop
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for (i in 2:niter){
    #going to try and use mvn to calculate alpha and beta
    theta.prev <- theta.vals[i-1,]

    if (i > 3e3){ #cov of the alpha and beta vectors
      # print("if statement")
      sigma.sd2 <- ((2.4)^2 /n.parnames) * cov(theta.vals[1:(i-1),]) #changed from theta.vals[1:i-1,]

      #cov <- matrix(c(step.size.alpha, 0, 0, step.size.beta), ncol=2)
      #covar.vals[i] <- sigma.sd2


      #mu.sigma.xi.dependent <- revd(n=1, loc = theta.prev[1], scale = theta.prev[2], shape =theta.prev[3] , type=c("GEV"))

      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = theta.prev, sigma = sigma.sd2) #sigma = covariance matrix

    }
    else{

      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = theta.prev, sigma = sigma.sd)


    }

    current.theta <- mu.sigma.xi.dependent

    theta.star.prob <- log.post.final2(p = current.theta,
                                       parnames = parnames,
                                       data = data,
                                       temps = temps$values,
                                       min=priors$min,
                                       max=priors$max)

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

#make some fake data , use real temp data

#---------------------- lets start w non stationary location-------------------------------------------------

mu0 <- 400

mu1 <- 200

sigma  <- 150

xi <- .02

x.hgt <- seq (0,3000, by=10)

x.gev <- seq(from=0,to=75,by=1)

simple.curve <- rep(0, length(temps$values))
#the true value curve
for (i in 1:length(temps$values)){
  simple.curve[i] <- sum(devd(x.gev, loc = mu0  + mu1*temps$values[i], scale = sigma, shape = xi, log=TRUE, type=c('GEV')))
}

hist(simple.curve)

plot(simple.curve, type='l')

#time to make some data ---- this is what I will use

epsilon.gev <- rep(0, length(temps$values))

for(i in 1:length(temps$values)){ #synthetic data
  epsilon.gev[i] <- revd(1, loc=mu0 + mu1*temps$values[i], scale = sigma, shape = xi, type=c("GEV"))
}

plot(epsilon.gev)
hist(epsilon.gev)

param.vals <- c(400,200,150,.02)
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
steps <- c(.1, .02, .5, .02)

run1.mu.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = epsilon.gev, initial.steps = steps)

param.vals <- c(300,150, 200,.5)
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
steps <- c(.3, .2, .05, .2)

run2.mu.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = epsilon.gev, initial.steps = steps)

param.vals <- c(600,175, 300,.09)
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
steps <- c(.03, .2, .5, .01)

run3.mu.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = epsilon.gev, initial.steps = steps)

mcmc1 <- as.mcmc(run1.mu.nonstat$chains)
mcmc2 <- as.mcmc(run2.mu.nonstat$chains)
mcmc3 <- as.mcmc(run3.mu.nonstat$chains)

#compare the different values between the chains of alpha and beta
mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
plot(mcmc.chain.list)

gelman.plot(mcmc.chain.list)

plot(run1.mu.nonstat$chains[,1], type= 'l')
plot(run2.mu.nonstat$chains[,1], type = 'l')
plot(run3.mu.nonstat$chains[,1], type = 'l')

hist(run1.mu.nonstat$chains[1000:1e5,1])
abline(v=mean(run1.mu.nonstat$chains[1000:1e5,1]), col='red')

hist(run2.mu.nonstat$chains[1000:1e5,1])
abline(v=mean(run2.mu.nonstat$chains[1000:1e5,1]), col='red')

hist(run3.mu.nonstat$chains[1000:1e5,1])
abline(v=mean(run3.mu.nonstat$chains[1000:1e5,1]), col='red') #5% and 95% quantile lines also , also true value

#-----------now to try it with data!!!!

#-----------non stationary mu----------------------
param.vals <- c(400,200,150,.02)
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
steps <- c(.1, .02, .5, .02)

run1.mu.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,150, 200,.5)
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
steps <- c(.3, .2, .05, .2)

run2.mu.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600,175, 300,.09)
param.names <- c('mu0', 'mu1', 'sigma', 'xi')
steps <- c(.03, .2, .5, .01)

run3.mu.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

mcmc1 <- as.mcmc(run1.mu.nonstat$chains)
mcmc2 <- as.mcmc(run2.mu.nonstat$chains)
mcmc3 <- as.mcmc(run3.mu.nonstat$chains)

#compare the different values between the chains of alpha and beta
mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
plot(mcmc.chain.list)

gelman.plot(mcmc.chain.list)

plot(run1.mu.nonstat$chains[,1], type= 'l')
plot(run2.mu.nonstat$chains[,1], type = 'l')
plot(run3.mu.nonstat$chains[,1], type = 'l')

plot(run1.mu.nonstat$chains[,1], type= 'l')
plot(run2.mu.nonstat$chains[,1], type = 'l')
plot(run3.mu.nonstat$chains[,1], type = 'l')

hist(run1.mu.nonstat$chains[1000:1e5,1])
abline(v=mean(run1.mu.nonstat$chains[1000:1e5,1]), col='red')

hist(run2.mu.nonstat$chains[1000:1e5,1])
abline(v=mean(run2.mu.nonstat$chains[1000:1e5,1]), col='red')

hist(run3.mu.nonstat$chains[1000:1e5,1])
abline(v=mean(run3.mu.nonstat$chains[1000:1e5,1]), col='red')


#---------------nonstationary sigma------------------------------

param.vals <- c(400,125,150,.02)
param.names <- c('mu', 'sigma0', 'sigma1', 'xi')
steps <- c(.1, .02, .5, .02)

run1.sigma.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,150, 200,.5)
steps <- c(.3, .2, .05, .2)

run2.sigma.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600,175, 300,.09)
steps <- c(.03, .2, .5, .01)

run3.sigma.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)


#---------------non-stationary xi--------------------------------

param.vals <- c(400,125,1,.02)
param.names <- c('mu', 'sigma','xi0', 'xi1')
steps <- c(.1, .02, .25, .02)

run1.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,150, .1,.5)
steps <- c(.3, .2, .05, .2)

run2.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600,175, .1,.09)
steps <- c(.03, .2, .05, .01)

run3.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

#---------------non stationary mu and sigma----------------------

param.vals <- c(400,125,100,150,.02)
param.names <- c('mu0', 'mu1', 'sigma0','sigma1', 'xi1')
steps <- c(.1, .2,.02, .25, .02)

run1.mu.sigma.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,400,125, .1,.5)
steps <- c(.3, .2, .1, .2)

run2.mu.sigma.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600, 175, 400,100,.09)
steps <- c(.03, .2, .05, .2, .01)

run3.mu.sigma.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)


#---------------non stationary mu and xi ------------------------

param.vals <- c(400,125,100,150,.02)
param.names <- c('mu0', 'mu1', 'sigma0','sigma1', 'xi1')
steps <- c(.1, .2,.02, .25, .02)

run1.mu.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,400,150, .1,.5)
steps <- c(.3, .2, .05, .1, .2)

run2.mu.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600, 175, 400,.1,.09)
steps <- c(.03, .2, .05, .2, .01)

run3.mu.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)


#---------------non stationary sigma and xi----------------------
param.vals <- c(400,125,100,.1,.02)
param.names <- c('mu0', 'mu1', 'sigma0','sigma1', 'xi1')
steps <- c(.1, .2,.02, .25, .02)

run1.sigma.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,400,150, .1,.5)
steps <- c(.3, .2, .05, .1, .2)

run2.sigma.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600, 175, 400,.1,.09)
steps <- c(.03, .2, .05, .2, .01)

run3.sigma.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

#---------------non stationary mu sigma xi-----------------------

param.vals <- c(400,125,100,150,.02, .5)
param.names <- c('mu0', 'mu1', 'sigma0','sigma1','xi0', 'xi1')
steps <- c(.1, .2,.02, .25, .02, .5)

run1.mu.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(300,400,150,125, .1,.5)
steps <- c(.3, .2, .05, .1, .2)

run2.mu.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)

param.vals <- c(600, 175, 400,100,.1,.09)
steps <- c(.03, .2, .05, .2, .5,.01)

run3.mu.xi.nonstat <- mcmc.function.gev(niter = 1e5, params = param.vals, parnames = param.names, data = tide.data$max, initial.steps = steps)
