

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data()

source('read_tide_data.R')
tide.data <- read.tide.data()

source('gev_fitting_MCMC.R')

#doing GEV MCMC fitting for non stationary models 

#leggo

log.like.calc.nonstat <- function(p, parnames, data, temps){ #non stationary 
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
      
      p.sigma <- exp(sigma0+ sigma1*temps)
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
  
  #let alpha be normally distrtibuted 
  p.mu <- sum(dunif(x=mu, min= 0, max=3000, log = TRUE))
  #let beta be uniformly distributed 
  p.sigma <- sum(dunif(x=sigma, min= 0, max=3000, log = TRUE))
  
  p.xi <- sum(dunif(x=xi,min= -10, max=10, log = TRUE))
  
  #add together b/c logs
  p.mu.sigma.xi <- p.mu  + p.sigma + p.xi
  return(p.mu.sigma.xi)
}

log.post.final2 <- function(p, parnames, data, min, max, temps){
  
  log.like.pri.nonstat <- log.like.pri.nonstat(p, parnames, min, max)
  
  #log.like.pri <- log.like.pri2(p, loc, sd, xi, min, max)
  
  if (is.finite(log.like.pri)){
    log.like <- log.like.calc.nonstat(p, parnames, min, max, temps)
    
    log.like.final <- (log.like + log.like.pri)
  }
  else{
    log.like.final <-log.like.pri  
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

priors$max <- 600

mcmc.function.gev <- function(niter, params, parnames, data, temps, initial.steps){
  #theta.vals <- cbind(rep(NA, niter), rep(NA, niter), rep(NA, niter)) 
  n.parnames <- length(parnames)
  
  #creating matrix that is niter by a certain number of parameters large 
  theta.vals <- matrix(nrow=niter, ncol=n.parnames)
  
  #theta.vals <- mat.or.vec( niter, 3)
  
  #making array to store the probabilities of theta params 
  theta.prob.vals <- rep(0, niter)
  
  #start with theta.1
  for (i in 1:n.parnames){
    theta.vals[1,i] <- params[i]
  }
  #theta.1 <- c(initial.mu, initial.sigma, initial.xi)
  theta.1 <- theta.vals[,1]
  #calculate probability of this given theta set I chose 
  theta.1.prob <- log.post.final2(p=theta.1, 
                                  parnames = parnames,
                                  data = tide.data$max, 
                                  temps = temps$values, 
                                  min=priors$min, 
                                  max=priors$max) #use y.func as parameter when using this function 
  
  #place theta values and probability values into vectors to old
  theta.prob.vals[1] <- theta.1.prob
  
  param.changes <- rep(0,niter)
  
  sigma.sd <- ((2.4)^2 /n.parnames) *matrix.maker(n.params, initial.steps)
  
  
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
      
      # covar.vals[i] <- sigma.sd
     # mu <- c(theta.prev[1],theta.prev[2], theta.prev[3])
      
      #mu.sigma.xi.dependent <- revd(n=1, loc = theta.prev[1], scale = theta.prev[2], shape =theta.prev[3] , type=c("GEV"))
      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = theta.prev, sigma = sigma.sd) 
      #print(mu.sigma.xi.dependent)
      
      
    }
    
    current.theta <- mu.sigma.xi.dependent
    
    theta.star.prob <- log.post.final2(p = current.theta, 
                                       parnames = parnames,
                                       data = tide.data$max, 
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

mu <- 400

mu0 <- 200 

sigma  <- 150

xi <- .02 

x.gev <- seq(from=0,to=75,by=1)

simple.curve <- rep(0, length(temps$values))
#the true value curve
for (i in 1:length(temps$values)){
  simple.curve[i] <- sum(devd(x.gev, loc = mu + mu0*temps$values[i], scale = sigma, shape = xi, log=FALSE, type=c('GEV')))
}

plot(simple.curve)

#time to make some data

epsilon.gev <- rep(0, length(temps$values))

for(i in 1:length(temps$values)){
  epsilon.gev[i] <- sum(revd(length(x.gev), loc=mu + mu0*temps$values[i], scale = sigma, shape = xi, type=c("GEV")))
}

plot(x.gev, epsilon.gev)

epsilon.gev.devd <- rep(0, length(temps$values))

for(i in 1:length(temps$values)){

  epsilon.gev.devd[i] <- devd(epsilon.gev, loc = mu + mu0*temps$values[i], scale = sigma, shape = xi, log=FALSE, type=c('GEV'))

}
plot(epsilon.gev, epsilon.gev.devd )
