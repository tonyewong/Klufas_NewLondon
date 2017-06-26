#Alex Klufas
#linear_fitting_MCMC.R
#Written on June 14, 2017 
#Modified on June 14, 2017
#
#
install.packages("mvtnorm")

install.packages("coda")

library(mvtnorm)
#Step 1: Fix alpha, beta, and vector x

alpha <- 2 

beta <- -1

x <- seq(from = 0, to = 100, by = 1)

#Step 2: Generate y star 

y.star <- alpha * x + beta 

plot(x, y.star, type='l')

#Step 3: Fix sigma (standard deviation) and generate data

sigma <- 5

epsilon <- rnorm(x, mean=0, sd = sigma)

y.func <- epsilon + y.star

points(x, y.func, col='red')

#Step 4: Find Likelihood for Theta-o using 'data points' that were just made 

#chose theta-o

#theta.0 <- c(4, -2)

#Need to find likelihood for this 

#log.like.theta.0 <- sum(dnorm(theta.0[1] * y.func + theta.0[2] , mean = 0, sd=sigma ,log = TRUE))

#theta.1 <- c(2,-1)

#log.like.theta.1 <- sum(dnorm(theta.1[1] * y.func + theta.1[2] , mean = 0, sd=sigma ,log = TRUE))

#need to find probability of certain theta 0 ? 
library(DEoptim)

neg.log.like <- function(p, sigma){ # p = theta.1, data = y.func 
  val <- sum(dnorm(p[1] * x + p[2] , mean = y.func, sd=sigma ,log = TRUE)) #make mu = 0 here 
  return(val)
}

log.like.pri <- function(p, mu, sigma, min, max){
  #should be for loop for all of the parameters in the function, will just start w two - alpha and beta 
  alpha <- p[1]
  beta <- p[2]
  
  #let alpha be normally distrtibuted 
  p.alpha <- sum(dnorm(x=alpha, mean = mu, sd = sigma, log = TRUE))
  #let beta be uniformly distributed 
  p.beta <- sum(dunif(x=beta, min = min, max = max, log = TRUE))
  
  #add together b/c logs
  p.alpha.beta <- (p.alpha + p.beta)
  return(p.alpha.beta)
}


log.post.final <- function(p, data, mu, sigma, min, max){
  
  log.like <- neg.log.like(p, sigma)
  
  log.like.pri <- log.like.pri(p, mu, sigma, min, max)
  
  log.like.final <- (log.like + log.like.pri)
  
  return(log.like.final)
}

#true alpha = 2, true beta = -1 

lower.bound <- c(0, -5)
upper.bound <- c(5, 10)

#optimization 
log.like.deoptim <- DEoptim(log.prior.final, lower=lower.bound, upper=upper.bound, data = y.func, mu = 2, sigma = 5)

#plotting 
lines(x, log.like.deoptim$optim$bestmem[1]*x + log.like.deoptim$optim$bestmem[2], col='blue')

#MCMC Stuff

#making vector that will be filled in with values 

#making this a function to use #note that mu, sigma, min, max are hard coded atm, will need to fix this for later 
mcmc.function <- function(niter, initial.alpha, initial.beta, data, initial.alpha.step, initial.beta.step){
  theta.vals <- cbind(rep(NA, niter), rep(NA, niter)) 

  n.theta.vals <- length(theta.vals[,1])

  theta.prob.vals <- rep(0, niter)

  #start with theta.1

  theta.1 <- c(initial.alpha, initial.beta)

  #calculate probability of this given theta set I chose 
  theta.1.prob <- log.post.final(p=theta.1, data=data, mu=2, sigma =5 , min=-5, max=5) #use y.func as parameter when using this function 

  #place theta values and probability values into vectors to old
  theta.vals[1,] <- cbind(theta.1[1], theta.1[2]) 
  theta.prob.vals[1] <- theta.1.prob

  param.changes <- rep(0,niter)
  #set stepsizes for alpha and beta 

  step.size.alpha <- initial.alpha.step

  step.size.beta <-  initial.beta.step
  
 # covar.vals <- rep(0,niter)

  sigma <- ((2.4)^2 /2) *matrix(c(step.size.alpha, 0, 0, step.size.beta), ncol=2)


#for loop 
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for (i in 2:niter){
  #going to try and use mvn to calculate alpha and beta
    theta.prev <- theta.vals[i-1,]
  
   if (i > 3e3){ #cov of the alpha and beta vectors 
      sigma2 <- ((2.4)^2 /2) * cov(theta.vals[1:(i-1),])
    #cov <- matrix(c(step.size.alpha, 0, 0, step.size.beta), ncol=2)
      covar.vals[i] <- sigma2

      mu <- c(theta.prev[1],theta.prev[2])
    
      alpha.beta.dependent <- rmvnorm(n = 1, mean = mu, sigma = sigma2) #sigma = covariance matrix
    
     alpha.star <- alpha.beta.dependent[1]
     beta.star <- alpha.beta.dependent[2]
     
    
   }else{
     covar.vals[i] <- sigma
     mu <- c(theta.prev[1],theta.prev[2])
     alpha.beta.dependent <- rmvnorm(n = 1, mean = mu, sigma = sigma) 
     alpha.star <- alpha.beta.dependent[1]
     beta.star <- alpha.beta.dependent[2]
    #propose theta.star values - need to chose an alpha and a beta 
    #alpha.star <- rnorm(1, mean =  theta.prev[1], sd = step.size.alpha)
    #beta.star <- runif(1, min = theta.prev[2] - step.size.beta, max = theta.prev[2] + step.size.beta)
    }
  
    current.theta <- c(alpha.star, beta.star)
  
    #calculate probability of new params
    theta.star.prob <- log.post.final(p=current.theta, data= y.func, mu= 0, sigma=5, min = -5, max =5 )
    
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
    #theta.vals[i,2] <- current.theta[2]
    #the probability that came out of this for loop 
   theta.prob.vals[i] <- current.theta.prob
  
  setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #need to return list of a lot of values from this function 
  mcmc.out <- vector('list', 4)
  names(mcmc.out) <- c('chains', 'acceptances', 'lpost', 'covar')
  
  mcmc.out$chains <- theta.vals
  mcmc.out$acceptances <- param.changes
  mcmc.out$lpost <- theta.prob.vals
  mcmc.out$covar <- covar.vals
  return(mcmc.out)
}


run1 <- mcmc.function(niter = 1e5, initial.alpha = 2, initial.beta = 4, data = y.func, initial.alpha.step =.5, initial.beta.step=.1)
run2 <- mcmc.function(niter= 1e5, initial.alpha = 1, initial.beta =-2, data= y.func, initial.alpha.step=.02, initial.beta.step = .34)
run3 <- mcmc.function(niter= 1e5, initial.alpha = 0, initial.beta =3 , data= y.func, initial.alpha.step = .75 , initial.beta.step = .05)

library(coda)

mcmc1 <- as.mcmc(run1$chains)
mcmc2 <- as.mcmc(run2$chains)
mcmc3 <- as.mcmc(run3$chains)

#compare the different values between the chains of alpha and beta 
mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
plot(mcmc.chain.list)

#from this graph, can see that we can igmore values before 1000
gelman.plot(mcmc.chain.list)

#--------i_burnin----------------------------------------------------
niter <- 1e5
my.seq <- seq(from = 1e3, to = niter, by = 1e3)
list.of.gelman.diags <- rep(0, length(my.seq))

for (i in 1:length(my.seq)){
  mcmc1 <- as.mcmc(run1$chains[1:my.seq[i],])
  mcmc2 <- as.mcmc(run2$chains[1:my.seq[i],])
  mcmc3 <- as.mcmc(run3$chains[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags[i] <- gelman.diag(mcmc.chain.list)[2]
}

#from this for loop - can tell that there is convergence after the 6000th iteration 

#plot the different alphas in the mcmcs between runs and compare 
plot(run1$chains[,1], type= 'l')
plot(run2$chains[,1], type = 'l')
plot(run3$chains[,1], type = 'l')


#alpha histograms 
hist(run1$chains[1000:1e5,1])
abline(v=mean(run1$chains[1000:1e5,1]), col='red')

hist(run2$chains[1000:1e5,1])
abline(v=mean(run2$chains[1000:1e5,1]), col='red')

hist(run3$chains[1000:1e5,1])
abline(v=mean(run3$chains[1000:1e5,1]), col='red') #5% and 95% quantile lines also , also true value 

hist(run1$chains[6000:1e5,1])
hist(run2$chains[6000:1e5,1])
hist(run3$chains[6000:1e5,1])

#comparing different chunks of convergence pieces 
hist(run1$chains[8e4:1e5,1])
abline(v=mean(run1$chains[8e4:1e5,1]), col='red')
hist(run1$chains[4e4:9e4,1])
abline(v=mean(run1$chains[4e4:9e4,1]), col='red')
hist(run1$chains[6e4:6.5e4,1])
abline(v=mean(run1$chains[6e4:6.5e4,1]), col='red')

hist(run1$chains[6000:1e5,1])
hist(run2$chains[6000:1e5,1])
hist(run3$chains[6000:1e5,1])

#plot the different betas in the mcmcs between runs and compare 
plot(run1$chains[,2], type= 'l')
plot(run2$chains[,2], type = 'l')
plot(run3$chains[,2], type = 'l')

#beta histograms

hist(run1$chains[1000:1e5,2])
abline(v=mean(run1$chains[1000:1e5,2]), col='red')

hist(run2$chains[1000:1e5,2])
abline(v=mean(run2$chains[1000:1e5,2]), col='red')

hist(run3$chains[1000:1e5,2])
abline(v=mean(run3$chains[1000:1e5,2]), col='red')

#plot(seq(from= 0,to = 1, by =1), gelman.diag.vals)

run1.alpha <- density(x=run1$chains[1e4:1e5,1])
run2.alpha <- density(x=run2$chains[1e4:1e5,1])
run3.alpha <- density(x=run3$chains[1e4:1e5,1])

plot(run1.alpha, col='red')
lines(run2.alpha, col='blue')
lines(run3.alpha, col='green')
legend(x=10,c("run1", "run2", "run3"), col = c('red', 'blue', 'green'))



#---------comparing density kernels of different chunks of convergence of alpha------------------------

run1.alpha <- density(x=run1$chains[8e4:1e5,1])
run2.alpha <- density(x=run2$chains[4e4:9e4,1])
run3.alpha <- density(x=run3$chains[6e4:6.5e4,1])

plot(run1.alpha, col='red')
lines(run2.alpha, col='blue')
lines(run3.alpha, col='green')

#------------------------------------------------------------------------------------------------------

run1.beta <- density(x=run1$chains[1e4:1e5, 2])
run2.beta <- density(x=run2$chains[1e4:1e5, 2])
run3.beta <- density(x=run3$chains[1e4:1e5, 2])

plot(run1.beta, col='red')
lines(run2.beta, col = 'blue')
lines(run3.beta, col='green')

#---------comparing density kernels of different chunks of convergence of alpha------------------------

run1.beta <- density(x=run1$chains[8e4:1e5,2])
run2.beta <- density(x=run2$chains[4e4:9e4, 2])
run3.beta <- density(x=run3$chains[6e4:6.5e4, 2])

plot(run1.beta, col='red')
lines(run2.beta, col = 'blue')
lines(run3.beta, col='green')

#------------------------------------------------------------------------------------------------------

#make a fake GEV 
library(extRemes)

#set params - true values 
mu <- 400

sigma  <- 150

xi <- .02 

x.gev <- seq(from=0,to=1000,by=5)

#the true value curve
stationary.curve <- devd(x.gev, loc = mu, scale = sigma, shape = xi, log=FALSE, type=c('GEV'))
plot(stationary.curve, type='l')
                      
#time to make some data

epsilon.gev <- revd(length(x.gev), loc=mu, scale = sigma, shape = xi, type=c("GEV"))

epsilon.gev.devd <- devd(epsilon.gev, loc = mu, scale = sigma, shape = xi, log=FALSE, type=c('GEV'))

x.gev.star <-  epsilon.gev.devd 

plot(x.gev.star)
lines(stationary.curve)

plot(epsilon.gev, epsilon.gev.devd )
lines(x.gev, stationary.curve)




neg.log.like2 <- function(p, data, loc, sigma, xi){ # p = theta.1, data = y.func 
  val <- sum(devd(data, loc= loc, scale = sigma, shape = xi, log=TRUE))#make mu = 0 here 
  return(val)
}

log.like.pri2 <- function(p, loc, sd, xi, min, max){
  #should be for loop for all of the parameters in the function, will just start w two - alpha and beta 
  mu <- p[1]
  sigma <- p[2]
  xi <- p[3]
  
  #let alpha be normally distrtibuted 
  p.mu <- sum(dnorm(x=mu, mean = loc, sd = sd, log = TRUE))
  #let beta be uniformly distributed 
  p.sigma <- sum(dunif(x=sigma, min = min, max = max, log = TRUE))
  
  p.xi <- sum(dunif(x=xi, min = min, max = max, log = TRUE))
  
  #add together b/c logs
  p.mu.sigma.xi <- p.mu  + p.sigma + p.xi
  return(p.mu.sigma.xi)
}


log.post.final2 <- function(p, data, loc, sd, xi, min, max){
  
  log.like <- neg.log.like2(p, data, loc, sd, xi)
  
  log.like.pri <- log.like.pri2(p, loc, sd, xi, min, max)
  
  log.like.final <- (log.like + log.like.pri)
  
  return(log.like.final)
}



mcmc.function.gev <- function(niter, initial.mu, initial.sigma, initial.xi, data, initial.mu.step , initial.sigma.step, initial.xi.step){
  theta.vals <- cbind(rep(NA, niter), rep(NA, niter), rep(NA, niter)) 
  
  n.theta.vals <- length(theta.vals[,1])
  
  theta.prob.vals <- rep(0, niter)
  
  #start with theta.1
  
  theta.1 <- c(initial.mu, initial.sigma, initial.xi)
  
  #calculate probability of this given theta set I chose 
  theta.1.prob <- log.post.final2(p=theta.1, data=data, loc=100, sd =50 , xi = .01, min=-5, max=5) #use y.func as parameter when using this function 
  
  #place theta values and probability values into vectors to old
  theta.vals[1,] <- cbind(theta.1[1], theta.1[2], theta.1[3]) 
  theta.prob.vals[1] <- theta.1.prob
  
  param.changes <- rep(0,niter)
  #set stepsizes for alpha and beta 
  
  step.size.mu <- initial.mu.step
  
  step.size.sigma <-  initial.sigma.step
  
  step.size.xi <- initial.xi.step 
  
  covar.vals <- rep(0,niter)
  
  sigma.sd <- ((2.4)^2 /2) *matrix(c(step.size.mu, 0, 0, 0, step.size.sigma , 0, 0, 0, step.size.beta), ncol=3)
  
  
  #for loop 
  pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
  for (i in 2:niter){
    #going to try and use mvn to calculate alpha and beta
    theta.prev <- theta.vals[i-1,]
    
    if (i > 3e3){ #cov of the alpha and beta vectors 
      sigma.sd2 <- ((2.4)^2 /2) *cov(theta.vals[1:(i-1),])

      #cov <- matrix(c(step.size.alpha, 0, 0, step.size.beta), ncol=2)
      covar.vals[i] <- sigma.sd2
      
      mu <- c(theta.prev[1],theta.prev[2], theta.prev[3])
      
      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = mu, sigma = sigma.sd2) #sigma = covariance matrix
      
      mu.star <- mu.sigma.xi.dependent[1]
      sigma.star <- mu.sigma.xi.dependent[2]
      xi.star <- mu.sigma.xi.dependent[3]
      
      
    }else{
      covar.vals[i] <- sigma.sd
      mu <- c(theta.prev[1],theta.prev[2], theta.prev[3])
      
      mu.sigma.xi.dependent <- rmvnorm(n = 1, mean = mu, sigma = sigma.sd) 
      
      mu.star <- mu.sigma.xi.dependent[1]
      sigma.star <- mu.sigma.xi.dependent[2]
      xi.star <- mu.sigma.xi.dependent[3]
    }
    
    current.theta <- c(mu.star, sigma.star, xi.star)
    
    #calculate probability of new params
    theta.star.prob <- log.post.final2(p=current.theta, data= data, loc=100, sd =50 , xi = .01, min = -5, max =5 )
    
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
    #theta.vals[i,2] <- current.theta[2]
    #the probability that came out of this for loop 
    theta.prob.vals[i] <- current.theta.prob
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #need to return list of a lot of values from this function 
  mcmc.out <- vector('list', 4) #chains, acceptances, lpost, covar
  names(mcmc.out) <- c('chains', 'acceptances', 'lpost', 'covar')
  
  mcmc.out$chains <- theta.vals
  mcmc.out$acceptances <- param.changes
  mcmc.out$lpost <- theta.prob.vals
  mcmc.out$covar <- covar.vals
  return(mcmc.out)
}
                        
gev.run1 <- mcmc.function.gev(niter = 1e5, initial.mu = 50, initial.sigma = 50, initial.xi = 1, data =epsilon.gev.devd , initial.mu.step = .2 , initial.sigma.step = 1, initial.xi.step = .05)                        



print("Number of New Params Accepted")
print(sum(param.changes))

print("Percent Accepted")
print(sum(param.changes) / niter)

plot(1:niter, theta.vals[,1], type='l')
plot(1:niter, theta.vals[,2], type='l')
plot(1:niter, theta.prob.vals, type='l')

hist(theta.vals[(niter / 2):niter,1])
hist(theta.vals[,2])
