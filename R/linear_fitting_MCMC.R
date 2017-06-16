#Alex Klufas
#linear_fitting_MCMC.R
#Written on June 14, 2017 
#Modified on June 14, 2017
#
#

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
niter <- 1e6
theta.vals <- cbind(rep(NA, niter), rep(NA, niter)) 

n.theta.vals <- length(theta.vals[,1])

theta.prob.vals <- rep(0, niter)

#start with theta.1

theta.1 <- c(4.0, 2.0)

#calculate probability of this given theta set I chose 
theta.1.prob <- log.post.final(p=theta.1, data=y.func, mu=2, sigma =5 , min=-5, max=5)

#place theta values and probability values into vectors to old
theta.vals[1,] <- cbind(theta.1[1], theta.1[2]) 
theta.prob.vals[1] <- theta.1.prob

param.changes <- rep(0,niter)
#set stepsizes for alpha and beta 

step.size.alpha <- .025

step.size.beta <- .1

#for loop 
pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
for (i in 2:niter){
  theta.prev <- theta.vals[i-1,]
  #propose theta.star values - need to chose an alpha and a beta 
  alpha.star <- rnorm(1, mean = 2, sd = step.size.alpha)
  beta.star <- runif(1, min = theta.prev[2] - step.size.beta, max = theta.prev[2] + step.size.beta)
  
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
  }else
  #if greater than or equal to alpha.prob, prev theta value assigned, theta params from before also replace current theta values 
   {
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

print("Number of New Params Accepted")
print(sum(param.changes))

print("Percent Accepted")
print(sum(param.changes) / niter)

plot(1:niter, theta.vals[,1], type='l')
plot(1:niter, theta.vals[,2], type='l')
plot(1:niter, theta.prob.vals, type='l')

hist(theta.vals[(niter / 2):niter,1])
hist(theta.vals[,2])
