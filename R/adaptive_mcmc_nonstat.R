#Alex Klufas 
#June 28, 2017 
#adapt MCMC Stuff

library(adaptMCMC)
library(mvtnorm)
library(extRemes)
library(coda)

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data()

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
tide.data <- read.tide.data()

setwd('~/codes/Klufas_NewLondon/R/')
source('gev_nonstationary_MCMC.R')

niter <- 1e6
initial.params <- c(500,100,.2)
parnames <- c('mu', 'sigma', 'xi')

#took about 75 sec 
mcmc.test.stationary <- MCMC(log.post.final2, niter, initial.params, 
                             adapt=TRUE,
                             acc.rate = .44, 
                             n.start=3e3, 
                             parnames = parnames, 
                             data = tide.data$max , 
                             temps = temps$values)

#niter <- 1e4
mcmc.stationary.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                               n.chain = 3, 
                                               n.cpu = 1, 
                                               scale =c(.1,.3,.2), 
                                               adapt = TRUE, 
                                               acc.rate = .44, 
                                               gamma = .75, 
                                               list = TRUE, 
                                               n.start = 3e3, 
                                               packages = 'extRemes',
                                               parnames = parnames, 
                                               data = tide.data$max, 
                                               temps = temps$values)

#pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.stationary.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.stationary.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.stationary.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#53e3

#remove everything before burn-in val 
mcmc.stationary.parallel.post.burnin.1 <- mcmc.stationary.parallel[[1]]$samples[53e3:1e6,] 
best.params.stat <- which.max(mcmc.stationary.parallel[[1]]$log.p)


#niter <- 5e4
initial.params <- c(1000,300,100,.2)
parnames <- c('mu0', 'mu1', 'sigma', 'xi')
mcmc.mu.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                          n.chain = 3, 
                                          n.cpu = 1, 
                                          scale =c(.1,.2,.3,.2), 
                                          adapt = TRUE, 
                                          acc.rate = .44, 
                                          gamma = .75, 
                                          list = TRUE, 
                                          n.start = 3e3, 
                                          packages = 'extRemes',
                                          parnames = parnames, 
                                          data = tide.data$max, 
                                          temps = temps$values)

my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags.mu <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.mu.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.mu.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.mu.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.mu[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#53e3 - same place as one above 

mcmc.mu.nonstat.parallel.post.burnin.1 <- mcmc.mu.nonstat.parallel[[1]]$samples[53e3:1e6,] 
best.params.mu.nonstat <- which.max(mcmc.mu.nonstat.parallel[[1]]$log.p)

#niter <- 5e4
initial.params <- c(500,100,100,.2)
parnames <- c('mu', 'sigma0', 'sigma1', 'xi')
mcmc.sigma.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                          n.chain = 3, 
                                          n.cpu = 1, 
                                          scale =c(.1,.2,.3,.2), 
                                          adapt = TRUE, 
                                          acc.rate = .44, 
                                          gamma = .75, 
                                          list = TRUE, 
                                          n.start = 3e3, 
                                          packages = 'extRemes',
                                          parnames = parnames, 
                                          data = tide.data$max, 
                                          temps = temps$values)
#niter <- 1e6
my.seq <- seq(from = 3e3, to = niter, by = 1e3)
list.of.gelman.diags.sigma <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#173e3
mcmc.sigma.nonstat.parallel.post.burnin.1 <- mcmc.sigma.nonstat.parallel[[1]]$samples[173e3:1e6,] 
best.params.sigma.nonstat <- which.max(mcmc.sigma.nonstat.parallel[[1]]$log.p)


initial.params <- c(500,100,.01,.2)
parnames <- c('mu', 'sigma', 'xi0', 'xi1')
mcmc.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                             n.chain = 3, 
                                             n.cpu = 1, 
                                             scale =c(.1,.2,.3,.2), 
                                             adapt = TRUE, 
                                             acc.rate = .44, 
                                             gamma = .75, 
                                             list = TRUE, 
                                             n.start = 3e3, 
                                             packages = 'extRemes',
                                             parnames = parnames, 
                                             data = tide.data$max, 
                                             temps = temps$values)

my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags.xi <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.xi[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#173e3

mcmc.xi.nonstat.parallel.post.burnin.1 <- mcmc.xi.nonstat.parallel[[1]]$samples[173e3:1e6,] 
best.params.xi.nonstat <- which.max(mcmc.xi.nonstat.parallel[[1]]$log.p)

initial.params <- c(500, 1000,100,75,.2)
parnames <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi')
mcmc.mu.sigma.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                          n.chain = 3, 
                                          n.cpu = 1, 
                                          scale =c(.1,.05,.2,.3,.2), 
                                          adapt = TRUE, 
                                          acc.rate = .44, 
                                          gamma = .75, 
                                          list = TRUE, 
                                          n.start = 3e3, 
                                          packages = 'extRemes',
                                          parnames = parnames, 
                                          data = tide.data$max, 
                                          temps = temps$values)

my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags.mu.sigma <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.mu.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#363e3

mcmc.mu.sigma.nonstat.parallel.post.burnin.1 <- mcmc.mu.sigma.nonstat.parallel[[1]]$samples[363e3:1e6,] 
best.params.mu.sigma.nonstat <- which.max(mcmc.mu.sigma.nonstat.parallel[[1]]$log.p)

initial.params <- c(500, 1000,100,.075,.2)
parnames <- c('mu0', 'mu1', 'sigma', 'xi0','xi1')
mcmc.mu.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                                n.chain = 3, 
                                                n.cpu = 1, 
                                                scale =c(.1,.05,.2,.3,.2), 
                                                adapt = TRUE, 
                                                acc.rate = .44, 
                                                gamma = .75, 
                                                list = TRUE, 
                                                n.start = 3e3, 
                                                packages = 'extRemes',
                                                parnames = parnames, 
                                                data = tide.data$max, 
                                                temps = temps$values)

my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags.mu.xi <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.mu.xi[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#453e3

mcmc.mu.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[453e3:1e6,] 
best.params.mu.xi.nonstat <- which.max(mcmc.mu.xi.nonstat.parallel[[1]]$log.p)

#save.image('mcmc.params.RData')


#need to re run this one---------------------------------

initial.params <- c(500,100,75,.2, .05)
parnames <- c('mu', 'sigma0', 'sigma1', 'xi0', 'xi1')
mcmc.sigma.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                                n.chain = 3, 
                                                n.cpu = 1, 
                                                scale =c(.1,.05,.2,.3,.2), 
                                                adapt = TRUE, 
                                                acc.rate = .44, 
                                                gamma = .75, 
                                                list = TRUE, 
                                                n.start = 3e3, 
                                                packages = 'extRemes',
                                                parnames = parnames, 
                                                data = tide.data$max, 
                                                temps = temps$values)

#this is throwing an error I need to deal with later - Error in chol.default(W) : 
#the leading minor of order 1 is not positive definite 
my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags.sigma.xi <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.sigma.xi[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#453e3

mcmc.sigma.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[453e3:1e6,] 
best.params.sigma.xi.nonstat <- which.max(mcmc.sigma.xi.nonstat.parallel[[1]]$log.p)


initial.params <- c(500, 1000,100,75,.2, .05)
parnames <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi0', 'xi1')
mcmc.mu.sigma.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                                n.chain = 3, 
                                                n.cpu = 1, 
                                                scale =c(.1,.05,.2,.3,.2, .07), 
                                                adapt = TRUE, 
                                                acc.rate = .44, 
                                                gamma = .75, 
                                                list = TRUE, 
                                                n.start = 3e3, 
                                                packages = 'extRemes',
                                                parnames = parnames, 
                                                data = tide.data$max, 
                                                temps = temps$values)

my.seq <- seq(from = 3e3, to = niter, by = 1e4)
list.of.gelman.diags.mu.sigma.xi <- rep(0, length(my.seq))
pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
for (i in 1:length(my.seq)){   #53e3
  mcmc1 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
  mcmc2 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
  mcmc3 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
  list.of.gelman.diags.mu.sigma.xi[i] <- gelman.diag(mcmc.chain.list)[2]
  setTxtProgressBar(pb, i)
}
close(pb)
#453e3

mcmc.sigma.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[453e3:1e6,] 
best.params.mu.sigma.xi.nonstat <- which.max(mcmc.mu.sigma.xi.nonstat.parallel[[1]]$log.p)

save.image('mcmc.rerun.RData')


#save.image('mcmc.allnonstat.for.real.RData')

#save.image('mcmc.parallel.longer.iter.RData')

#setwd('~/codes/Klufas_NewLondon/R/')
#load('mcmc.allnonstat.for.real.RData')



#setwd('~/codes/Klufas_NewLondon/R/')
#load('mcmc.parallel.RData')


#niter <- 1e4
#initial.params <- c(500, 1000,100,75,.2, .05)
#parnames <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi0', 'xi1')
#mcmc.mu.sigma.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                               #    n.chain = 3, 
                                                 #  n.cpu = 1, 
                                                 #  scale =c(.1,.5,.2,.3,.2, .7), 
                                                 #  adapt = TRUE, 
                                                 #  acc.rate = .44, 
                                                 #  gamma = .75, 
                                                 #  list = TRUE, 
                                                  # n.start = 3e3, 
                                                 #  packages = 'extRemes',
                                                 #  parnames = parnames, 
                                                  # data = tide.data$max, 
                                                 #  temps = temps$values)

#dev.off()
#niter <- 1e4
#initial.params <- c(500,1000,100,75,.2,.4)
#parnames <- c('mu0','mu1','sigma0', 'sigma1', 'xi0', 'xi1')
#mcmc.text.all.nonstat <- MCMC.parallel(log.post.final2, niter, initial.params, 
                                    #   n.chain = 3, 
                                      # n.cpu = 1, 
                                    #   scale =c(.1,.5,.2,.3,.2, .7), 
                                     #  adapt = TRUE, 
                                     #  acc.rate = .44, 
                                     #  gamma = .75, 
                                     #  list = TRUE, 
                                     #  n.start = 3e3, 
                                     #  packages = 'extRemes',
                                      # parnames = parnames, 
                                      # data = tide.data$max, 
                                      # temps = temps$values)

#mcmc.all.nonstat.one.chain <- MCMC(log.post.final2, niter, initial.params, 
                                #  adapt=TRUE,
                                #  acc.rate = .44, 
                                 # n.start=3e3, 
                                 # parnames = parnames, 
                                 # data = tide.data$max , 
                                 # temps = temps$values)