# Alex Klufas
# Written on June 28, 2017
# Modified on August 3, 2017
# 
# Code to perform MCMC calculations to get posterior distributions 
# Takes a long time to run 


library(adaptMCMC)
library(mvtnorm)
library(extRemes)
library(coda)

setwd('~/codes/Klufas_NewLondon/R/')
source('read_temp_data.R')

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')

setwd('~/codes/Klufas_NewLondon/R/')
source('gev_nonstationary_MCMC.R')

niter <- 1e6
set.seed(10)
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
                                               acc.rate = .235,
                                               gamma = .75,
                                               list = TRUE,
                                               n.start = 3e3,

                                               packages = 'extRemes',
                                               parnames = parnames,
                                               data = tide.data$max,
                                               temps = temps$values)

initial.params <- c(1000,300,100,.2)
parnames <- c('mu0', 'mu1', 'sigma', 'xi')
mcmc.mu.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
                                          n.chain = 3,
                                          n.cpu = 1,
                                          scale =c(.1,.2,.3,.2),
                                          adapt = TRUE,
                                          acc.rate = .235,
                                          gamma = .75,
                                          list = TRUE,
                                          n.start = 3e3,
                                          packages = 'extRemes',
                                          parnames = parnames,
                                          data = tide.data$max,
                                          temps = temps$values)


initial.params <- c(500,100,100,.2)
parnames <- c('mu', 'sigma0', 'sigma1', 'xi')
mcmc.sigma.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
                                             n.chain = 3,
                                             n.cpu = 1,
                                             scale =c(.1,.2,.3,.2),
                                             adapt = TRUE,
                                             acc.rate = .235,
                                             gamma = .75,
                                             list = TRUE,
                                             n.start = 3e3,
                                             packages = 'extRemes',
                                             parnames = parnames,
                                             data = tide.data$max,
                                             temps = temps$values)

initial.params <- c(500,100,.01,.2)
parnames <- c('mu', 'sigma', 'xi0', 'xi1')
mcmc.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
                                          n.chain = 3,
                                          n.cpu = 1,
                                          scale =c(.1,.2,.3,.2),
                                          adapt = TRUE,
                                          acc.rate = .235,
                                          gamma = .75,
                                          list = TRUE,
                                          n.start = 3e3,
                                          packages = 'extRemes',
                                          parnames = parnames,
                                          data = tide.data$max,
                                          temps = temps$values)

initial.params <- c(500, 1000,100,75,.2)

parnames <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi')
mcmc.mu.sigma.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
                                                n.chain = 3,
                                                n.cpu = 1,
                                                scale =c(.1,.05,.2,.3,.2),
                                                adapt = TRUE,
                                                acc.rate = .235,
                                                gamma = .75,
                                                list = TRUE,
                                                n.start = 3e3,
                                                packages = 'extRemes',
                                                parnames = parnames,
                                                data = tide.data$max,
                                                temps = temps$values)

initial.params <- c(500,20,75,.2, .05)
parnames <- c('mu', 'sigma0', 'sigma1', 'xi0', 'xi1')
mcmc.sigma.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
                                                n.chain = 3,
                                                n.cpu = 1,
                                                scale =c(.1,.05,.2,.3,.2),
                                                adapt = TRUE,
                                                acc.rate = .235,
                                                gamma = .75,
                                                list = TRUE,
                                                n.start = 3e3,
                                                packages = 'extRemes',
                                                parnames = parnames,
                                                data = tide.data$max,
                                                temps = temps$values)


initial.params <- c(500, 1000,100,75,.2, .05)
parnames <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi0', 'xi1')
mcmc.mu.sigma.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
                                                   n.chain = 3,
                                                   n.cpu = 1,
                                                   scale =c(.1,.05,.2,.3,.2, .07),
                                                   adapt = TRUE,
                                                   acc.rate = .235,
                                                   gamma = .75,
                                                   list = TRUE,
                                                   n.start = 3e3,
                                                   packages = 'extRemes',
                                                   parnames = parnames,
                                                   data = tide.data$max,
                                                   temps = temps$values)


# Code below is structure to get burnin value-----------------------------------------

#pb <- txtProgressBar(min=0,max=niter,initial=0,style=3)
#my.seq <- seq(from = 3e3, to = niter, by = 1e3)
#list.of.gelman.diags <- rep(0, length(my.seq))
#pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
#for (i in 1:length(my.seq)){   #53e3
# mcmc1 <- as.mcmc(mcmc.stationary.parallel[[1]]$samples[1:my.seq[i],])
#  mcmc2 <- as.mcmc(mcmc.stationary.parallel[[2]]$samples[1:my.seq[i],])
##  mcmc3 <- as.mcmc(mcmc.stationary.parallel[[3]]$samples[1:my.seq[i],])
#  mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
 # list.of.gelman.diags[i] <- gelman.diag(mcmc.chain.list)[2]
 # setTxtProgressBar(pb, i)
#}
#close(pb)
#49e4


#remove everything before burn-in val 
#mcmc.stationary.parallel.post.burnin.1 <- mcmc.stationary.parallel[[1]]$samples[49e3:5e4,] 
#best.params.stat <- which.max(mcmc.stationary.parallel[[1]]$log.p)
#save.image('mcmc.stat.post.burnin.RData')

#remove everything before burn-in val
# mcmc.stationary.parallel.post.burnin.1 <- mcmc.stationary.parallel[[1]]$samples[53e3:1e6,]
# best.params.stat <- which.max(mcmc.stationary.parallel[[1]]$log.p)




#my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
#53e3 - same place as one above 

# mcmc.mu.nonstat.parallel.post.burnin.1 <- mcmc.mu.nonstat.parallel[[1]]$samples[53e3:1e6,] 
# best.params.mu.nonstat <- which.max(mcmc.mu.nonstat.parallel[[1]]$log.p)

# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #53e3 - same place as one above
# 
# mcmc.mu.nonstat.parallel.post.burnin.1 <- mcmc.mu.nonstat.parallel[[1]]$samples[53e3:1e6,]
# best.params.mu.nonstat <- which.max(mcmc.mu.nonstat.parallel[[1]]$log.p)



# my.seq <- seq(from = 3e3, to = niter, by = 1e3)
# list.of.gelman.diags.sigma <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #173e3
# mcmc.sigma.nonstat.parallel.post.burnin.1 <- mcmc.sigma.nonstat.parallel[[1]]$samples[173e3:1e6,] 
# best.params.sigma.nonstat <- which.max(mcmc.sigma.nonstat.parallel[[1]]$log.p)

# my.seq <- seq(from = 3e3, to = niter, by = 1e3)
# list.of.gelman.diags.sigma <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #173e3
# mcmc.sigma.nonstat.parallel.post.burnin.1 <- mcmc.sigma.nonstat.parallel[[1]]$samples[173e3:1e6,]
# best.params.sigma.nonstat <- which.max(mcmc.sigma.nonstat.parallel[[1]]$log.p)

# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #173e3
# 
# mcmc.xi.nonstat.parallel.post.burnin.1 <- mcmc.xi.nonstat.parallel[[1]]$samples[173e3:1e6,] 
# best.params.xi.nonstat <- which.max(mcmc.xi.nonstat.parallel[[1]]$log.p)

#initial.params <- c(500, 100,50,25,.2)
# =======
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #173e3
# 
# mcmc.xi.nonstat.parallel.post.burnin.1 <- mcmc.xi.nonstat.parallel[[1]]$samples[173e3:1e6,]
# best.params.xi.nonstat <- which.max(mcmc.xi.nonstat.parallel[[1]]$log.p)


# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.sigma <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
#363e3



# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #453e3
# 
# mcmc.mu.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[453e3:1e6,] 
# best.params.mu.xi.nonstat <- which.max(mcmc.mu.xi.nonstat.parallel[[1]]$log.p)
# 
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #453e3
# 
# mcmc.mu.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[453e3:1e6,]
# best.params.mu.xi.nonstat <- which.max(mcmc.mu.xi.nonstat.parallel[[1]]$log.p)


#save.image('mcmc.params.RData')


#need to re run this one---------------------------------
#niter <- 1e6
# image.save('mcmc.sigma.xi.RData')


# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.sigma.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.sigma.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# 
# mcmc.sigma.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[453e3:1e6,]
# best.params.sigma.xi.nonstat <- which.max(mcmc.sigma.xi.nonstat.parallel[[1]]$log.p)

# initial.params <- c(500, 50,100,.075,.2)
# parnames <- c('mu0', 'mu1', 'sigma', 'xi0','xi1')
# mcmc.mu.xi.nonstat.parallel <- MCMC.parallel(log.post.final2, niter, initial.params,
#                                              n.chain = 3,
#                                              n.cpu = 1,
#                                              scale =c(.1,.05,.2,.3,.2),
#                                              adapt = TRUE,
#                                              acc.rate = .235,
#                                              gamma = .75,
#                                              list = TRUE,
#                                              n.start = 3e3,
#                                              packages = 'extRemes',
#                                              parnames = parnames,
#                                              data = tide.data$max,
#                                              temps = temps$values)


# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.sigma.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.sigma.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
#453e3


mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[763e3:1e6,]
best.params.mu.sigma.xi.nonstat <- which.max(mcmc.mu.sigma.xi.nonstat.parallel[[1]]$log.p)


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

mcmc.mu.sigma.nonstat.parallel.post.burnin.1 <- mcmc.mu.sigma.nonstat.parallel[[1]]$samples[383e3:1e6,]
best.params.mu.sigma.xi.nonstat <- which.max(mcmc.mu.sigma.nonstat.parallel[[1]]$log.p)





#stationary
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.stationary.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.stationary.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.stationary.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #53e3
# 
# #remove everything before burn-in val 
# mcmc.stationary.parallel.post.burnin.1 <- mcmc.stationary.parallel[[1]]$samples[33e3:1e6,] 
# best.params.stat <- which.max(mcmc.stationary.parallel[[1]]$log.p)
# 
# #mu non stat
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #73e3 - same place as one above 
# 
# mcmc.mu.nonstat.parallel.post.burnin.1 <- mcmc.mu.nonstat.parallel[[1]]$samples[53e3:1e6,] 
# best.params.mu.nonstat <- which.max(mcmc.mu.nonstat.parallel[[1]]$log.p)
# 
# #sigma non stat 
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.sigma <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #293e3
# mcmc.sigma.nonstat.parallel.post.burnin.1 <- mcmc.sigma.nonstat.parallel[[1]]$samples[283e3:1e6,] 
# best.params.sigma.nonstat <- which.max(mcmc.sigma.nonstat.parallel[[1]]$log.p)
# 
# #xi nonstat 
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #113e3
# 
# mcmc.xi.nonstat.parallel.post.burnin.1 <- mcmc.xi.nonstat.parallel[[1]]$samples[93e3:1e6,] 
# best.params.xi.nonstat <- which.max(mcmc.xi.nonstat.parallel[[1]]$log.p)
# 
# #mu & sigma non stat #NEEED TO COME BACK HERE
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.sigma <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.sigma.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.sigma[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #433e3
#
# mcmc.mu.sigma.nonstat.parallel.post.burnin.1 <- mcmc.mu.sigma.nonstat.parallel[[1]]$samples[433e3:1e6,]
# best.params.mu.sigma.nonstat <- which.max(mcmc.mu.sigma.nonstat.parallel[[1]]$log.p)

#mu & xi non stat 
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #123e3
# 
# mcmc.mu.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[63e3:1e6,] 
# best.params.mu.xi.nonstat <- which.max(mcmc.mu.xi.nonstat.parallel[[1]]$log.p)
# 
# #sigma and xi non stat 
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.sigma.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.sigma.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.sigma.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #453e3
# 
# mcmc.sigma.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.xi.nonstat.parallel[[1]]$samples[93e3:1e6,] 
# best.params.sigma.xi.nonstat <- which.max(mcmc.sigma.xi.nonstat.parallel[[1]]$log.p)
# 
# #all non stat 
# my.seq <- seq(from = 3e3, to = niter, by = 1e4)
# list.of.gelman.diags.mu.sigma.xi <- rep(0, length(my.seq))
# pb <- txtProgressBar(min=0,max=length(my.seq),initial=0,style=3)
# for (i in 1:length(my.seq)){   #53e3
#   mcmc1 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[1:my.seq[i],])
#   mcmc2 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[2]]$samples[1:my.seq[i],])
#   mcmc3 <- as.mcmc(mcmc.mu.sigma.xi.nonstat.parallel[[3]]$samples[1:my.seq[i],])
#   mcmc.chain.list <- mcmc.list(list(mcmc1, mcmc2, mcmc3))
#   list.of.gelman.diags.mu.sigma.xi[i] <- gelman.diag(mcmc.chain.list)[2]
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# #573e3
# 
# mcmc.mu.sigma.xi.nonstat.parallel.post.burnin.1 <- mcmc.mu.sigma.xi.nonstat.parallel[[1]]$samples[573e3:1e6,] 
# best.params.mu.sigma.xi.nonstat <- which.max(mcmc.mu.sigma.xi.nonstat.parallel[[1]]$log.p)
