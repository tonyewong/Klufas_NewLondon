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

niter <- 1e5
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
initial.params <- c(500,100,100,.2)
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

initial.params <- c(500,100,75,.2, .05)
parnames <- c('mu0', 'sigma0', 'sigma1', 'xi0', 'xi1')
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

initial.params <- c(500, 1000,100,75,.2, .05)
parnames <- c('mu0', 'mu1', 'sigma0', 'sigma1', 'xi')
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


save.image('mcmc.parallel.RData')
