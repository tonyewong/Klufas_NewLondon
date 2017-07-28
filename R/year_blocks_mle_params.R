#working with different year blocks for estimating params 

setwd('~/codes/Klufas_NewLondon/R/')

source('read_temp_data.R')
temps <- read.temp.data(1939, 2014)

setwd('~/codes/Klufas_NewLondon/R/')
source('read_tide_data.R')
tide.data <- read.tide.data()

setwd('~/codes/Klufas_NewLondon/R/')
source('optimization_sf.R')


all.p.names <- list(c('mu', 'sigma', 'xi'), 
                    c('mu0', 'mu1', 'sigma', 'xi'),
                    c('mu', 'sigma0', 'sigma1', 'xi'),
                    c('mu', 'sigma', 'xi0', 'xi1'),
                    c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi'),
                    c('mu', 'sigma0', 'sigma1' , 'xi0', 'xi1'),
                    c('mu0', 'mu1', 'sigma', 'xi0', 'xi1'),
                    c('mu0', 'mu1', 'sigma0', 'sigma1' , 'xi0', 'xi1'))

other.p.names <- c('stat', 'mu','sigma', 'xi', 'mu.sigma', 'sigma.xi', 'mu.xi', 'mu.sigma.xi')
lower.bound <- list(c(0,0,-5),
                    c(0,-100, 0, -5),
                    c(0, 0, 0, -200),
                    c(0, 0, -5, -5), 
                    c(0,-100, -100, -100, -5),
                    c(0, -100, -100, -1, -1),
                    c(0, -100, 0, -1, -10),
                    c(0,-100, 0, -100, -1, -1)
)
upper.bound <- list(c(3000,3000,5),
                    c(3000,1000, 1000, 5),
                    c(3000, 1000, 1000, 1000),
                    c(3000, 1000 , 100, 100),
                    c(3000,100, 1000, 10 , 5),
                    c(3000, 1000, 10 , 1, 1),
                    c(3000, 100, 1000, 10, 10),
                    c(3000,100, 1000, 10 , 1, 1)
)



get.prior.estimates <- function(city, city.temps){
  
  sublist <- vector('list', 8)
  
  names(sublist) <- other.p.names
  for (j in 1: length(all.p.names)){
    print(j)
    optim.like <- DEoptim(neg.log.like.calc,
                          lower=lower.bound[[j]], 
                          upper=upper.bound[[j]], 
                          data = city, 
                          temps = city.temps, 
                          parnames = all.p.names[[j]])
    sublist[[other.p.names[j]]] <- optim.like$optim$bestmem
    
  }

  return(sublist)
}

#-----------------------------------------------------------------------------------
new.london.10 <- get.prior.estimates(new.london$max[56:66], new.london.temps$values[56:66])
new.london.20 <- get.prior.estimates(new.london$max[46:66], new.london.temps$values[46:66])
new.london.30 <- get.prior.estimates(new.london$max[36:66], new.london.temps$values[36:66])
new.london.40 <- get.prior.estimates(new.london$max[26:66], new.london.temps$values[26:66])
new.london.50 <- get.prior.estimates(new.london$max[16:66], new.london.temps$values[16:66])
new.london.all <- get.prior.estimates(new.london$max, new.london.temps$values)

boston.10 <- get.prior.estimates(boston$max[81:91], boston.temps$values[81:91])
boston.20 <- get.prior.estimates(boston$max[71:91], boston.temps$values[71:91])
boston.30 <- get.prior.estimates(boston$max[61:91], boston.temps$values[61:91])
boston.40 <- get.prior.estimates(boston$max[51:91], boston.temps$values[51:91])
boston.50 <- get.prior.estimates(boston$max[41:91], boston.temps$values[41:91])
boston.60 <- get.prior.estimates(boston$max[31:91], boston.temps$values[31:91])
boston.70 <- get.prior.estimates(boston$max[21:91], boston.temps$values[21:91])
boston.all <- get.prior.estimates(boston$max, boston.temps$values)

atlantic.city.10 <-  get.prior.estimates(atlantic.city$max[81:91], atlantic.city.temps$values[81:91])
atlantic.city.20 <-  get.prior.estimates(atlantic.city$max[71:91], atlantic.city.temps$values[71:91])
atlantic.city.30 <-  get.prior.estimates(atlantic.city$max[61:91], atlantic.city.temps$values[61:91])
atlantic.city.40 <-  get.prior.estimates(atlantic.city$max[51:91], atlantic.city.temps$values[51:91])
atlantic.city.50 <-  get.prior.estimates(atlantic.city$max[41:91], atlantic.city.temps$values[41:91])
atlantic.city.60 <-  get.prior.estimates(atlantic.city$max[31:91], atlantic.city.temps$values[31:91])
atlantic.city.70 <-  get.prior.estimates(atlantic.city$max[21:91], atlantic.city.temps$values[21:91])
atlantic.city.all <-  get.prior.estimates(atlantic.city$max, atlantic.city.temps$values)

portland.10 <- get.prior.estimates(portland$max[80:90], portland.temps$values[80:90])
portland.20 <- get.prior.estimates(portland$max[70:90], portland.temps$values[70:90])
portland.30 <- get.prior.estimates(portland$max[60:90], portland.temps$values[60:90])
portland.40 <- get.prior.estimates(portland$max[50:90], portland.temps$values[50:90])
portland.50 <- get.prior.estimates(portland$max[40:90], portland.temps$values[40:90])
portland.60 <- get.prior.estimates(portland$max[30:90], portland.temps$values[30:90])
portland.70 <- get.prior.estimates(portland$max[20:90], portland.temps$values[20:90])
portland.all <- get.prior.estimates(portland$max, portland.temps$values)
#save.image('10.year.gaps2.RData')

