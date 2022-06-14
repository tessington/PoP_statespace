# Script to load and set up data for analysis
thedata <- readRDS("data/compiled_data.RDS")
make_host_par <- function(x) paste(x[1],"-", x[2])
thedata$hostpar <- apply(FUN = make_host_par, MAR = 1, X = cbind(as.character(thedata$fish.spc), as.character(thedata$para.spc)))
thedata$hostpar <- as.factor(thedata$hostpar)                          
makezeroone <- function(x) min(1,x)
thedata$foo <- sapply(X = thedata$count, FUN = makezeroone)

#### Optional - remove rare taxa ######
foo.min <- 0.04
mean.hostpar <- thedata %>%
  group_by(hostpar) %>%
  summarise(ntotal = n(), npresent = sum(foo), meancount = mean(count))
mean.hostpar$pfoo <- mean.hostpar$npresent / mean.hostpar$ntotal

hostpar.2.keep <- mean.hostpar$hostpar[mean.hostpar$pfoo>=foo.min]

thedata <- thedata %>%
  filter(hostpar %in% hostpar.2.keep)

thedata$hostpar <- as.factor(as.character(thedata$hostpar))
thedata$para.spc <- as.factor(thedata$para.spc)
thedata$fish.spc <- as.factor(thedata$fish.spc)
thedata$par.type <- as.factor(thedata$par.type)
thedata$lifecycle <- as.factor(thedata$lifecycle)
thedata$fish.id <- as.factor(thedata$fish.id)


# pool hosts 3 and 4
thedata$nhosts[thedata$nhosts == 4] <- 3
thedata$nhosts <- as.factor(thedata$nhosts)

thedata$nhost_term <- apply(FUN = make_host_par, MAR = 1, X = cbind(as.character(thedata$nhosts), as.character(thedata$terminal)))
thedata$nhost_term <- as.factor(thedata$nhost_term)
thedata$slat <- scale(thedata$lat)
thedata$slat[is.na(thedata$slat)] <- 0
thedata$slong <- scale(thedata$long)
thedata$slong[is.na(thedata$slong)] <- 0
# get the SD of year, used to convert slopes to per year
year.scale <- sd(thedata$year)

