library(tidyr)
library(TMB)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)

### Setup Data ####
source("src/load_data.R")
source("src/configure_TMB_runs.R")        

spc.data <- dplyr::filter(thedata, para.spc == "contracaecum_sp_all")
### Select host fish species commonly infected ####
spc.data$fish.spc <- as.character(spc.data$fish.spc)

foo.min <- 0.04
mean.host <- thedata %>%
  group_by(fish.spc) %>%
  summarise(ntotal = n(), npresent = sum(foo), meancount = mean(count))
mean.host$pfoo <- mean.host$npresent / mean.host$ntotal

host.2.keep <- mean.host$fish.spc[mean.host$pfoo>=foo.min]

spc.data <- spc.data %>%
  filter(fish.spc %in% host.2.keep)

spc.data$fish.spc <- as.factor(spc.data$fish.spc)


#### Do year indexing ####
min.year <- min(spc.data$year)
max.year <- max(spc.data$year)
yearlist <- min.year:max.year
year.lookup <- function(x,y) which(y == x)
yearindex <- sapply(X = spc.data$year, FUN = year.lookup, y = yearlist)



#### Create X and U ####
Xij <- model.matrix(~ -1 + slat + slong + length, data = spc.data)
Uij <- model.matrix(~ -1 + fish.spc, data = spc.data)


### Setup TMB ####
ndata <- nrow(Xij)
tmb_data <- list(n = ndata,
                 nyears = length(yearlist),
                 y = spc.data$count,
                 Xij = Xij,
                 Uij = Uij,
                 yearindex = yearindex - 1
)

tmb_pars <- list( beta = rep(0, 3),
                  gamma = rep(0, length(host.2.keep)),
                  logsigma_gamma = 0,
                  vt = rep(0, ndata -1),
                  logitrho = 0,
                  logsigma_vt = 0,
                  logno = 0,
                  logphi = 0
)

model <- "state_space"
compile(paste0("src/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)

#### Extract fixed effects ####
fixef <- summary(rep, "fixed")
raef <- summary(rep, "random")
reef <- summary(rep, "report")

get_est <- function(ef, parname) return(ef[grep(rownames(ef), pattern = paste0("\\b",parname, "\\b")),])

beta_est <- get_est(fixef, "beta")


#### Extract random effects ####
gamma_est <- get_est(raef, "gamma")
print(gamma_est)

#### Exrtract population size ####
ntEst <- get_est(reef, "nt")
plot(yearlist, ntEst[,1])

#### Extract population growth rate ####
rtEst <- get_est(reef, "rt")



