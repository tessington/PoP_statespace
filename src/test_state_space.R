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

### Make fake data for testing ####
rho <- 0.5
phi <- 1
sigma_v <- 0.1
logno <- 0


# make fake betas for length, lat and long
betas <- runif(n = 3, -1, 1)
# make fake gammas for taxa
gammas <- runif(n = length(host.2.keep), -1, 1)

#### Do year indexing ####
min.year <- min(spc.data$year)
max.year <- max(spc.data$year)
yearlist <- min.year:max.year
year.lookup <- function(x,y) which(y == x)
yearindex <- sapply(X = spc.data$year, FUN = year.lookup, y = yearlist)

# Initialize parameters to store lognt and rt
lognt <- NULL
r_t <- NULL

#### Calculate population dynamics ####
lognt[1] <- logno
r_t[1] <- rnorm(1,0, sigma_v)

for (i in 2:length(yearlist)) {
  v_t <- rnorm(1, 0, sigma_v)
  r_t[i] <- rho * r_t[i-1] + sqrt(1 - rho^2) * v_t
  lognt[i] <- lognt[i-1] + r_t[i- 1]
}

#### Create X and U ####
Xij <- model.matrix(~ -1 + slat + slong + length, data = spc.data)
Uij <- model.matrix(~ -1 + fish.spc, data = spc.data)

#### Create simulated counts from model ####
eta_ij <- Xij %*% betas + Uij %*% gammas + lognt[yearindex]
ndata <- length(eta_ij)
count <- rnbinom(n = ndata,  mu = exp(eta_ij), size = phi)

### Setup TMB ####

tmb_data <- list(n = ndata,
                 nyears = length(yearlist),
                 y = count,
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

#### Compare fixed effects to true ####
fixef <- summary(rep, "fixed")
raef <- summary(rep, "random")
reef <- summary(rep, "report")

get_est <- function(ef, parname) return(ef[grep(rownames(ef), pattern = paste0("\\b",parname, "\\b")),])

beta_est <- get_est(fixef, "beta")
cbind(beta_est, betas)

#### Compare random effects to true
gamma_est <- get_est(raef, "gamma")
cbind(gamma_est, gammas)

#### Compare population size to true
ntEst <- get_est(reef, "nt")
cbind(ntEst, exp(lognt))


