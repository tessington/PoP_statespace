### File to setup and configure data for TMB model runs
#### Create Design Matrixes #####
# for convenience, make separate design matrix for intercepts (all par / fish combos)
Xint <- model.matrix(~ -1 + hostpar, data = thedata)
Xg <- model.matrix(~ -1 + nhosts, data = thedata)
Xa <- model.matrix(~ -1 + par.type, data = thedata)
Xgt <- model.matrix(~ -1 + nhost_term, data = thedata)

# Now make the fixed effects, need to do this in two steps to add in parasite-level fish dependence

# separate fish length effects for each parasite / fish combo
Xlength <- Xint * as.numeric(thedata$length)

Xlat <- model.matrix(~ -1 +  para.spc, data = thedata) * as.numeric(thedata$slat) 
Xlong <- model.matrix(~ -1 + para.spc, data = thedata) * as.numeric(thedata$slong)


X <- cbind(Xlength, Xlat)

# random effects design matrix - parasite specific year effects
U <- model.matrix(~ -1 + para.spc, data = thedata) * as.numeric(scale(thedata$year))
A <- model.matrix(~ -1 + par.type, data = thedata)* as.numeric(scale(thedata$year))
Z <- model.matrix(~ -1 + fish.id, data = thedata)

#### Set up other data ####
spc <- levels(thedata$para.spc)
nspc <- length(spc)
ngrp <- ncol(A)

# create a vector of parasite species ID for each observation
Xspc <- (as.numeric(thedata$para.spc) -1)
spc.groups <- (levels(thedata$par.type))
# look up Life history (LH) for each group
LH.tmp <- rep(NA, length(spc.groups))
for (i in 1:length(spc.groups)) {
  LH.tmp[i] <- dplyr::filter(thedata, par.type == spc.groups[i])$lifecycle[1]
}
y <- as.integer(thedata$count)
# subtract 1 for C++ indexing
LH <- LH.tmp -1

# get index of taxonomic group for each species, and create species specific LH type
# and number of hosts in life cycle
tax.spc.tmp <-LH.spc.tmp <- NH.spc.tmp <- NHterm.spc.tmp <- rep(NA, nspc)
for (i in 1:length(spc)) {
  tax.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$par.type[1]
  LH.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$lifecycle[1]
  NH.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$nhosts[1]
  NHterm.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$nhost_term[1]
}

# subtract 1 for C++ indexing

LHspc <- as.integer(LH.spc.tmp -1)
NHspc <- as.integer(NH.spc.tmp -1)
NHtermspc <- as.integer(NHterm.spc.tmp - 1)