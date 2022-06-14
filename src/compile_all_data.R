### Code to compile all data from all fish species and place in a single data frame
library(tidyr)
library(dplyr)
select <- dplyr::select
# load density from MEPS paper
load("data/raw/fish_density_tim_essington/year_effects.rdata")

fishnames <-
  c("sole",
    "rockfish",
    "hake",
    "pollock",
    "smelt",
    "perch",
    "ratfish",
    "herring")
fish.w.density <-
  c("sole", "hake", "pollock", "smelt", "perch",  "ratfish", "herring")
density.names <-
  c(
    "English sole",
    "rockfish",
    "Pacific whiting",
    "Walleye pollock",
    "smelt",
    "Shiner perch",
    "Spotted ratfish",
    "herring"
  )
# Column names in data frames that are NOT parasite counts

non.par.names <- c(
  "Year.Collected",
  "longjitt",
  "latjitt",
  "site",
  "Fish.ID",
  "Standard.Length..mm.",
  "Total.Length..mm.",
  "long",
  "lat",
  "temp_na_rm",
  "stemp",
  "Pb",
  "As",
  "Zn",
  "Ni",
  "V",
  "Cr",
  "Cu",
  "Ba",
  "Be",
  "Sig8_lignin",
  "Lamb8",
  "Bd.V_soil_biomarker",
  "fish.spc",
  "fish_density",
  "sfish_density"
)

non.par.names <- tolower(non.par.names)


# loop through fishnames, load CSV file for each
### Compile these into a single dataframe ####
# loop through fish species

data.ver <- "2022.02.21"
for (i in 1:length(fishnames)) {
  # load CSV file containing compiled data
  filename <-
    paste0("data/processed/for_analysis/",
           fishnames[i],
           "_",
           data.ver,
           ".csv")
  tmp.dat <- read.csv(file = filename,
                      header = TRUE)
  # remove annoying "X" column
  tmp.dat <- dplyr::select(tmp.dat,!any_of(c("X", "X.1")))
  
  
  # create a column indicating fish species name
  tmp.dat$fish.spc <- fishnames[i]
  # put in long format
  colnames(tmp.dat) <- tolower(colnames(tmp.dat))
  # need to correct two things.  One the Contracaecum_larvae_PV, needs to be
  # renamed Contracaecum_sp_all
  
  clPVindex <-
    which(grepl(colnames(tmp.dat), pattern = "contracaecum_larvae_pv"))
  if (length(clPVindex) > 0) {
    new.df.names <- colnames(tmp.dat)
    new.df.names[clPVindex] <- "contracaecum_sp_all"
    names(tmp.dat) <- new.df.names
  }
  
  ### Next, if there is any "unknown" subtaxa, make them species specific ####
  unknindex <- which(grepl(colnames(tmp.dat), pattern = "_unkn"))
  if (length(unknindex) > 0) {
    new.df.names <- colnames(tmp.dat)
    new.df.names[unknindex] <-
      paste0(new.df.names[unknindex], "_", fishnames[i])
    names(tmp.dat) <- new.df.names
  }
  
  ### replace fish_density with a standardized fish density ####
  if (fishnames[i] %in% fish.w.density) {
    if (!fishnames[i] %in% c("smelt", "herring")) {
      eval.text <- paste0("fish.dens <- y.list$`", density.names[i], "`")
      eval(parse(text = eval.text))
      fish.dens <- as.data.frame(fish.dens)
      fish.dens$year.collected <- seq(1948L, 1977L, by = 1L)
      fish.dens$sfish_density <- as.numeric(scale(fish.dens$median))
      fish.dens <-
        dplyr::select(fish.dens, year.collected, sfish_density)
    }
    if (fishnames[i] == "smelt") {
      fish.dens <-
        read.csv(
          "data/raw/fish_density_correigh_greene/cpuedata.csv",
          sep = ",",
          header = TRUE
        )
      fish.dens <- fish.dens %>%
        group_by(Year) %>%
        summarize_at(vars(Surf.smelt.CPUE), funs(mean = mean(., na.rm =
                                                               TRUE)))

      colnames(fish.dens) <- c("year.collected", "fish_density")
      fish.dens$sfish_density <-
        as.numeric(scale(fish.dens$fish_density))
      fish.dens <-
        dplyr::select(fish.dens, year.collected, sfish_density)
    }
    if (fishnames[i] == "herring") {
      fish.dens <- read.csv("data/raw/fish_density_correigh_greene/cpuedata.csv",sep=",",header=TRUE)
      fish.dens <- fish.dens %>%
        group_by(Year) %>%
        summarize_at(vars(Herring.CPUE), funs(mean = mean(., na.rm =
                                                               TRUE)))
      
      colnames(fish.dens) <- c("year.collected", "fish_density")
      fish.dens$sfish_density <-
        as.numeric(scale(fish.dens$fish_density))
      fish.dens <-
        dplyr::select(fish.dens, year.collected, sfish_density)
    }
    # merge this into tmp.data
    tmp.dat <-
      left_join(x = tmp.dat, y = fish.dens, by = "year.collected")
  }
  
  dfnames <- (colnames(tmp.dat))
  cols.2.remove <- c(non.par.names)
  
  parasite.names <- dfnames[!dfnames %in% cols.2.remove]
  long.tmp.dat <- pivot_longer(
    data = tmp.dat,
    cols = all_of(parasite.names),
    names_to = "para.spc",
    values_to = "count"
  )
  # scale length separately for each fish species
  long.tmp.dat$length <- scale(long.tmp.dat$total.length..mm.)
  
  
  # add this to existing data frame (or if i = 1, initialize data frame)
  if (i == 1) {
    all.dat.long <- long.tmp.dat
  } else {
    # annoying thing to fix issue where some data frames don't have a standard length column
    col.names <- colnames(long.tmp.dat)
    missing.cols <- setdiff(colnames(all.dat.long), col.names)
    # create this column if needed
    if (length(missing.cols) > 0) {
      for (j in 1:length(missing.cols)) {
        col.2.add <- missing.cols[j]
        eval.text <- paste0("long.tmp.dat$", col.2.add, "<-NA")
        eval(parse(text = eval.text))
      }
    }
    
    all.dat.long <- rbind(all.dat.long, long.tmp.dat)
  }
}

  all.dat.long <- all.dat.long %>%
    mutate(year = year.collected) %>%
    select(
      site,
      year,
      fish.id,
      length,
      long,
      lat,
      fish.spc,
      para.spc,
      count,
      sfish_density,
      temp_na_rm,
      pb,
      as,
      zn,
      ni,
      v,
      cr,
      cu,
      ba,
      be,
      sig8_lignin,
      lamb8,
      bd.v_soil_biomarker
    )
  
  
  
  ### Add species group and life history strategies ####
  
  # load parasite types file
  psite.types <-
    read.csv(file = "data/raw/relating_psite_names_to_codes/corresponding_psite_types_master.csv", header = T)
  all.dat.long$par.type <- NA
  
  # loop through parasite names, for each look up the parasite type, and place
  # that into data frame
  unique.para.spc <- tolower(unique(all.dat.long$para.spc))
  
  for (i in 1:length(unique.para.spc)) {
    tmp.para <- unique.para.spc[i]
    # find match in psite.types
    para.type.index <-
      which(sapply(
        X = tolower(psite.types$Model),
        FUN = grepl,
        pattern = tmp.para
      ))
    all.dat.long$par.type[which(all.dat.long$para.spc == tmp.para)] <-
      psite.types$psite_type[para.type.index[1]]
  }
  
  # there will be some with no matches because of unkn_all_fishspecies which don't appear
  na.list <- which(is.na(all.dat.long$par.type))
  phyl <- c("cest", "nem", "trem")
  phyl.types <- c("cestode", "nematode", "trematode")
  
  for (i in 1:length(na.list)) {
    # get parasite name
    par.name <- all.dat.long$para.spc[na.list[i]]
    for (j in 1:3) {
      if (grepl(par.name, pattern = phyl[j]))
        all.dat.long$par.type[na.list[i]] <- phyl.types[j]
    }
  }
  
  
  # finally, add in the life history groupings
  direct.transmit <- c("copepod", "leech", "monogenean")
  complex.life.cycle <-
    c("trematode", "cestode", "nematode", "acanthocephalan")
  
  all.dat.long$lifecycle <- NA
  direct.index <- which(all.dat.long$par.type %in% direct.transmit)
  complex.index <-
    which(all.dat.long$par.type %in% complex.life.cycle)
  all.dat.long$lifecycle[direct.index] <- "direct"
  all.dat.long$lifecycle[complex.index] <- "complex"
  
  # Now get the number of hosts
  all.dat.long$nhosts <- NA
  all.dat.long$terminal <- NA
  # load data
  host.data <- read.csv("data/raw/psite_life_cycle_complexity/psite_lc_complexity.csv", header = T)
  host.data$terminal[is.na(host.data$terminal)] <- 2
  for (i in 1:length(unique.para.spc)) {
    tmp.para <- unique.para.spc[i]
    # find match in psite.types
    para.type.index <-which(host.data$tim_code == tmp.para)
    all.dat.long$nhosts[which(all.dat.long$para.spc == tmp.para)] <-
     host.data$n_hosts[para.type.index[1]]
    all.dat.long$terminal[which(all.dat.long$para.spc == tmp.para)] <- host.data$terminal[para.type.index[1]]
  }
  
  
  ### Ordination on Contaminant Data ####
  ## Here, use PCoA or PCA (try both) on annual data on contaminants.
  yearly.dat <- all.dat.long %>%
    group_by(year) %>%
    summarise(
      pb = mean(pb),
      as = mean(as),
      zn = mean(zn),
      ni = mean(ni),
      v = mean(v),
      cr = mean(cr),
      cu = mean(cu),
      ba = mean(ba),
      be = mean(be),
      sig8_lignin = mean(sig8_lignin),
      lamb8 = mean(lamb8)
    ) %>%
    filter(!is.na(pb), !year %in% c(1977, 1978, 1991))
  

  # compare with PCA
  ord_predictor_PCA <-
    princomp(x = yearly.dat[, -1], scores = T, cor = T)
  biplot(ord_predictor_PCA, main = "ordination on contaminants")
  
  # PCA seems totally sensible here.  First two axes account for 74 percent of variance.
  
  ord_1 <- (ord_predictor_PCA$scores[, 1])
  ord_2 <- ord_predictor_PCA$scores[, 2]
  plot(yearly.dat$year, ord_1, type = "p", pch = 21, bg = "black")
  plot(yearly.dat$year, ord_2, type = "p", pch = 21, bg = "black")
  points(yearly.dat$year, smooth.spline(ord_1)$y, pch = 21, bg="red")
  loess_predict <- predict(loess(y ~ x, data.frame(x = yearly.dat$year, y = ord_1), span = 0.25))
  points(yearly.dat$year, loess_predict, pch = 21, bg = "blue")
  ord_1_loess <- loess_predict <- predict(loess(y ~ x, data.frame(x = yearly.dat$year, y = ord_1), span = 0.25))
  ord_2_loess <- loess_predict <- predict(loess(y ~ x, data.frame(x = yearly.dat$year, y = ord_2), span = 0.25))
  ## Now, need to put these ordination values back into long data. Use brute force
  
  all.dat.long$ord1 <- all.dat.long$ord2 <- all.dat.long$ord1raw <- all.dat.long$ord2raw <- NA
  for (i in 1:nrow(yearly.dat)) {
    year.2.use <- yearly.dat$year[i]
    year.index <- which(all.dat.long$year == year.2.use)
    all.dat.long$ord1[year.index] <- ord_1_loess[i]
    all.dat.long$ord2[year.index] <- ord_2_loess[i]
    all.dat.long$ord1raw[year.index] <- ord_1[i]
    all.dat.long$ord2raw[year.index] <- ord_2[i]
  }
  
  # now remove the contaminant data from long dataframe
  all.dat.long <- dplyr::select(all.dat.long,-pb,-as,-zn,-ni,-v ,-cr,-cu,-ba,-be,-sig8_lignin,-lamb8)
  
  
  # calculate rolling average temperature over decade time period with linearly declining weights
  
  library(smooth)
  yearly.dat <- all.dat.long %>%
    group_by(year) %>%
    summarise(
      temp = mean(temp_na_rm)
    ) %>%
    filter(!is.na(temp))
  y <- structure(yearly.dat$temp, 
                 .tsp = c(min(yearly.dat$year), max(yearly.dat$year), 1),
                 class = "ts")
  smoothtemp <- sma(y, h=0, silent=FALSE)
  yearly.dat$smooth <- smoothtemp$fitted
  all.dat.long$temp_smooth <- NA
  for (i in 1:nrow(yearly.dat)) {
    year.2.use <- yearly.dat$year[i]
    year.index <- which(all.dat.long$year == year.2.use)
    all.dat.long$temp_smooth[year.index] <- yearly.dat$smooth[i]
  }
  
  
  ### save into an RDS file for easy loading ####
  saveRDS(all.dat.long, file = "data/compiled_data.RDS")
  