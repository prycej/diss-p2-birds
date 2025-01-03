library(data.table)
library(reshape2)
library(stats4)
library(rjags)
library(runjags)
library(dplyr)
library(BayesTools)
library(coda)
library(remotes)
library(fitR)
library(MCMCvis)
library(dplyr)
library(tidyr)
library(stringr)
library(bayesplot)
library(ggplot2)
library(zoo) # for converting to date object
library(lubridate) # for extracting month
library(jagsUI)

blackcap_data_raw <- fread("CR.blackcap_FixRing.csv", header = TRUE)
chiffchaff_data_raw <- fread("CR.chifchaf_FixRing.csv", header = TRUE)
robin_data_raw <- fread("CR.robin_FixRing.csv", header = TRUE)


clean_dataset <- function(df, species) {
  df <- df %>% mutate(rn = row_number(), species = species) %>%
    mutate(id = paste0(species, "_", rn)) %>% select(-rn) # add id_column made up of the species
  # and the row number from the originating dataset
  df <- df[,5:(ncol(df))]  # remove first 4 months (as they are only a partial
  # winter)
  df <- df %>% filter(rowSums(select(., -c(id, age_at_ringing, species))) > 0) # remove
  # rows where there are no capture information
  df
}

(blackap_data_cleaned <- clean_dataset(blackcap_data_raw, 'blackcap'))
(chiffchaff_data_cleaned <- clean_dataset(chiffchaff_data_raw, 'chiffchaff'))
(robin_data_cleaned <- clean_dataset(robin_data_raw, 'robin'))

birds <- rbind(blackap_data_cleaned, chiffchaff_data_cleaned,
               robin_data_cleaned)

##############################
# EDA
##############################

# custom theme for plotting purposes
custom_theme <- theme(
  axis.text.x = element_text(size = size / 2),
  axis.text.y = element_text(size = size / 2),
  axis.title.x = element_text(size = size / 2, margin = margin(t = 10)), # 10 units
  # of space above the x-axis title
  axis.title.y = element_text(size = size / 2, margin = margin(r = 10)), # 10 units
  # of space to the right of the y-axis title
  plot.title = element_text(size = size / 2, hjust = 0.5)
)

###### Age at first ringing proportions by species
ggplot(birds, aes(fill=age_at_ringing, x=species)) +
  geom_bar(position="fill") +
  coord_flip() +
  scale_fill_manual(values = c("Unknown" = "grey", "adult" = "dark red",
                               "juvenile" = "pink")) +
  labs(x = "Species", y = "Proportion", fill = "Age at first ringing") +
  custom_theme +
  theme(legend.text = element_text(size = size / 2),
        legend.title = element_text(size = size / 2))
ggsave("report/images/age_props.png")


#### Capture count proportions by species
birds$row_sum <- rowSums(select_if(birds, is.numeric))
# Convert the row_sum to a factor
birds$row_sum <- as.factor(birds$row_sum)
# Plot the data
ggplot(birds, aes(fill=row_sum, x=species)) +
  geom_bar(position="fill") +
  scale_fill_discrete(name = "Row Sum") +
  coord_flip() +
  labs(x = "Species", y = "Proportion") +
  custom_theme +
  theme(legend.text = element_text(size = size / 2),
        legend.title = element_text(size = size / 2))
ggsave("report/images/capture_count_proportions.png")


# Pivot and preprocess data
birds_long <- birds %>%
  pivot_longer(cols = starts_with("20"), names_to = "date", values_to = "count") %>%
  mutate(
    year = as.numeric(substr(date, 1, 4)),
    month = as.numeric(substr(date, 5, 6)),
    winter_period = ifelse(month >= 10, year, year - 1)
  ) %>%
  mutate(
    month_name = factor(ifelse(month >= 10, month - 9, month + 3),
                        levels = c(1:7),
                        labels = c("October", "November", "December",
                                   "January", "February", "March", "April"))
  )


# Capture count by winter period
for (species_idx in c("blackcap", "chiffchaff", "robin")) {
  winter_sums <- birds_long %>%
    filter(species == species_idx) %>%
    group_by(winter_period, age_at_ringing) %>%
    summarise(count = sum(count, na.rm = TRUE))

  p_winter <- ggplot(winter_sums, aes(x = factor(winter_period), y = count,
                                      fill = age_at_ringing)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(name = "Winter Period", labels = function(x) paste0("W", as.numeric(x) - 1)) +
    labs(y = "Count", fill = "Age Category") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom", legend.text = element_text(size = size / 2),
          legend.title = element_text(size = size / 2))  # Move legend to the bottom

  print(p_winter)
  ggsave(paste0("report/images/", species_idx, "_capture_counts_per_year.png"), p_winter)
}

# Capture count across months by age at ringing and for each species
for (species_idx in c("blackcap", "chiffchaff", "robin")) {
  monthly_sums <- birds_long %>%
    filter(species == species_idx) %>%
    group_by(month_name, age_at_ringing) %>%
    summarise(count = sum(count, na.rm = TRUE))

  p_month <- ggplot(monthly_sums, aes(x = month_name, y = count,
                                      fill = age_at_ringing)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Month", y = "Count", fill = "Age Category") +
    scale_x_discrete(limit = c("October", "November", "December", "January",
                               "February", "March", "April")) +
    custom_theme +
    theme(legend.text = element_text(size = size / 2),
          legend.title = element_text(size = size / 2), axis.text.x =
            element_text(angle = 45, hjust = 1), legend.position = "bottom")  # Move legend to the bottom

  print(p_month)
  ggsave(paste0("report/images/", species_idx, "_capture_count_per_month.png"),
         plot = p_month)
}

##############################
# INTER-WINTER MODELLING
##############################

create_inter_winter_dataset <- function(df) {
    # collapses datset by creating a column winter_presence which calcualtes
    # if a bird was capture at all during that winter
    df %>% pivot_longer(-c(id, age_at_ringing, species, row_sum), names_to = "month",
                 values_to = "value") %>%
    mutate(year = as.integer(str_sub(month, start = 1, end = 4)),
           month_num = as.integer(str_sub(month, start = 5, end = 6))) %>%
    mutate(winter = ifelse(month_num >= 10, year, year - 1)) %>%
    group_by(id, winter, age_at_ringing, species) %>%
    summarise(winter_presence = ifelse(any(value == 1), 1, 0)) %>%
    pivot_wider(names_from = winter, values_from = winter_presence) %>%
    # sort by id column
    mutate(numeric_id = as.numeric(str_extract(id, "\\d+"))) %>%
    arrange(species, numeric_id) %>%
    select(-numeric_id) %>% ungroup
}

birds_inter <- create_inter_winter_dataset(birds)

create_intra_winter_dataset <- function(df) {
  df_long <- birds %>%
    pivot_longer(cols = starts_with("20"),
                 names_to = "year_month",
                 values_to = "value") %>%
    separate(year_month, into = c("year", "month"), sep = 4)

  df_long$year <- as.numeric(df_long$year)
  df_long$month <- as.numeric(df_long$month)

  # Creating a new variable for winter
  df_long <- df_long %>%
    mutate(winter = ifelse(month >= 10, year, year - 1)) %>%
    mutate(winter_month = ifelse(month >= 10, month - 9, month + 3)) %>%
    select(-c(month, year))

  # Pivoting data so each winter has its own column
  df_wide <- df_long %>%
    pivot_wider(names_from = winter_month, values_from = value)  %>%
    filter(rowSums(.[,c('1', '2', '3', '4', '5', '6', '7')]) != 0) %>%
    group_by(id) %>% ## where the id column is duplicated because of a capture
    # appearing multiple times for different years, take the first year
    slice_min(row_number(), n = 1) %>%
    ungroup() %>%
    # sort by id column
    mutate(numeric_id = as.numeric(str_extract(id, "\\d+"))) %>%
    arrange(species, numeric_id) %>%
    select(-numeric_id)

  # Renaming columns 1..7 to be Oct..Apr
  df_wide <- df_wide %>%
    rename(`Oct` = `1`,
           `Nov` = `2`,
           `Dec` = `3`,
           `Jan` = `4`,
           `Feb` = `5`,
           `Mar` = `6`,
           `Apr` = `7`)

  df_wide
}

birds_intra <- create_intra_winter_dataset(birds)

######### survival probabilities inter
blackcap_inter <- birds_inter %>% filter(species == 'blackcap') %>% ungroup
robin_inter <- birds_inter %>% filter(species == 'robin') %>% ungroup
chiffchaff_inter <- birds_inter %>% filter(species == 'chiffchaff') %>% ungroup


blackcap_intra <- birds_intra %>% filter(species == 'blackcap') %>% ungroup
robin_intra <- birds_intra %>% filter(species == 'robin') %>% ungroup
chiffchaff_intra <- birds_intra %>% filter(species == 'chiffchaff') %>% ungroup

###modelling
#nz <-5000
#yaug <-rbind(blackcap_inter[,4:ncol(blackcap_inter)],
#             array(0,dim=c(nz,ncol(blackcap_inter[,4:ncol(blackcap_inter)]))),
#             se.names = FALSE)
##yaug <-rbind(blackcap_data,array(0,dim=c(nz,ncol(blackcap_data))), use.names = FALSE)
#
#
#jags.data <-list(yaug=yaug,M=nrow(yaug),T=ncol(yaug))
##inits_m0 <- function() list(z=rep(1,nrow(yaug)),p=runif(1,0,1))
#inits_mt <- function() list(z=rep(1,nrow(yaug)),p=runif(ncol(yaug),0,1))
#params <-c("N","p","omega")
#
#
#model <- jags.model("blackcap_mt.jags", data = jags.data, inits=inits_mt,
#                    n.chains = 3)
#update(model, nb)
#saveRDS(model, "jag_runs/blackcap_mt_all_months_model")
#samples <- coda.samples(model, variable.names = params, n.iter = ni,  thin = nt)
#saveRDS(samples, paste("jag_runs/blackcap_mt_all_months"))
#plot(samples)
#
############# using binomial likelhood on row sums
#blackcap_intra <- birds %>% filter(species == 'blackcap')
#create_abundance_jags_model <- function(df) {
# y <- rowSums(df)
# naug <- 5000
# yaug <- c(y, rep(0, naug))
# ( M <- length(yaug) )

# jagsData <- list(y = yaug, n = ncol(df), M = M,
#                  w = ifelse(yaug > 0, 1, NA))
# wanted <- c("N")

# outM0 <- jags(jagsData, NULL, wanted,
#               model="closescaptures.jags", DIC=FALSE,
#               n.chains=nc, n.iter=ni*5, n.thin=nt, n.burnin=nb)

  #library(mcmcOutput)
  #outB <- mcmcOutput(outM0, default='N')
  #diagPlot(outB)

  #plot(outB)

#  outM0

#}
#
#model_w2007 <- create_abundance_jags_model(birds[,1:7])
#model_w2008 <- create_abundance_jags_model(birds[,8:14])
#model_w2009 <- create_abundance_jags_model(birds[,15:21])
#model_w2010 <- create_abundance_jags_model(birds[,22:28])

#MCMCplot(blackcap_m0_age, params = c("phi.juv0", "phi.ad0", "p0", "mean.phi.diff"))
#MCMCplot(chiffchaff_m0_age, params = c("phi.juv0", "phi.ad0", "p0", "mean.phi.diff"))
#MCMCplot(robin_m0_age, params = c("phi.juv0", "phi.ad0", "p0", "mean.phi.diff"))
#
#MCMCplot(blackcap_mt_age, params = c("phi.diff"))
#MCMCplot(chiffchaff_mt_age, params = c("phi.diff"))
#MCMCplot(robin_mt_age, params = c("phi.diff"))
#
#library(ggplot2)
#
#
#
#

set.seed(42)

ni <-10000
nt <-3
nb <-2000
nc <-3
na <- 1000

marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow =
                      n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] -
      sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

## Choose inits at 3 difference parts of the parameter space to minimise chance
## of lucky starting values
inits_m0_inter <- function(chain){
  list(
    phi.juv0 = runif(1),
    phi.ad0 = runif(1),
    p0 = runif(1))
}

inits_mt_age_inter <- function(chain){
  list(
    phi.juv = runif(10, 0, 1),
    phi.ad = runif(10, 0, 1),
    p = runif(10, 0, 1))
}

inits_m0_intra <- function(chain){
  list(
    phi0 = runif(1),
    p0 = runif(1))
}

inits_mt_intra <- function(chain){
  list(
    phi = runif(6, 0, 1),
    p = runif(6, 0, 1))
}

jags_model_multinomial_age_inter <- function(df, hyper_par, inits,
                                       model_type = "mt",
                                       params =  c("phi.juv", "phi.ad", "p")) {

  CH.J <- df %>% filter(age_at_ringing == "juvenile") %>%
    select(c(-id, -age_at_ringing, -species))
  CH.A <- df %>% filter(age_at_ringing == "adult") %>%
    select(c(-id, -age_at_ringing, -species))

  cap <- apply(CH.J, 1, sum)
  ind <- which(cap >= 2)
  CH.J.R <- CH.J[ind,] # Juvenile CH recaptured at least once
  CH.J.N <- CH.J[-ind,] # Juvenile CH never recaptured
  # Remove first capture
  first <- numeric()
  for (i in 1:dim(CH.J.R)[1]){
    first[i] <- min(which(CH.J.R[i,]==1))
  }
  CH.J.R1 <- CH.J.R
  for (i in 1:dim(CH.J.R)[1]){
    CH.J.R1[i,first[i]] <- 0
  }
  # Add grown-up juveniles to adults and create m-array
  CH.A.m <- rbind(CH.A, CH.J.R1)
  marr.a <- marray(CH.A.m)
  # Create CH matrix for juveniles, ignoring subsequent recaptures
  second <- numeric()
  for (i in 1:dim(CH.J.R1)[1]){
    second[i] <- min(which(CH.J.R1[i,]==1))
  }
  CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
  for (i in 1:dim(CH.J.R)[1]){
    CH.J.R2[i,first[i]] <- 1
    CH.J.R2[i,second[i]] <- 1
  }
  # Create m-array for these
  CH.J.R.marray <- marray(CH.J.R2)
  CH.J.R.marray[,dim(CH.J)[2]] <- 0
  # Create the m-array for juveniles never recaptured and add it to the
  #previous m-array
  CH.J.N.marray <- marray(CH.J.N)
  marr.j <- CH.J.R.marray + CH.J.N.marray

  releases.juv <- rowSums(marr.j)
  releases.ad <- rowSums(marr.a)
  jdata <- list(marr.j = marr.j, marr.a = marr.a,
                R.juv = releases.juv, R.ad = releases.ad,
                n.occasions=ncol(marr.j),
                alpha_phi_j = hyper_par$alpha_phi_j,
                beta_phi_j = hyper_par$beta_phi_j,
                alpha_phi_ad = hyper_par$alpha_phi_j,
                beta_phi_ad = hyper_par$beta_phi_ad,
                alpha_p = hyper_par$alpha_p,
                beta_p = hyper_par$beta_p)


  out <- jags(jdata, inits, params,
              paste0("models/multinomial_",model_type,"_age_inter.jags"), DIC=TRUE,
              n.chains=nc, n.adapt=na, n.iter=ni, parallel = TRUE,  seed = 42)
  out
}

jags_model_multinomial_intra <- function(df, hyper_par, inits,
                                       model_type = "mt",
                                       params =  c("phi.juv", "phi.ad", "p")) {

  marr <-  marray(df %>% select(-c(id, age_at_ringing, species, winter, row_sum)))
  releases <- rowSums(marr)
  jdata <- list(marr=marr,
                R = releases,
                n.occasions=ncol(marr),
                alpha_phi= hyper_par$alpha_phi,
                beta_phi = hyper_par$beta_phi,
                alpha_p = hyper_par$alpha_p,
                beta_p = hyper_par$beta_p)


  out <- jags(jdata, inits, params,
              paste0("models/multinomial_",model_type,"_intra_old.jags"),
              DIC=TRUE, n.chains=nc, n.adapt=na, n.iter=ni, parallel = TRUE,
              seed = 42)
  out
}


##############################
# PRIOR DISTRIBUTIONS
##############################
# Define parameters for the Beta distributions
params <- list(c(8, 11.58), c(8, 17.33), c(6, 3.93), c(3, 7))
names <- c("blackcaps", "chiffchaffs", "robins", "capture_probability")
colors <- c("red", "blue", "green", "black")

# Generate data and plots
for(i in 1:4){
  x <- seq(0, 1, length.out = 100)
  y <- dbeta(x, params[[i]][1], params[[i]][2])
  df <- data.frame(x = x, y = y)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(color = colors[i]) +
    theme_minimal() +
    custom_theme +
    #  ggtitle(names[i]) +
    xlab("Probability") +
    ylab("Density")

  # Save the plot
  ggsave(filename = paste0("report/images/prior_distribution_for_", names[i],
                           ".png"), plot = p, width = 6, height = 4)
}



###test for sensitivity analysis
alpha_phi_j <- 1
beta_phi_j <- 1
alpha_phi_ad <- 1
beta_phi_ad <-  1
alpha_p <- 1
beta_p <- 1

##############################
# INTER-WINTER MODELLING
##############################

alpha_phi_j_blackcap <- 8
beta_phi_j_blackcap <- 11.59
alpha_phi_ad_blackcap <- 8
beta_phi_ad_blackcap <-  11.59
alpha_p <- 3
beta_p <- 7


blackcap_m0_age_inter <- jags_model_multinomial_age_inter(blackcap_inter,
                                        inits = inits_m0_inter,
                                        model_type = "m0",
                                        hyper_par = list(
                                        alpha_phi_j = alpha_phi_j_blackcap,
                                        beta_phi_j = beta_phi_j_blackcap,
                                        lpha_phi_ad= alpha_phi_ad_blackcap,
                                        beta_phi_ad= beta_phi_ad_blackcap,
                                        alpha_p = alpha_p,
                                        beta_p = beta_p),
                                        params = c("phi.juv0",
                                                   "phi.ad0",
                                                   "p0", "mean.phi.diff",
                                                   "fit.juv", "fit.new.juv",
                                                   "fit.ad", "fit.new.ad"))

blackcap_m0_age_inter
MCMCplot(blackcap_m0_age_inter, params = c("phi.juv0", "phi.ad0", "p0",
                                           "mean.phi.diff"))
print(xtable(as.data.frame(blackcap_m0_age_inter$summary[,-11]),
   caption = "Blackcap $M_0$ inter-winter model diagnostic and summary output",
   label = "tab:blackcap_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/blackcap_m0_age_inter.tex")

alpha_phi_j_chiffchaff <- 8
beta_phi_j_chiffchaff  <- 17.33
alpha_phi_ad_chiffchaff  <- 8
beta_phi_ad_chiffchaff  <-  17.33
alpha_p_chiffchaff  <- 3
beta_p_chiffchaff  <- 7


chiffchaff_m0_age_inter <- jags_model_multinomial_age_inter(chiffchaff_inter,
                                      inits = inits_m0_inter,
                                      model_type = "m0",
                                      hyper_par = list(
                                      alpha_phi_j = alpha_phi_j_chiffchaff,
                                      beta_phi_j = beta_phi_j_chiffchaff,
                                      alpha_phi_ad= alpha_phi_ad_chiffchaff,
                                      beta_phi_ad= beta_phi_ad_chiffchaff,
                                      alpha_p = alpha_p,
                                      beta_p = beta_p),
                                      params = c("phi.juv0", "phi.ad0",
                                                 "p0", "mean.phi.diff",
                                                 "fit.juv",
                                                 "fit.new.juv",
                                                 "fit.ad",
                                                 "fit.new.ad"))

chiffchaff_m0_age_inter
MCMCplot(chiffchaff_m0_age_inter, params = c("phi.juv0", "phi.ad0", "p0",
                                             "mean.phi.diff"))
print(xtable(as.data.frame(chiffchaff_m0_age_inter$summary[,-11]),
 caption = "Chiffchaff $M_0$ inter-winter model diagnostic and summary output",
 label = "tab:chiffchaff_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/chiffchaff_m0_age_inter.tex")
alpha_phi_j_robin <- 6
beta_phi_j_robin <- 3.93
alpha_phi_ad_robin <- 6
beta_phi_ad_robin <- 3.93
alpha_p_robin <- 3
beta_p_robin <- 7


robin_m0_age_inter <- jags_model_multinomial_age_inter(robin_inter,
                                           inits = inits_m0_inter,
                                           model_type = "m0",
                                           hyper_par = list(
                                           alpha_phi_j = alpha_phi_j_robin,
                                           beta_phi_j = beta_phi_j_robin,
                                           lpha_phi_ad= alpha_phi_ad_robin,
                                           beta_phi_ad= beta_phi_ad_robin,
                                           alpha_p = alpha_p,
                                           beta_p = beta_p),
                                           params = c("phi.juv0",
                                                      "phi.ad0",
                                                      "p0", "mean.phi.diff",
                                                      "fit.juv", "fit.new.juv",
                                                      "fit.ad", "fit.new.ad"))

robin_m0_age_inter
print(xtable(as.data.frame(robin_m0_age_inter$summary[,-11]),
     caption = "Robin $M_0$ inter-winter model diagnostic and summary output",
     label = "tab:robin_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/robin_m0_age_inter.tex")
MCMCplot(robin_m0_age_inter, params = c("phi.juv0", "phi.ad0", "p0",
                                        "mean.phi.diff"))

blackcap_mt_age_inter <- jags_model_multinomial_age_inter(blackcap_inter,
                                          inits = inits_mt_age_inter,
                                          model_type = "mt",
                                          hyper_par = list(
                                          alpha_phi_j = alpha_phi_j_blackcap,
                                          beta_phi_j = beta_phi_j_blackcap,
                                          lpha_phi_ad= alpha_phi_ad_blackcap,
                                          beta_phi_ad= beta_phi_ad_blackcap,
                                          alpha_p = alpha_p,
                                          beta_p = beta_p),
                                          params = c("phi.juv", "phi.ad",
                                                     "p", "phi.diff",
                                                     "fit.juv", "fit.new.juv",
                                                     "fit.ad", "fit.new.ad"))
(model <- blackcap_mt_age_inter)
MCMCplot(blackcap_mt_age_inter, params = c("phi.diff"))
print(xtable(as.data.frame(blackcap_mt_age_inter$summary[,-11]),
   caption = "Blackcap $M_t$ inter-winter model diagnostic and summary output",
   label = "tab:blackcap_mt_age_inter_summary_output", type = "latex"),
      file = "report/model_outputs/blackcap_mt_age_inter.tex")


chiffchaff_mt_age_inter <- jags_model_multinomial_age_inter(chiffchaff_inter,
                                        inits = inits_mt_age_inter,
                                        model_type = "mt",
                                        hyper_par = list(
                                        alpha_phi_j = alpha_phi_j_chiffchaff,
                                        beta_phi_j = beta_phi_j_chiffchaff,
                                        lpha_phi_ad= alpha_phi_ad_chiffchaff,
                                        beta_phi_ad= beta_phi_ad_chiffchaff,
                                        alpha_p = alpha_p,
                                        beta_p = beta_p),
                                        params = c("phi.juv", "phi.ad", "p",
                                                   "phi.diff", "fit.juv",
                                                   "fit.new.juv", "fit.ad",
                                                   "fit.new.ad"))

MCMCplot(chiffchaff_mt_age_inter, params = c("phi.diff"))
print(xtable(as.data.frame(chiffchaff_mt_age_inter$summary[,-11]),
 caption = "Chiffchaff $M_t$ inter-winter model diagnostic and summary output",
 label = "tab:chiffchaff_mt_age_inter_summary_output", type = "latex"),
      file = "report/model_outputs/chiffchaff_mt_age_inter.tex")

robin_mt_age_inter <- jags_model_multinomial_age_inter(robin_inter,
                                             inits = inits_mt_age_inter,
                                             model_type = "mt",
                                             hyper_par = list(
                                             alpha_phi_j = alpha_phi_j_robin,
                                             beta_phi_j = beta_phi_j_robin,
                                             alpha_phi_ad = alpha_phi_ad_robin,
                                             beta_phi_ad = beta_phi_ad_robin,
                                             alpha_p = alpha_p,
                                             beta_p = beta_p),
                                             params = c("phi.juv", "phi.ad",
                                                        "p", "phi.diff",
                                                        "fit.juv", "fit.new.juv",
                                                        "fit.ad", "fit.new.ad"))

robin_mt_age_inter
MCMCplot(robin_mt_age_inter, params = c("phi.diff"))
print(xtable(as.data.frame(robin_mt_age_inter$summary[,-11]),
       caption = "Robin $M_t$ inter-winter model diagnostic and summary output",
            label = "tab:robin_mt_age_inter_summary_output", type = "latex"),
      file = "report/model_outputs/robin_mt_age_inter.tex")


##############################
# INTRA-WINTER MODELLING
##############################

blackcap_m0_intra <- jags_model_multinomial_intra(blackcap_intra,
                                                  inits = inits_m0_intra,
                                                  model_type = "m0",
                                                  hyper_par = list(
                                                  alpha_phi = alpha_phi_j_robin,
                                                  beta_phi = beta_phi_j_robin,
                                                  alpha_p = alpha_p,
                                                  beta_p = beta_p),
                                                  params = c("mean.phi",
                                                             "mean.p",
                                                             "fit", "fit.new"))
print(xtable(as.data.frame(blackcap_m0_intra$summary[,-11]),
   caption = "Blackcap $M_t$ intra-winter model diagnostic and summary output",
   label = "tab:blackcap_m0_intra", type = "latex"),
      file = "report/model_outputs/blackcap_mt_intra.tex")

chiffchaff_m0_intra <- jags_model_multinomial_intra(chiffchaff_intra,
                                                  inits = inits_m0_intra,
                                                  model_type = "m0",
                                                  hyper_par = list(
                                                  alpha_phi = alpha_phi_j_robin,
                                                  beta_phi = beta_phi_j_robin,
                                                  alpha_p = alpha_p,
                                                  beta_p = beta_p),
                                                  params = c("mean.phi",
                                                             "mean.p",
                                                             "fit", "fit.new"))
print(xtable(as.data.frame(chiffchaff_m0_intra$summary[,-11]),
   caption = "Chiffchaff $M_t$ intra-winter model diagnostic and summary output",
   label = "tab:blackcap_m0_intra", type = "latex"),
      file = "report/model_outputs/chiffchaff_mt_intra.tex")

robin_m0_intra <- jags_model_multinomial_intra(robin_intra,
                                               inits = inits_m0_intra,
                                               model_type = "m0",
                                               hyper_par = list(
                                               alpha_phi = alpha_phi_j_robin,
                                               beta_phi = beta_phi_j_robin,
                                               alpha_p = alpha_p,
                                               beta_p = beta_p),
                                               params = c("mean.phi",
                                                          "mean.p",
                                                          "fit", "fit.new"))
print(xtable(as.data.frame(robin_m0_intra$summary[,-11]),
      caption = "Robin $M_t$ intra-winter model diagnostic and summary output",
             label = "tab:robin_m0_intra", type = "latex"),
      file = "report/model_outputs/robin_mt_intra.tex")

blackcap_mt_intra <- jags_model_multinomial_intra(blackcap_intra,
                                                  inits = inits_mt_intra,
                                                  model_type = "mt",
                                                  hyper_par = list(
                                                  alpha_phi = alpha_phi_j_robin,
                                                  beta_phi = beta_phi_j_robin,
                                                  alpha_p = alpha_p,
                                                  beta_p = beta_p),
                                                  params = c("phi",
                                                             "p",
                                                             "fit", "fit.new"))
print(xtable(as.data.frame(blackcap_mt_intra$summary[,-11]),
   caption = "Blackcap $M_t$ intra-winter model diagnostic and summary output",
   label = "tab:blackcap_m0_intra", type = "latex"),
      file = "report/model_outputs/blackcap_mt_intra.tex")

chiffchaff_mt_intra <- jags_model_multinomial_intra(chiffchaff_intra,
                                                inits = inits_mt_intra,
                                                model_type = "mt",
                                                hyper_par = list(
                                                alpha_phi = alpha_phi_j_robin,
                                                beta_phi = beta_phi_j_robin,
                                                alpha_p = alpha_p,
                                                beta_p = beta_p),
                                                params = c("phi", "p", "fit",
                                                           "fit.new"))

print(xtable(as.data.frame(chiffchaff_mt_intra$summary[,-11]),
 caption = "Chiffchaff $M_t$ intra-winter model diagnostic and summary output",
 label = "tab:chiffchaff_m0_intra", type = "latex"),
      file = "report/model_outputs/chiffchaff_mt_intra.tex")

robin_mt_intra <- jags_model_multinomial_intra(robin_intra,
                                               inits = inits_mt_intra,
                                               model_type = "mt",
                                               hyper_par = list(
                                               alpha_phi = alpha_phi_j_robin,
                                               beta_phi = beta_phi_j_robin,
                                               alpha_p = alpha_p,
                                               beta_p = beta_p),
                                               params = c("phi", "p", "fit",
                                                          "fit.new"))
print(xtable(as.data.frame(robin_mt_intra$summary[,-11]),
     caption = "Robin $M_t$ intra-winter model diagnostic and summary output",
     label = "tab:robin_m0_intra", type = "latex"),
      file = "report/model_outputs/robin_mt_intra.tex")

overall_CJS(as.matrix((birds_inter %>%
                         filter(age_at_ringing == "juvenile"))[, 4:14]),
            rep(1, nrow((birds_inter %>%
                           filter(age_at_ringing == "juvenile"))[, 4:14])))

overall_CJS(as.matrix(blackcap_inter[, 4:14]),
            rep(1, nrow(blackcap_inter[, 4:14])))
overall_CJS(as.matrix(chiffchaff_inter[, 4:14]),
            rep(1, nrow(chiffchaff_inter[, 4:14])))
overall_CJS(as.matrix(robin_inter[, 4:14]),
            rep(1, nrow(robin_inter[, 4:14])))

overall_CJS(as.matrix((blackcap_inter %>%
                         filter(age_at_ringing == "juvenile"))[, 4:14]),
            rep(1, nrow((blackcap_inter %>%
                           filter(age_at_ringing == "juvenile"))[, 4:14])))
overall_CJS(as.matrix((chiffchaff_inter %>%
                         filter(age_at_ringing == "juvenile"))[, 4:14]),
            rep(1, nrow((chiffchaff_inter %>%
                           filter(age_at_ringing == "juvenile"))[, 4:14])))
overall_CJS(as.matrix((robin_inter %>%
                         filter(age_at_ringing == "juvenile"))[, 4:14]),
            rep(1, nrow((robin_inter %>%
                           filter(age_at_ringing == "juvenile"))[, 4:14])))

overall_CJS(as.matrix((blackcap_inter %>%
                         filter(age_at_ringing == "adult"))[, 4:14]),
            rep(1, nrow((blackcap_inter %>%
                           filter(age_at_ringing == "adult"))[, 4:14])))
overall_CJS(as.matrix((chiffchaff_inter %>%
                         filter(age_at_ringing == "adult"))[, 4:14]),
            rep(1, nrow((chiffchaff_inter %>%
                           filter(age_at_ringing == "adult"))[, 4:14])))
overall_CJS(as.matrix((robin_inter %>%
                         filter(age_at_ringing == "adult"))[, 4:14]),
            rep(1, nrow((robin_inter %>%
                           filter(age_at_ringing == "adult"))[, 4:14])))


overall_CJS(as.matrix((birds_intra %>%
                         filter(age_at_ringing == "juvenile"))[, 6:12]),
            rep(1, nrow((birds_intra %>%
                           filter(age_at_ringing == "juvenile"))[, 6:12])))

overall_CJS(as.matrix(blackcap_intra[, 6:12]),
            rep(1, nrow(blackcap_intra[, 6:12])))
overall_CJS(as.matrix(chiffchaff_intra[, 6:12]),
            rep(1, nrow(chiffchaff_intra[, 6:12])))
overall_CJS(as.matrix(robin_intra[, 6:12]),
            rep(1, nrow(robin_intra[, 6:12])))



model_comparison <- data.frame(
  Species = c("Blackcap", "Chiffchaff", "Robin"),
  M_0 = c(blackcap_m0_age_inter$DIC, chiffchaff_m0_age_inter$DIC, robin_m0_age_inter$DIC),
  M_t = c(blackcap_mt_age_inter$DIC, chiffchaff_mt_age_inter$DIC, robin_mt_age_inter$DIC)
)

print(xtable(model_comparison, type = "latex", caption = "Model comparison",
             label = "tab:model_comparison"),
      file = "report/model_outputs/model_comparison.tex")


# Function to calculate p-value
calc_pvalue <- function(model, suffix = ""){
  return(sum(model$sims.list[[paste0("fit", suffix)]] >
               model$sims.list[[paste0("fit.new", suffix)]]) /
           length(model$sims.list[[paste0("fit.new", suffix)]]))
}

# Creating a dataframe of all p-values
p_values <- data.frame(
  species = c("Blackcap", "Chiffchaff", "Robin"),
  juveniles_m0 = c(calc_pvalue(blackcap_m0_age_inter, ".juv"),
                   calc_pvalue(chiffchaff_m0_age_inter, ".juv"),
                   calc_pvalue(robin_m0_age_inter, ".juv")),
  adults_m0 = c(calc_pvalue(blackcap_m0_age_inter, ".ad"),
                calc_pvalue(chiffchaff_m0_age_inter, ".ad"),
                calc_pvalue(robin_m0_age_inter, ".ad")),
  juveniles_mt = c(calc_pvalue(blackcap_mt_age_inter, ".juv"),
                   calc_pvalue(chiffchaff_mt_age_inter, ".juv"),
                   calc_pvalue(robin_mt_age_inter, ".juv")),
  adults_mt = c(calc_pvalue(blackcap_mt_age_inter, ".ad"),
                calc_pvalue(chiffchaff_mt_age_inter, ".ad"),
                calc_pvalue(robin_mt_age_inter, ".ad"))
)


# Creating a LaTeX table using xtable and saving it to a file
latex_pvalues <- xtable(p_values,
      caption="Inter-winter model goodness of fit results (bayesian p-values)",
      label="tab:gof_results_inter")
print(latex_pvalues, type="latex", file="report/model_outputs/gof_results_age_inter.tex")


model_comparison_intra <- data.frame(
  Species = c("Blackcap", "Chiffchaff", "Robin"),
  M_0 = c(blackcap_m0_intra$DIC, chiffchaff_m0_intra$DIC,
          robin_m0_intra$DIC),
  M_t = c(blackcap_mt_intra$DIC, chiffchaff_mt_intra$DIC,
          robin_mt_intra$DIC)
)

print(xtable(model_comparison_intra, type = "latex",
             caption = "Intra-winter model comparisons",
             label = "tab:model_comparison_intra"),
      file = "report/model_outputs/model_comparison_intra.tex")

# Creating a dataframe of all p-values
p_values_intra <- data.frame(
  species = c("Blackcap", "Chiffchaff", "Robin"),

  m0 = c(calc_pvalue(blackcap_m0_intra, ""),
         calc_pvalue(chiffchaff_m0_intra, ""),
         calc_pvalue(robin_m0_intra, "")),

  mt = c(calc_pvalue(blackcap_mt_intra, ""),
         calc_pvalue(chiffchaff_mt_intra, ""),
         calc_pvalue(robin_mt_intra, ""))
)

latex_pvalues_intra <- xtable(p_values_intra,
    caption="Intra-winter model goodness of fit results (bayesian p-values)",
                              label="tab:gof_results_intra")
print(latex_pvalues_intra, type="latex",
      file="report/model_outputs/gof_results_intra.tex")


species_list <- list(blackcap_mt_age_inter, chiffchaff_mt_age_inter, robin_mt_age_inter)
species_names <- c("Blackcap", "Chiffchaff", "Robin")

# Path to save the plots
path_to_save <- "report/images/"

for (i in 1:3) {
  leis.result <- species_list[[i]]
  T <- length(leis.result$mean$phi.ad)
  lower_ad <- upper_ad <- lower_juv <- upper_juv <- lower_diff <- upper_diff <- numeric(T)
  for (t in 1:T) {
    lower_ad[t] <- quantile(leis.result$sims.list$phi.ad[,t], 0.025)
    upper_ad[t] <- quantile(leis.result$sims.list$phi.ad[,t], 0.975)
    lower_juv[t] <- quantile(leis.result$sims.list$phi.juv[,t], 0.025)
    upper_juv[t] <- quantile(leis.result$sims.list$phi.juv[,t], 0.975)
    lower_diff[t] <- quantile(leis.result$sims.list$phi.diff[,t], 0.025)
    upper_diff[t] <- quantile(leis.result$sims.list$phi.diff[,t], 0.975)
  }

  years <- seq(2007, 2016)

  for (j in c("ad", "juv", "diff")) {
    y_values <- leis.result$mean[[paste0("phi.", j)]]
    lower_values <- get(paste0("lower_", j))
    upper_values <- get(paste0("upper_", j))

    df <- data.frame(Year = years, Value = y_values, Lower = lower_values,
                     Upper = upper_values)

    p <- ggplot(df, aes(x = Year, y = Value)) +
      geom_line(color = "blue") +
      geom_point(shape = 16, color = "blue") +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "red", alpha = 0.2) +
      ylim(range(c(lower_values, upper_values))) +
      labs(y = "Probability", x = "") +
      custom_theme +
      theme(axis.text.x = element_text(angle = 45, margin = margin(t = 15)),
            panel.grid.major.x = element_line(color = "grey"),
            panel.grid.minor.x = element_blank()) +
      scale_x_continuous(breaks = years) # Optional, to explicitly set the x breaks

    print(p)
    ggsave(filename = paste0(path_to_save, species_names[i], "_", j,
                             "_inter.png"), plot = p, width = 6, height = 4)
  }
}

species_list_intra <- list(blackcap_mt_intra, chiffchaff_mt_intra,
                           robin_mt_intra)

for (i in 1:3) {
  leis.result <- species_list_intra[[i]]
  T <- length(leis.result$mean$phi)
  lower_phi <- upper_phi <- numeric(T)
  for (t in 1:T) {
    lower_phi[t] <- quantile(leis.result$sims.list$phi[,t], 0.025)
    upper_phi[t] <- quantile(leis.result$sims.list$phi[,t], 0.975)
  }

  y_values <- leis.result$mean$phi

  months <- c("October", "November", "December", "January", "February", "March")
  df <- data.frame(Month = months, Value = y_values, Lower = lower_phi,
                   Upper = upper_phi)
  df$Month <- factor(df$Month, levels = months)


  p <- ggplot(df, aes(x = Month, y = Value, group = 1)) +
    geom_line(color = "blue") +
    geom_point(shape = 16, color = "blue") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "red", alpha = 0.2) +
    ylim(range(c(lower_phi, upper_phi))) +
    labs(y = "Probability", x = "") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, margin = margin(t = 15)),
          panel.grid.major.x = element_line(color = "grey"),
          panel.grid.minor.x = element_blank())



  print(p)
  ggsave(filename = paste0(path_to_save, species_names[i], "_", "intra.png"),
         plot = p, width = 6, height = 4)
}







