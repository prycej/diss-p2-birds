#------------------------------------------------------------------------------
# Script Name:        main.R
# Description:        r script to perform the exploratory data analysis,
#                     modelling for predicting inter-winter and intra-winter
#                     survival given capture history data
#
# Author:             Jacob Pryce
# Date:               2023-08-06
# Version:            1.0
# Dependencies: - data.table: Enhanced, efficient data manipulation
#               - stats4: Basic statistical functions, estimation
#               - dplyr: Grammar of data manipulation
#               - tidyr: Data cleaning, organisation
#               - stringr: Easy string manipulation
#               - ggplot2: Declarative graphics creation system
#               - jagsUI: Interface for JAGS Bayesian models
#               - xtable: Exporting tables to LaTeX
#
# Data Sources: - "CR.blackcap_FixRing.csv"
#               - "CR.chifchaf_FixRing.csv"
#               - "CR.robin_FixRing.csv"
#------------------------------------------------------------------------------

library(data.table) # Enhanced, efficient data manipulation
library(stats4) # Basic statistical functions, estimation
library(dplyr) # Grammar of data manipulation
library(tidyr) # Data cleaning, organisation
library(stringr) # Easy string manipulation
library(ggplot2) # Declarative graphics creation system
library(jagsUI)  # Interface for JAGS Bayesian models
library(xtable) # Exporting tables to LaTeX

# Check and create directory for images
if (!dir.exists("report/images")) {
  dir.create("report/images", recursive = TRUE)
}

# Check and create directory for model outputs
if (!dir.exists("report/model_outputs")) {
  dir.create("report/model_outputs", recursive = TRUE)
}


# Define file paths
file_paths <- c("CR.blackcap_FixRing.csv",
                "CR.chifchaf_FixRing.csv",
                "CR.robin_FixRing.csv")

# Check if files exist and read them in if they do
if(all(file.exists(file_paths))){
  blackcap_data_raw <- tryCatch({
    fread("CR.blackcap_FixRing.csv", header = TRUE)
  }, error = function(e) {
    message("Error reading CR.blackcap_FixRing.csv: ", e)
    return(NULL)
  })
  
  chiffchaff_data_raw <- tryCatch({
    fread("CR.chifchaf_FixRing.csv", header = TRUE)
  }, error = function(e) {
    message("Error reading CR.chifchaf_FixRing.csv: ", e)
    return(NULL)
  })
  
  robin_data_raw <- tryCatch({
    fread("CR.robin_FixRing.csv", header = TRUE)
  }, error = function(e) {
    message("Error reading CR.robin_FixRing.csv: ", e)
    return(NULL)
  })
} else {
  stop("One or more files do not exist. Please check the file paths.")
}


##############################
# DATA PREPARATION
##############################

clean_dataset <- function(df, species) {
  #' Clean a dataset by adding an ID column, removing unwanted columns, and
  #' filtering rows
  #'
  #' @param df Dataframe containing the raw data.
  #' @param species Character string indicating the species name.
  #' @return A cleaned dataframe with a unique ID for each row, unnecessary
  #'  columns removed, and rows with no capture information filtered out.
  #'
  #'
  if (!is.data.frame(df)) stop("df is not a dataframe.")
  if (!is.character(species) | length(species) != 1) {
    stop("species must be asingle character string.")
  }
  
  # Add row number and species name columns
  df <- df %>% mutate(rn = row_number(), species = species) %>%
    # Create a unique ID for each row using species name and row number, then
    # remove the temporary rn column
    mutate(id = paste0(species, "_", rn)) %>% select(-rn)
  
  # Remove the first four months (as they represent only a partial winter)
  df <- df[, 5:(ncol(df))]
  
  # Filter out rows with no capture information, by ensuring the sum of values
  #  (excluding certain columns) is greater than 0
  df <- df %>% filter(rowSums(select(., -c(id, age_at_ringing, species))) > 0)
  
  # Return the cleaned dataset
  df
}

blackap_data_cleaned <- clean_dataset(blackcap_data_raw, 'blackcap')
chiffchaff_data_cleaned <- clean_dataset(chiffchaff_data_raw, 'chiffchaff')
robin_data_cleaned <- clean_dataset(robin_data_raw, 'robin')

birds <- rbind(blackap_data_cleaned, chiffchaff_data_cleaned,
               robin_data_cleaned)

##############################
# EDA
##############################
size <- 40
# custom theme for plotting purposes
custom_theme <- theme(
  axis.text.x = element_text(size = size),
  axis.text.y = element_text(size = size),
  axis.title.x = element_text(size = size, margin = margin(t = 10)),
  axis.title.y = element_text(size = size, margin = margin(r = 10)),
  plot.title = element_text(size = size, hjust = 0.5)
)

###### Age at first ringing proportions by species #####
# Order the levels of the 'species' variable
birds$species <- factor(birds$species, levels = c("robin", "chiffchaff",
                                                  "blackcap"))

ggplot(birds, aes(fill=age_at_ringing, x=species)) +
  geom_bar(position="fill", width=0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Unknown" = "grey", "adult" = "dark red",
                               "juvenile" = "pink")) +
  labs(x = "Species", y = "Proportion", fill = "Age at first ringing") +
  custom_theme +
  theme(legend.text = element_text(size = size / (4/3)),
        legend.title = element_text(size = size / (4/3)),  
        legend.position = "bottom",
        axis.text.x = element_text(size = size / (4/3)),
        axis.title.x = element_text(size = size / (4/3)))
ggsave("report/images/age_props.png")


#### Capture count proportions by species ####
birds$row_sum <- rowSums(select_if(birds, is.numeric))
# Convert the row_sum to a factor
birds$row_sum <- as.factor(birds$row_sum)
# Plot the data
ggplot(birds, aes(fill=row_sum, x=species)) +
  geom_bar(position="fill", width = 0.7) +
  scale_fill_discrete(name = "Capture count") +
  coord_flip() +
  labs(x = "Species", y = "Proportion") +
  custom_theme +
  theme(legend.text = element_text(size = size/ (4/3)),
        legend.title = element_text(size = size / (4/3)),  
        legend.position = "bottom",
        axis.text.x = element_text(size = size / (4/3)),
        axis.title.x = element_text(size = size / (4/3)))
ggsave("report/images/capture_count_proportions.png")


birds_long <- birds %>%
  # Pivot the dataset from wide to long, taking columns that start with "20"
  # as date information and storing the corresponding values in "count"
  pivot_longer(cols = starts_with("20"), names_to = "date",
               values_to = "count") %>%
  # Extract year and month from the "date" column as separate numerical variables
  mutate(
    year = as.numeric(substr(date, 1, 4)),
    month = as.numeric(substr(date, 5, 6)),
    # Determine the winter period based on the month, considering months from
    #  October as belonging to the next winter period
    winter_period = ifelse(month >= 10, year, year - 1)
  ) %>%
  # Transform month number into winter month names, starting from October as
  # month 1, November as month 2, etc.
  mutate(
    month_name = factor(ifelse(month >= 10, month - 9, month + 3),
                        levels = c(1:7),
                        labels = c("Oct", "Nov", "Dec",
                                   "Jan", "Feb", "Mar", "Apr"))
  )


###### Capture count by winter period ####
for (species_idx in c("blackcap", "chiffchaff", "robin")) {
  winter_sums <- birds_long %>%
    filter(species == species_idx) %>%
    group_by(winter_period, age_at_ringing) %>%
    summarise(count = sum(count, na.rm = TRUE))
  
  p_winter <- ggplot(winter_sums, aes(x = factor(winter_period), y = count,
                                      fill = age_at_ringing)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(name = "Winter Period", 
                     labels = function(x) paste0("W", as.numeric(x) - 1)) +
    labs(y = "Count", fill = "Age") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom", legend.text = element_text(size = size),
          legend.title = element_text(size = size / (4/3)))
  
  print(p_winter)
  ggsave(paste0("report/images/", species_idx, "_capture_counts_per_year.png"),
         p_winter, width = 10, height = 8)
}

##### Capture count across months by age at ringing and for each species ####
for (species_idx in c("blackcap", "chiffchaff", "robin")) {
  monthly_sums <- birds_long %>%
    filter(species == species_idx) %>%
    group_by(month_name, age_at_ringing) %>%
    summarise(count = sum(count, na.rm = TRUE))
  
  p <- ggplot(monthly_sums, aes(x = month_name, y = count,
                                fill = age_at_ringing)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Month", y = "Count", fill = "Age") +
    scale_x_discrete(limit = c("Oct", "Nov", "Dec", "Jan",
                               "Feb", "Mar", "Apr")) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom", legend.text = element_text(size = size),
          legend.title = element_text(size = size / (4/3)))
  print(p)
  ggsave(paste0("report/images/", species_idx, "_capture_count_per_month.png"),
         plot = p, width = 10, height = 8)
}

##############################
# MODELLING SETUP
##############################

create_inter_winter_dataset <- function(df) {
  #' Create Inter-Winter Dataset
  #'
  #' This function takes a dataset containing bird observation data and
  #' collapses it by creating a new column 'winter_presence', which calculates
  #' if a bird was captured at all during a particular winter season.
  #'
  #' The dataset is restructured by identifying the years and months of
  #' observation and grouping by winter season (October to April)
  #' For each winter season, a presence value is computed for each bird,
  #' indicating whether the bird was captured at least once during that winter.
  #'
  #' @param df A data frame containing bird observation data. It must include
  #' columns for 'id', 'age_at_ringing', 'species', and 'row_sum', plus
  #' additional columns  representing the months of observation.
  #' @return A data frame grouped by bird 'id', 'species', and 'age_at_ringing',
  #'  with columns for each winter season, representing the presence (1) or
  #'  absence (0)  of the bird during that winter. The resulting data frame is
  #'  sorted by species and numeric ID.
  
  # Check if df is a data.frame and has necessary columns
  if (!is.data.frame(df)) stop("df must be a dataframe.")
  if (!all(c('id', 'age_at_ringing', 'species', 'row_sum') %in% colnames(df)))
    stop("df must have 'id', 'age_at_ringing', 'species', and 'row_sum' columns.")
  
  df %>%
    pivot_longer(-c(id, age_at_ringing, species, row_sum),
                 names_to = "month", values_to = "value") %>% # Pivot to longer format
    mutate(year = as.integer(str_sub(month, start = 1, end = 4)),
           month_num = as.integer(str_sub(month, start = 5, end = 6))) %>% # Extract
    # year & month
    mutate(winter = ifelse(month_num >= 10, year, year - 1)) %>% # Calculate
    # winter period
    group_by(id, winter, age_at_ringing, species) %>%
    summarise(winter_presence = ifelse(any(value == 1), 1, 0)) %>% # Calculate
    # if a bird appeared in that winter
    pivot_wider(names_from = winter, values_from = winter_presence) %>%
    mutate(numeric_id = as.numeric(str_extract(id, "\\d+"))) %>% # Extract
    # numeric ID
    arrange(species, numeric_id) %>% # Arrange by species & numeric ID
    select(-numeric_id) %>%
    ungroup
  
}

birds_inter <- create_inter_winter_dataset(birds)

create_intra_winter_dataset <- function(df) {
  #' Create Intra-Winter Dataset
  #'
  #' This function transforms a given dataset into an intra-winter format,
  #' pivoting and reformatting the data so that each winter has its own column.
  #' The dataset is expected to contain bird observation data, with columns
  #' representing years and months.
  #'
  #' @param df Data frame with bird observation data, with columns starting
  #' with "20" representing year and month.
  #' @return A data frame where each row represents an observation, with columns
  #'   for winter months (Oct, Nov, Dec, Jan, Feb, Mar, Apr) and additional
  #'   information like species and ID. If there are duplicated IDs due to
  #'   multiple captures in different years, only the first year is considered.
  #' @note This function assumes that the input data frame contains specific
  #'   columns, including a unique identifier 'id' and a species identifier.
  #'   Any rows where the values for months 1 to 7 are all zero will be filtered
  #'   out.
  
  # Check if df is a data.frame and has necessary columns
  if (!is.data.frame(df)) stop("df must be a dataframe.")
  if (!"id" %in% colnames(df)) stop("df must have 'id' column.")
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
    group_by(id) %>%
    slice_min(row_number(), n = 1) %>% # Where the id column is duplicated
    # because of a capture appearing multiple times for different years, take
    # the first year
    ungroup() %>%
    # sort by id column
    mutate(numeric_id = as.numeric(str_extract(id, "\\d+"))) %>%
    arrange(species, numeric_id) %>%
    select(-numeric_id)
  
  # Renaming columns
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

######### Separate out by species #####
blackcap_inter <- birds_inter %>% filter(species == 'blackcap') %>% ungroup
robin_inter <- birds_inter %>% filter(species == 'robin') %>% ungroup
chiffchaff_inter <- birds_inter %>% filter(species == 'chiffchaff') %>% ungroup


blackcap_intra <- birds_intra %>% filter(species == 'blackcap') %>% ungroup
robin_intra <- birds_intra %>% filter(species == 'robin') %>% ungroup
chiffchaff_intra <- birds_intra %>% filter(species == 'chiffchaff') %>% ungroup
########
set.seed(42)

n_iter <-10000 # Number of iterations
n_burnin <-2000 # Number of burnin samples
n_chains <-3
n_adapt <- 1000

marray <- function(CH){
  #' Create an m-array from capture history (CH) data
  #'
  #' This function takes a capture history (CH) matrix, where rows represent
  #' individual animals and columns represent capture occasions. It returns an
  #' m-array, which is a matrix representing the transition of marked individuals
  #' between capture occasions, including those not recaptured.
  #'
  #' @param CH A matrix of capture history (CH). Rows correspond to individuals,
  #'   and columns correspond to occasions. The value is 1 if captured, 0
  #'   otherwise.
  #' @return A matrix representing the m-array, with columns for each occasion
  #'   and an additional column for never recaptured individuals.
  #'
  n_ind <- dim(CH)[1] # Number of individuals
  n_occasions <- dim(CH)[2] # Number of occasions (time periods)
  
  # Initialise the m-array with zeros, with columns for each occasion and an
  # additional one for never recaptured
  m_array <- matrix(data = 0, ncol = n_occasions+1, nrow = n_occasions)
  
  # Iterate through each occasion, summing the captures for that time period
  for (t in 1:n_occasions){
    m_array[t,1] <- sum(CH[,t])
  }
  
  # Iterate through each individual
  for (i in 1:n_ind){
    pos <- which(CH[i,] != 0) # Find positions where captures occurred
    num_captures <- length(pos) # Count the number of captures for this
    # individual
    
    # Iterate through captures, incrementing the counts in the corresponding
    # m-array cells
    for (j in 1:(num_captures-1)){
      m_array[pos[j], pos[j+1]] <- m_array[pos[j], pos[j+1]] + 1
    }
  }
  
  # Iterate through each occasion again, calculating the number of individuals
  # that were never recaptured
  for (t in 1:n_occasions){
    m_array[t,n_occasions+1] <- m_array[t,1] - sum(m_array[t,2:n_occasions])
  }
  
  # Return the m-array, excluding the last row and first column as they contain
  # redundant information
  res <- m_array[1:(n_occasions-1),2:(n_occasions+1)]
  return(res)
}

## Choose random initial values
inits_m0_inter <- function(chain){
  list(
    phi_juv0 = runif(1),
    phi_adult0 = runif(1),
    p0 = runif(1))
}

inits_mt_age_inter <- function(chain){
  list(
    phi_juv = runif(10, 0, 1),
    phi_adult = runif(10, 0, 1),
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
                                             model_type, params) {
  #' Run inter-winter JAGS model
  #'
  #' This function processes capture history data, separates juveniles and adults,
  #' and runs a JAGS model to estimate survival and recapture probabilities for
  #' different age classes. The JAGS model is defined in an external file.
  #'
  #' @param df A data frame containing capture history data, with columns for
  #'   individual ID, age at ringing (juvenile or adult), species, and capture
  #'   occasions.
  #' @param hyper_par A list containing hyperparameters for the JAGS model.
  #' @param inits A list containing initial values for the JAGS model parameters.
  #' @param model_type A character string indicating the model type
  #' @param params A character vector containing the names of parameters to be
  #'  returned from the JAGS model
  #' @return An object containing the JAGS model output.
  #'
  
  
  # Check if df is a data.frame and has necessary columns
  if (!is.data.frame(df)) stop("df must be a dataframe.")
  if (!all(c("id", "age_at_ringing", "species") %in% colnames(df))) {
    stop("df must have 'id', 'age_at_ringing' and 'species' columns.")
  }
  
  # Check if hyper_par and inits are lists
  if (!is.list(hyper_par)) stop("hyper_par must be a list.")
  if (!is.function(inits)) stop("inits must be a function")
  
  # Check if model_type is a character string
  if (!is.character(model_type)) stop("model_type must be a character string.")
  
  # Check if params is a character vector
  if (!is.character(params)) stop("params must be a character vector.")
  
  CH_juv <- df %>% filter(age_at_ringing == "juvenile") %>%
    select(c(-id, -age_at_ringing, -species)) # Juvenile capture history
  CH_adult <- df %>% filter(age_at_ringing == "adult") %>%
    select(c(-id, -age_at_ringing, -species)) # Adult capture history
  
  # Calculate the sum of captures for each juvenile
  cap_sum <- apply(CH_juv, 1, sum)
  
  # Identify juveniles that have been recaptured at least once
  idx <- which(cap_sum >= 2)
  
  # Separate juveniles into two categories: those recaptured at least once,
  # and those never recaptured
  CH_juv_recap <- CH_juv[idx,]  # Juveniles recaptured at least once
  CH_juv_never_recap <- CH_juv[-idx,] # Juveniles never recaptured
  
  # Process juveniles that have been recaptured at least once to remove the
  # first capture
  first <- numeric()
  for (i in 1:dim(CH_juv_recap)[1]) {
    first[i] <- min(which(CH_juv_recap[i,] == 1)) # Identify the first recapture
    # for each juvenile
    CH_juv_recap[i, first[i]] <- 0            # Set the first recapture to 0
  }
  
  # Combine grown-up juveniles with adults and convert to m-array format
  CH_adult_m <- rbind(CH_adult, CH_juv_recap)
  marr_adult <- marray(CH_adult_m)
  
  # Create a matrix to store the first and second recaptures for juveniles
  CH_juv_recap2 <- matrix(0, nrow = dim(CH_juv_recap)[1],
                          ncol = dim(CH_juv_recap)[2])
  second <- numeric()
  for (i in 1:dim(CH_juv_recap)[1]) {
    second[i] <- min(which(CH_juv_recap[i,] == 1)) # Identify the second recapture
    # for each juvenile
    CH_juv_recap2[i, first[i]] <- 1         # Mark the first recapture
    CH_juv_recap2[i, second[i]] <- 1         # Mark the second recapture
  }
  
  # Convert the processed juvenile recapture data into m-array format
  CH_juv_recap_marray <- marray(CH_juv_recap2)
  CH_juv_recap_marray[, dim(CH_juv)[2]] <- 0 # Ensure that the last column is zero
  # since all of them are released as adults
  
  # Convert the juveniles that were never recaptured into m-array format
  CH_juv_never_recap_marray <- marray(CH_juv_never_recap)
  
  # Combine the m-arrays for recaptured and never-recaptured juveniles
  marr_juv <- CH_juv_recap_marray + CH_juv_never_recap_marray
  
  # Calculate the sum of releases for each capture occasion
  releases_juv <- rowSums(marr_juv)
  releases_adult <- rowSums(marr_adult)
  
  # prepare JAGS data
  jdata <- list(marr_juv = marr_juv, marr_adult = marr_adult,
                releases_juv = releases_juv, releases_adult = releases_adult,
                n_occasions=ncol(marr_juv),
                alpha_phi_juv = hyper_par$alpha_phi_juv,
                beta_phi_juv = hyper_par$beta_phi_juv,
                alpha_phi_adult = hyper_par$alpha_phi_adult,
                beta_phi_adult = hyper_par$beta_phi_adult,
                alpha_p = hyper_par$alpha_p,
                beta_p = hyper_par$beta_p)
  
  # run JAGS model
  out <- jags(jdata, inits, params,
              paste0("multinomial_",model_type,"_age_inter.jags"),
              DIC=TRUE, n.chains=n_chains, n.adapt=n_adapt, n.iter=n_iter,
              parallel = TRUE, seed = 42)
  out
}

jags_model_multinomial_intra <- function(df, hyper_par, inits,
                                         model_type, params) {
  #' Run intra-winter JAGS model
  #'
  #' This function prepares capture history data for intra-seasonal analysis,
  #' converting it into a suitable format and running a JAGS model to estimate
  #' survival and recapture probabilities. The JAGS model is defined in an
  #' external file.
  #'
  #' @param df A data frame containing capture history data, with columns for
  #'   individual ID, age at ringing, species, winter, and capture occasions.
  #' @param hyper_par A list containing hyperparameters for the JAGS model.
  #' @param inits A list containing initial values for the JAGS model parameters.
  #' @param model_type A character string indicating the model type
  #' @param params A character vector containing the names of parameters to be
  #'  returned from the JAGS model
  #'   by JAGS
  #' @return An object containing the JAGS model output.
  #'
  
  # Check if df is a data.frame and has necessary columns
  if (!is.data.frame(df)) stop("df must be a dataframe.")
  if (!all(c("id", "age_at_ringing", "species", "winter")
           %in% colnames(df))) {
    stop("df must have 'id', 'age_at_ringing', 'species' and 'winter' columns.")
  }
  
  # Check if hyper_par and inits are lists
  if (!is.list(hyper_par)) stop("hyper_par must be a list.")
  if (!is.function(inits)) stop("inits must be a function")
  
  # Check if model_type is a character string
  if (!is.character(model_type)) stop("model_type must be a character string.")
  
  # Check if params is a character vector
  if (!is.character(params)) stop("params must be a character vector.")
  
  # Create m-array
  marr <-  marray(df %>% select(-c(id, age_at_ringing, species, winter, 
                                   row_sum)))
  
  # Calculate the sum of releases for each capture occasion
  releases <- rowSums(marr)
  
  # prepare JAGS data
  jdata <- list(marr=marr,
                R = releases,
                n_occasions=ncol(marr),
                alpha_phi= hyper_par$alpha_phi,
                beta_phi = hyper_par$beta_phi,
                alpha_p = hyper_par$alpha_p,
                beta_p = hyper_par$beta_p)
  
  # run JAGS model
  out <- jags(jdata, inits, params,
              paste0("multinomial_",model_type,"_intra.jags"),
              DIC=TRUE, n.chains=n_chains, n.adapt=n_adapt, n.iter=n_iter,
              parallel = TRUE, seed = 42)
  out
}


##############################
# PRIOR DISTRIBUTIONS
##############################
# Define parameters for the Beta distributions
alpha_phi_juv_blackcap <- alpha_phi_adult_blackcap <- 4
beta_phi_juv_blackcap <- beta_phi_adult_blackcap <- 5.69

alpha_phi_juv_chiffchaff <- alpha_phi_adult_chiffchaff <- 2.24
beta_phi_juv_chiffchaff  <- beta_phi_adult_chiffchaff <- 3.9

alpha_phi_juv_robin <- alpha_phi_adult_robin <- 3.4
beta_phi_juv_robin <- beta_phi_adult_robin <- 2.41

alpha_p <- 1.5
beta_p <- 3.5

params <- list(c(alpha_phi_juv_blackcap, beta_phi_juv_blackcap),
               c(alpha_phi_juv_chiffchaff, beta_phi_juv_chiffchaff),
               c(alpha_phi_juv_robin, beta_phi_juv_robin),
               c(alpha_p, beta_p))
names <- c("blackcaps", "chiffchaffs", "robins", "capture_probability")
colors <- c("#8B0000", "#00008B", "#006400", "black")

# Generate and plot prior distributions
for(i in 1:4){
  x <- seq(0, 1, length.out = 100)
  y <- dbeta(x, params[[i]][1], params[[i]][2])
  df <- data.frame(x = x, y = y)
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(color = colors[i], linewidth = 3) +
    theme_minimal() +
    custom_theme +
    
    xlab("Probability") +
    ylab("Density")
  
  ggsave(filename = paste0("report/images/prior_distribution_for_", names[i],
                           ".png"))
  print(p)
}


##############################
# INTER-WINTER MODELLING
##############################
########## M_0 blackcap model ############
blackcap_m0_age_inter <- jags_model_multinomial_age_inter(blackcap_inter,
                                      inits = inits_m0_inter,
                                      model_type = "m0",
                                      hyper_par = list(
                                      alpha_phi_juv = alpha_phi_juv_blackcap,
                                      beta_phi_juv = beta_phi_juv_blackcap,
                                      alpha_phi_adult= alpha_phi_adult_blackcap,
                                      beta_phi_adult= beta_phi_adult_blackcap,
                                      alpha_p = alpha_p,
                                      beta_p = beta_p),
                                      params = c("phi_juv0",
                                                 "phi_adult0",
                                                 "p0", "mean_phi_diff",
                                                 "disc_obs_juv",
                                                 "disc_rep_juv",
                                                 "disc_obs_adult",
                                                 "disc_rep_adult"))

blackcap_m0_age_inter

print(xtable(as.data.frame(blackcap_m0_age_inter$summary[,-11]),
   caption = "Blackcap $M_0$ inter-winter model diagnostic and summary output",
   label = "tab:blackcap_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/blackcap_m0_age_inter.tex")

########## M_0 chiffchaff model ############
chiffchaff_m0_age_inter <- jags_model_multinomial_age_inter(chiffchaff_inter,
                                    inits = inits_m0_inter,
                                    model_type = "m0",
                                    hyper_par = list(
                                    alpha_phi_juv = alpha_phi_juv_chiffchaff,
                                    beta_phi_juv = beta_phi_juv_chiffchaff,
                                    alpha_phi_adult_= alpha_phi_adult_chiffchaff,
                                    beta_phi_adult = beta_phi_adult_chiffchaff,
                                    alpha_p = alpha_p,
                                    beta_p = beta_p),
                                    params = c("phi_juv0", "phi_adult0",
                                               "p0", "mean_phi_diff",
                                               "disc_obs_juv",
                                               "disc_rep_juv",
                                               "disc_obs_adult",
                                               "disc_rep_adult"))

chiffchaff_m0_age_inter

print(xtable(as.data.frame(chiffchaff_m0_age_inter$summary[,-11]),
 caption = "Chiffchaff $M_0$ inter-winter model diagnostic and summary output",
 label = "tab:chiffchaff_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/chiffchaff_m0_age_inter.tex")

########## M_0 robin model ############
robin_m0_age_inter <- jags_model_multinomial_age_inter(robin_inter,
                                         inits = inits_m0_inter,
                                         model_type = "m0",
                                         hyper_par = list(
                                         alpha_phi_juv = alpha_phi_juv_robin,
                                         beta_phi_juv = beta_phi_juv_robin,
                                         alpha_phi_adult= alpha_phi_adult_robin,
                                         beta_phi_adult= beta_phi_adult_robin,
                                         alpha_p = alpha_p,
                                         beta_p = beta_p),
                                         params = c("phi_juv0",
                                                    "phi_adult0",
                                                    "p0", "mean_phi_diff",
                                                    "disc_obs_juv",
                                                    "disc_rep_juv",
                                                    "disc_obs_adult",
                                                    "disc_rep_adult"))

robin_m0_age_inter
print(xtable(as.data.frame(robin_m0_age_inter$summary[,-11]),
       caption = "Robin $M_0$ inter-winter model diagnostic and summary output",
       label = "tab:robin_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/robin_m0_age_inter.tex")

########## M_t blackcap model ############
blackcap_mt_age_inter <- jags_model_multinomial_age_inter(blackcap_inter,
                                      inits = inits_mt_age_inter,
                                      model_type = "mt",
                                      hyper_par = list(
                                      alpha_phi_juv = alpha_phi_juv_blackcap,
                                      beta_phi_juv = beta_phi_juv_blackcap,
                                      alpha_phi_adult= alpha_phi_adult_blackcap,
                                      beta_phi_adult= beta_phi_adult_blackcap,
                                      alpha_p = alpha_p,
                                      beta_p = beta_p),
                                      params = c("phi_juv", "phi_adult",
                                                 "p", "phi_diff",
                                                 "disc_obs_juv",
                                                 "disc_rep_juv",
                                                 "disc_obs_adult",
                                                 "disc_rep_adult"))
blackcap_mt_age_inter

print(xtable(as.data.frame(blackcap_mt_age_inter$summary[,-11]),
   caption = "Blackcap $M_t$ inter-winter model diagnostic and summary output",
   label = "tab:blackcap_mt_age_inter_summary_output", type = "latex"),
      file = "report/model_outputs/blackcap_mt_age_inter.tex")

########## M_t chiffchaff model ############
chiffchaff_mt_age_inter <- jags_model_multinomial_age_inter(chiffchaff_inter,
                                  inits = inits_mt_age_inter,
                                  model_type = "mt",
                                  hyper_par = list(
                                  alpha_phi_juv = alpha_phi_juv_chiffchaff,
                                  beta_phi_juv = beta_phi_juv_chiffchaff,
                                  alpha_phi_adult = alpha_phi_adult_chiffchaff,
                                  beta_phi_adult = beta_phi_adult_chiffchaff,
                                  alpha_p = alpha_p,
                                  beta_p = beta_p),
                                  params = c("phi_juv", "phi_adult", "p",
                                             "phi_diff",
                                             "disc_obs_juv",
                                             "disc_rep_juv",
                                             "disc_obs_adult",
                                             "disc_rep_adult"))

chiffchaff_mt_age_inter

print(xtable(as.data.frame(chiffchaff_mt_age_inter$summary[,-11]),
 caption = "Chiffchaff $M_t$ inter-winter model diagnostic and summary output",
 label = "tab:chiffchaff_mt_age_inter_summary_output", type = "latex"),
      file = "report/model_outputs/chiffchaff_mt_age_inter.tex")

########## M_t robin model ############
robin_mt_age_inter <- jags_model_multinomial_age_inter(robin_inter,
                                       inits = inits_mt_age_inter,
                                       model_type = "mt",
                                       hyper_par = list(
                                       alpha_phi_juv = alpha_phi_juv_robin,
                                       beta_phi_juv = beta_phi_juv_robin,
                                       alpha_phi_adult = alpha_phi_adult_robin,
                                       beta_phi_adult = beta_phi_adult_robin,
                                       alpha_p = alpha_p,
                                       beta_p = beta_p),
                                       params = c("phi_juv", "phi_adult",
                                                  "p", "phi_diff",
                                                  "disc_obs_juv",
                                                  "disc_rep_juv",
                                                  "disc_obs_adult",
                                                  "disc_rep_adult"))

robin_mt_age_inter

print(xtable(as.data.frame(robin_mt_age_inter$summary[,-11]),
     caption = "Robin $M_t$ inter-winter model diagnostic and summary output",
     label = "tab:robin_mt_age_inter_summary_output", type = "latex"),
      file = "report/model_outputs/robin_mt_age_inter.tex")

##############################
# INTER-MODEL SENSITIVITY ANALYSIS
##############################

sensitivity_analysis_inter_model <- jags_model_multinomial_age_inter(robin_inter,
                                     inits = inits_mt_age_inter,
                                     model_type = "mt",
                                     hyper_par = list(
                                     alpha_phi_juv = alpha_phi_juv_robin + 1,
                                     beta_phi_juv = beta_phi_juv_robin + 1,
                                     alpha_phi_adult= alpha_phi_adult_robin + 1,
                                     beta_phi_adult= beta_phi_adult_robin + 1,
                                     alpha_p = alpha_p + 1,
                                     beta_p = beta_p + 1),
                                     params = c("phi_juv", "phi_adult", "p",
                                                "phi_diff",
                                                "disc_obs_juv",
                                                "disc_rep_juv",
                                                "disc_obs_adult",
                                                "disc_rep_adult"))
sensitivity_analysis_inter_model
print(xtable(as.data.frame(sensitivity_analysis_inter_model$summary[,-11]),
     caption = "Sensitivity analysis check on robin inter-winter $M_t$ model",
     label = "tab:sensitivity_analysis_robin_age_inter", type = "latex"),
      file = "report/model_outputs/sensitivity_analysis_robin_age_inter.tex")

##############################
# INTRA-WINTER MODELLING
##############################
########## M_0 blackcap model ############
blackcap_m0_intra <- jags_model_multinomial_intra(blackcap_intra,
                                                inits = inits_m0_intra,
                                                model_type = "m0",
                                                hyper_par = list(
                                                alpha_phi = alpha_phi_juv_robin,
                                                beta_phi = beta_phi_juv_robin,
                                                alpha_p = alpha_p,
                                                beta_p = beta_p),
                                                params = c("mean_phi",
                                                           "mean_p",
                                                           "disc_obs",
                                                           "disc_rep"))
blackcap_m0_intra
print(xtable(as.data.frame(blackcap_m0_intra$summary[,-11]),
   caption = "Blackcap $M_t$ intra-winter model diagnostic and summary output",
   label = "tab:blackcap_mt_intra", type = "latex"),
      file = "report/model_outputs/blackcap_mt_intra.tex")

########## M_0 chiffchaff model ############
chiffchaff_m0_intra <- jags_model_multinomial_intra(chiffchaff_intra,
                                                inits = inits_m0_intra,
                                                model_type = "m0",
                                                hyper_par = list(
                                                alpha_phi = alpha_phi_juv_robin,
                                                beta_phi = beta_phi_juv_robin,
                                                alpha_p = alpha_p,
                                                beta_p = beta_p),
                                                params = c("mean_phi",
                                                           "mean_p",
                                                           "disc_obs",
                                                           "disc_rep"))
chiffchaff_m0_intra
print(xtable(as.data.frame(chiffchaff_m0_intra$summary[,-11]),
 caption = "Chiffchaff $M_t$ intra-winter model diagnostic and summary output",
 label = "tab:blackcap_m0_intra", type = "latex"),
      file = "report/model_outputs/chiffchaff_mt_intra.tex")

########## M_0 robin model ############
robin_m0_intra <- jags_model_multinomial_intra(robin_intra,
                                               inits = inits_m0_intra,
                                               model_type = "m0",
                                               hyper_par = list(
                                               alpha_phi = alpha_phi_juv_robin,
                                               beta_phi = beta_phi_juv_robin,
                                               alpha_p = alpha_p,
                                               beta_p = beta_p),
                                               params = c("mean_phi",
                                                          "mean_p",
                                                          "disc_obs",
                                                          "disc_rep"))
print(xtable(as.data.frame(robin_m0_intra$summary[,-11]),
       caption = "Robin $M_t$ intra-winter model diagnostic and summary output",
       label = "tab:robin_mt_intra", type = "latex"),
      file = "report/model_outputs/robin_mt_intra.tex")

########## M_t blackcap model ############
blackcap_mt_intra <- jags_model_multinomial_intra(blackcap_intra,
                                            inits = inits_mt_intra,
                                            model_type = "mt",
                                            hyper_par = list(
                                            alpha_phi = alpha_phi_juv_blackcap,
                                            beta_phi = beta_phi_juv_blackcap,
                                            alpha_p = alpha_p,
                                            beta_p = beta_p),
                                            params = c("phi",
                                                       "p",
                                                       "disc_obs",
                                                       "disc_rep"))
blackcap_mt_intra
print(xtable(as.data.frame(blackcap_mt_intra$summary[,-11]),
   caption = "Blackcap $M_t$ intra-winter model diagnostic and summary output",
   label = "tab:blackcap_m0_intra", type = "latex"),
      file = "report/model_outputs/blackcap_mt_intra.tex")

########## M_t chiffchaff model ############
chiffchaff_mt_intra <- jags_model_multinomial_intra(chiffchaff_intra,
                                          inits = inits_mt_intra,
                                          model_type = "mt",
                                          hyper_par = list(
                                          alpha_phi = alpha_phi_juv_chiffchaff,
                                          beta_phi = beta_phi_juv_chiffchaff,
                                          alpha_p = alpha_p,
                                          beta_p = beta_p),
                                          params = c("phi", "p",
                                                     "disc_obs",
                                                     "disc_rep"))
chiffchaff_mt_intra
print(xtable(as.data.frame(chiffchaff_mt_intra$summary[,-11]),
 caption = "Chiffchaff $M_t$ intra-winter model diagnostic and summary output",
 label = "tab:chiffchaff_mt_intra", type = "latex"),
      file = "report/model_outputs/chiffchaff_mt_intra.tex")

########## M_t robin model ############
robin_mt_intra <- jags_model_multinomial_intra(robin_intra,
                                               inits = inits_mt_intra,
                                               model_type = "mt",
                                               hyper_par = list(
                                               alpha_phi = alpha_phi_juv_robin,
                                               beta_phi = beta_phi_juv_robin,
                                               alpha_p = alpha_p,
                                               beta_p = beta_p),
                                               params = c("phi", "p",
                                                          "disc_obs",
                                                          "disc_rep"))
robin_mt_intra
print(xtable(as.data.frame(robin_mt_intra$summary[,-11]),
       caption = "Robin $M_t$ intra-winter model diagnostic and summary output",
       label = "tab:robin_m0_intra", type = "latex"),
      file = "report/model_outputs/robin_mt_intra.tex")

##############################
# INTRA-WINTER SENSITIVITY ANALYSYS
##############################
sensitivity_analysis_intra_model <- jags_model_multinomial_intra(robin_intra,
                                             inits = inits_mt_intra,
                                             model_type = "mt",
                                             hyper_par = list(
                                             alpha_phi = alpha_phi_juv_robin + 1,
                                             beta_phi = beta_phi_juv_robin + 1,
                                             alpha_p = alpha_p + 1,
                                             beta_p = beta_p + 1),
                                             params = c("phi", "p",
                                                        "disc_obs",
                                                        "disc_rep"))
sensitivity_analysis_intra_model
print(xtable(as.data.frame(sensitivity_analysis_intra_model$summary[,-11]),
       caption = "Sensitivity analysis check on robin intra-winter $M_t$ model",
       label = "tab:sensitivity_analysis_robin_age_intra", type = "latex"),
      file = "report/model_outputs/sensitivity_analysis_robin_intra.tex")

##############################
# INTER-WINTER MODELLING CHECKS
##############################

# Create a data frame for model comparison between three species, storing the
# DIC values for two models (M_0 and M_t)
model_comparison <- data.frame(
  Species = c("Blackcap", "Chiffchaff", "Robin"),
  M_0 = c(blackcap_m0_age_inter$DIC,
          chiffchaff_m0_age_inter$DIC,
          robin_m0_age_inter$DIC),
  M_t = c(blackcap_mt_age_inter$DIC,
          chiffchaff_mt_age_inter$DIC,
          robin_mt_age_inter$DIC)
)
model_comparison
print(xtable(model_comparison, type = "latex",
             caption = "Inter-winter model comparison",
             label = "tab:model_comparison_inter"),
      file = "report/model_outputs/model_comparison_inter.tex")


# Function to calculate Bayesian p-values
calc_pvalue <- function(model, suffix = ""){
  return(sum(model$sims.list[[paste0("disc_obs", suffix)]] >
               model$sims.list[[paste0("disc_rep", suffix)]]) /
           length(model$sims.list[[paste0("disc_rep", suffix)]]))
}

# Creating a dataframe of all p-values
p_values <- data.frame(
  species = c("Blackcap", "Chiffchaff", "Robin"),
  juveniles_m0 = c(calc_pvalue(blackcap_m0_age_inter, "_juv"),
                   calc_pvalue(chiffchaff_m0_age_inter, "_juv"),
                   calc_pvalue(robin_m0_age_inter, "_juv")),
  adults_m0 = c(calc_pvalue(blackcap_m0_age_inter, "_adult"),
                calc_pvalue(chiffchaff_m0_age_inter, "_adult"),
                calc_pvalue(robin_m0_age_inter, "_adult")),
  juveniles_mt = c(calc_pvalue(blackcap_mt_age_inter, "_juv"),
                   calc_pvalue(chiffchaff_mt_age_inter, "_juv"),
                   calc_pvalue(robin_mt_age_inter, "_juv")),
  adults_mt = c(calc_pvalue(blackcap_mt_age_inter, "_adult"),
                calc_pvalue(chiffchaff_mt_age_inter, "_adult"),
                calc_pvalue(robin_mt_age_inter, "_adult"))
)

## Add pooled p-values
p_values <- cbind(p_values[1:3],
                  overall_m0 = rowMeans(p_values[, c("juveniles_m0",
                                                     "adults_m0")]),
                  p_values[4:5],
                  overall_mt = rowMeans(p_values[, c("juveniles_mt",
                                                     "adults_mt")]))

p_values

# Creating a LaTeX table using xtable and saving it to a file
latex_pvalues <- xtable(p_values,
      caption="Inter-winter model goodness of fit results (Bayesian p-values)",
      label="tab:gof_results_inter")
print(latex_pvalues, type="latex",
      file="report/model_outputs/gof_results_age_inter.tex")

##############################
# INTRA-WINTER MODEL CHECKS
##############################

# Create a data frame for model comparison between three species, storing the
# DIC values for two models (M_0 and M_t)
model_comparison_intra <- data.frame(
  Species = c("Blackcap", "Chiffchaff", "Robin"),
  M_0 = c(blackcap_m0_intra$DIC, chiffchaff_m0_intra$DIC,
          robin_m0_intra$DIC),
  M_t = c(blackcap_mt_intra$DIC, chiffchaff_mt_intra$DIC,
          robin_mt_intra$DIC)
)
model_comparison_intra
print(xtable(model_comparison_intra, type = "latex",
             caption = "Intra-winter model comparisons",
             label = "tab:model_comparison_intra"),
      file = "report/model_outputs/model_comparison_intra.tex")

# Creating a dataframe of all Bayesian p-values
p_values_intra <- data.frame(
  species = c("Blackcap", "Chiffchaff", "Robin"),
  
  m0 = c(calc_pvalue(blackcap_m0_intra, ""),
         calc_pvalue(chiffchaff_m0_intra, ""),
         calc_pvalue(robin_m0_intra, "")),
  
  mt = c(calc_pvalue(blackcap_mt_intra, ""),
         calc_pvalue(chiffchaff_mt_intra, ""),
         calc_pvalue(robin_mt_intra, ""))
)
p_values_intra
latex_pvalues_intra <- xtable(p_values_intra,
      caption="Intra-winter model goodness of fit results (Bayesian p-values)",
      label="tab:gof_results_intra")
print(latex_pvalues_intra, type="latex",
      file="report/model_outputs/gof_results_intra.tex")

##############################
# INTER-WINTER MODEL RESULTS
##############################

species_list <- list(blackcap_mt_age_inter, chiffchaff_mt_age_inter,
                     robin_mt_age_inter)
species_names <- c("Blackcap", "Chiffchaff", "Robin")

# Path to save the plots
path_to_save <- "report/images/"
# Iterate through the first three species in the list
for (i in 1:3) {
  # Extract the result for the current species
  model_result <- species_list[[i]]
  
  # Determine the number of time points in the adult survival probability
  T <- length(model_result$mean$phi_adult)
  
  # Initialize vectors to store lower and upper quantiles for adults,
  # juveniles, and their differences
  lower_adult <- upper_adult <- lower_juv <- upper_juv <- lower_diff <- upper_diff <- numeric(T)
  
  # Calculate 2.5% and 97.5% quantiles for adult, juvenile, and difference
  # in survival probabilities
  for (t in 1:T) {
    lower_adult[t] <- quantile(model_result$sims.list$phi_adult[,t], 0.025)
    upper_adult[t] <- quantile(model_result$sims.list$phi_adult[,t], 0.975)
    lower_juv[t] <- quantile(model_result$sims.list$phi_juv[,t], 0.025)
    upper_juv[t] <- quantile(model_result$sims.list$phi_juv[,t], 0.975)
    lower_diff[t] <- quantile(model_result$sims.list$phi_diff[,t], 0.025)
    upper_diff[t] <- quantile(model_result$sims.list$phi_diff[,t], 0.975)
  }
  
  # Define the years for the x-axis
  years <- seq(2007, 2016)
  
  # Iterate through the categories (adult, juvenile, difference) to plot
  # the survival probabilities
  for (j in c("adult", "juv", "diff")) {
    # Extract y-values and corresponding lower and upper bounds
    y_values <- model_result$mean[[paste0("phi_", j)]]
    lower_values <- get(paste0("lower_", j))
    upper_values <- get(paste0("upper_", j))
    
    # Create a data frame to store the values for plotting
    df <- data.frame(Year = years, Value = y_values, Lower = lower_values,
                     Upper = upper_values)
    
    # Create the plot using ggplot2
    p <- ggplot(df, aes(x = Year, y = Value)) +
      geom_line(color = "blue") + # Line plot
      geom_point(shape = 16, color = "blue") + # Points
      geom_ribbon(aes(ymin = Lower, ymax = Upper),
                  fill = "red", alpha = 0.2) + # Confidence band
      (if (j == "diff") ylim(-1, 1) else ylim(0, 1)) + # Set y limits
      labs(y = "Probability", x = "") + # Labels
      custom_theme +
      theme(axis.text.x = element_text(angle = 45,
                                       margin = margin(t = 20, r = 0.75), 
                                       size = size * (3.5/4)),
            panel.grid.major.x = element_line(color = "grey"),
            panel.grid.minor.x = element_blank(), 
            plot.margin = unit(c(0, 0.75, 0, 0),
                               "inches")) +
      scale_x_continuous(breaks = years, labels = function(x) paste0("W", x)) +
      xlab("Winter period")
    
    
    print(p)
    ggsave(filename = paste0(path_to_save, species_names[i], "_", j,
                            "_inter.png"), plot = p, width = 12, height = 8)
  }
}


##############################
# INTRA-WINTER MODEL RESULTS
##############################

species_list_intra <- list(blackcap_mt_intra, chiffchaff_mt_intra,
                           robin_mt_intra)
# Iterate through the first three species in the list
for (i in 1:3) {
  # Extract the result for the current species
  model_result <- species_list_intra[[i]]
  
  # Determine the number of time points in the survival probability
  T <- length(model_result$mean$phi)
  
  # Initialise vectors to store lower and upper quantiles for survival
  # probabilities
  lower_phi <- upper_phi <- numeric(T)
  
  # Calculate 2.5% and 97.5% quantiles for survival probabilities
  for (t in 1:T) {
    lower_phi[t] <- quantile(model_result$sims.list$phi[,t], 0.025)
    upper_phi[t] <- quantile(model_result$sims.list$phi[,t], 0.975)
  }
  
  y_values <- model_result$mean$phi
  
  # Define the months for the x-axis
  months <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
  
  # Create a data frame to store the values for plotting
  df <- data.frame(Month = months, Value = y_values, Lower = lower_phi,
                   Upper = upper_phi)
  
  # Order the months factor to follow the given sequence
  df$Month <- factor(df$Month, levels = months)
  
  p <- ggplot(df, aes(x = Month, y = Value, group = 1)) +
    geom_line(color = "blue") +
    geom_point(shape = 16, color = "blue") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper),
                fill = "red", alpha = 0.2) + # Confidence band
    ylim(0, 1) +
    labs(y = "Probability", x = "") +
    xlab("Month") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45,
                                     margin = margin(t = 15, r = 20), 
                                     size = size),
          panel.grid.major.x = element_line(color = "grey"),
          panel.grid.minor.x = element_blank())
  
  print(p)
  
  ggsave(filename = paste0(path_to_save, species_names[i], "_", "intra.png"),
         plot = p, width = 10, height = 8)
}
