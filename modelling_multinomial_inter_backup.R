library(gridExtra)



ni <-10000
nt <-3
nb <-2000
nc <-3

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
marr <- marray(blackcap_inter[,4:ncol(blackcap_inter)])

( releases <- rowSums(marr) )

jdata <- list(marr = marr, R = releases, n.occasions=ncol(marr))
str(jdata)

wanted <- c("phi", "p")

( out_phi._p. <- jags(jdata, NULL, c("phi0", "p0"),
                      "models/multinomial_m0.jags", DIC=FALSE,
                      n.chains=nc, n.adapt=1000, n.iter=ni, parallel = TRUE) )

( out_phi._p. <- jags(jdata, NULL, wanted,
                      "models/multinomial_mt.jags", DIC=FALSE,
                      n.chains=nc, n.adapt=1000, n.iter=ni, parallel = TRUE) )


jags_model_multinomial <- function(marr, hyper_par, inits, model_type,
                                   group = NULL, params =  c("phi", "p")) {
   releases <- rowSums(marr)
   if ((model_type == "mt" | model_type == "m0") & missing(group)) {
     jdata <- list(marr = marr, R = releases, n.occasions=ncol(marr),
                   alpha_phi = hyper_par$alpha_phi,
                   beta_phi = hyper_par$beta_phi,
                   alpha_p = hyper_par$alpha_p,
                   beta_p = hyper_par$beta_p)
    } else {
     jdata <- list(marr = marr, R = releases, n.occasions=ncol(marr),
                   alpha_phi = hyper_par$alpha_phi,
                   beta_phi = hyper_par$beta_phi,
                   alpha_p = hyper_par$alpha_p,
                   beta_p = hyper_par$beta_p)
    }

  out <- jags(jdata, NULL, params,
                paste0("models/multinomial_",model_type,".jags"), DIC=TRUE,
                n.chains=3, n.adapt=1000, n.iter=10000, parallel = TRUE,
              seed = 42)
  out
}

jags_model_multinomial_age <- function(marr.adult, marr.juv, hyper_par, inits,
                                       model_type = "mt",
                                       params =  c("phi.juv", "phi.ad", "p")) {
  releases.juv <- rowSums(marr.juv)
  releases.adult <- rowSums(marr.adult)
    jdata <- list(marr.juv = marr.juv, marr.adult = marr.adult,
                  releases.juv = releases.juv, releases.adult = releases.adult,
                  n.occasions=ncol(marr.juv),
                  alpha_phi_j = hyper_par$alpha_phi_j,
                  beta_phi_j = hyper_par$beta_phi_j,
                  alpha_phi_ad = hyper_par$alpha_phi_j,
                  beta_phi_ad = hyper_par$beta_phi_ad,
                  alpha_p = hyper_par$alpha_p,
                  beta_p = hyper_par$beta_p)


  out <- jags(jdata, NULL, params,
              paste0("models/multinomial_",model_type,"_age.jags"), DIC=TRUE,
              n.chains=3, n.adapt=1000, n.iter=10000, parallel = TRUE,  seed = 42)
  out
}

survival_inits_m0 <-function(){list(phi=runif(3,0,1), p =runif(3,0,1))}

birds_m0 <- jags_model_multinomial(
  marr = marray(birds_inter[,4:ncol(birds_inter)]),
  inits = survival_inits_m0(),
  model_type="m0",
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi0", "p0"))

#MCMCtrace(birds_samples_m0, type = 'density', pdf = TRUE, filename = "report/images/birds_m0.pdf")
png("report/images/birds_m0.png")
par(mfrow=c(1,2))
plot(density(as.matrix(birds_m0$samples)[,"phi0"]))
plot(density(as.matrix(birds_m0$samples)[,"p0"]))
dev.off()

survival_inits_m0 <-function(){list(phi=runif(3,0,1), p =runif(3,0,1))}

blackcap_m0 <- jags_model_multinomial(
  marr = marray(blackcap_inter[,4:ncol(blackcap_inter)]),
  inits = survival_inits_m0(),
  model_type="m0",
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi0", "p0"))

blackcap_m0
png("report/images/blackcap_m0.png")
par(mfrow=c(1,2))
plot(density(as.matrix(blackcap_m0$samples)[,"phi0"]))
plot(density(as.matrix(blackcap_m0$samples)[,"p0"]))
dev.off()


survival_inits_m0 <-function(){phi=runif(1,0,1), p =runif(1,0,1))}

chiffchaff_m0 <- jags_model_multinomial(
  marr = marray(chiffchaff_inter[,4:ncol(chiffchaff_inter)]),
  inits = survival_inits_m0(),
  model_type="m0",
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi0", "p0"))

chiffchaff_m0
png("report/images/chiffchaff_m0.png")
par(mfrow=c(1,2))
plot(density(as.matrix(chiffchaff_m0$samples)[,"phi0"]))
plot(density(as.matrix(chiffchaff_m0$samples)[,"p0"]))
dev.off()


df <- data.frame(value = as.matrix(birds_samples_m0$samples))
ggplot(df, aes(x = phi0)) +
  geom_density() +
  labs(x = "Value", y = "Density", title = "Density plot") +
  theme_minimal()
ggsave("report/images/birds_m0.png")


survival_inits_m0 <-function(){list(z=create_z_init_matrix(robin_inter,
                                                           create_f1_values(robin_inter)),
                                    phi=runif(1,0,1),
                                    p =runif(1,0,1))}

robin_m0 <- jags_model_multinomial(
  marr = marray(robin_inter[,4:ncol(robin_inter)]),
  inits = survival_inits_m0(),
  model_type="m0",
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi0", "p0"))

robin_m0
png("report/images/robin_m0.png")
par(mfrow=c(1,2))
plot(density(as.matrix(robin_m0$samples)[,"phi0"]))
plot(density(as.matrix(robin_m0$samples)[,"p0"]))
dev.off()


# Define the color for each species using hex codes
color_map <- c("Robin" = "#619CFF", "Blackcap" = "#F8766D",
               "Chiffchaff" = "#00BA38", "Birds" = "black")

# Convert matrices to data frames and add species column
df_phi <- bind_rows(
  data.frame(species = "Robin", value = as.matrix(robin_m0$samples)[,"phi0"]),
  data.frame(species = "Blackcap", value = as.matrix(blackcap_m0$samples)[,"phi0"]),
  data.frame(species = "Chiffchaff", value = as.matrix(chiffchaff_m0$samples)[,"phi0"]),
  data.frame(species = "Birds", value = as.matrix(birds_m0$samples)[,"phi0"]) # new data set
)

df_p <- bind_rows(
  data.frame(species = "Robin", value = as.matrix(robin_m0$samples)[,"p0"]),
  data.frame(species = "Blackcap", value = as.matrix(blackcap_m0$samples)[,"p0"]),
  data.frame(species = "Chiffchaff", value = as.matrix(chiffchaff_m0$samples)[,"p0"]),
  data.frame(species = "Birds", value = as.matrix(birds_m0$samples)[,"p0"]) # new data set
)

# Create the density plots
plot_phi <- ggplot(df_phi, aes(x = value, color = species)) +
  geom_density() +
  scale_color_manual(values = color_map) +
  xlim(0, 1) +
  labs(x = "Value", y = "Density", title = "Density Plot for phi0") +
  theme_minimal()

plot_p <- ggplot(df_p, aes(x = value, color = species)) +
  geom_density() +
  scale_color_manual(values = color_map) +
  xlim(0, 1) +
  labs(x = "Value", y = "Density", title = "Density Plot for p0") +
  theme_minimal()

# Arrange the plots side by side
combined_plot <- grid.arrange(plot_phi, plot_p, ncol = 2)
ggsave("report/images/m0_models.png", combined_plot)



blackcap_mt_inter <- jags_model_multinomial(
  marr = marray(blackcap_inter[,4:ncol(blackcap_inter)]),
  inits = survival_inits_mt(),
  model_type="mt",
  params = c("phi", "p"))

plot(blackcap_mt_inter)
MCMCsummary(blackcap_mt_inter)
mcmc_trace(blackcap_mt_inter)
mcmc_areas(blackcap_samples_mt, pars = c("p[1]", "p[2]", "p[3]", "p[4]", "p[5]",
                                         "p[6]", "p[7]", "p[8]", "p[9]", "p[10]"))

MCMCplot(blackcap_samples_mt,
         params = c("p"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE',
         HPD = TRUE)

MCMCplot(blackcap_samples_mt,
         params = c("phi"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE',
         )


survival_inits_mt < list()

chiffchaff_mt_inter <- jags_model_multinomial(marr = marray(chiffchaff_inter[,4:ncol(chiffchaff_inter)]),
                                                hyper_par = list(alpha_phi = 2,
                                                                 beta_phi = 4,
                                                                 alpha_p = 2,
                                                                 beta_p = 4),
                                                #inits = survival_inits_mt,
                                                model_type="mt",
                                                params = c("phi", "p"))
plot(chiffchaff_mt_inter)

MCMCplot(blackcap_samples_mt, chiffchaff_samples_mt,
         params = c("phi"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')

survival_inits_mt <-function(){list(z=create_z_init_matrix(robin_inter,
                                                           create_f1_values(robin_inter)))}

robin_mt_inter <- jags_model_multinomial(marr = marray(robin_inter[,4:ncol(robin_inter)]),
                                                hyper_par = list(alpha_phi = 2,
                                                                 beta_phi = 4,
                                                                 alpha_p = 2,
                                                                 beta_p = 4),
                                                #inits = survival_inits_mt,
                                                model_type="mt",
                                                params = c("phi", "p"))
plot(robin_mt_inter)

MCMCplot(blackcap_mt_inter, robin_mt_inter,
         params = c("phi"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')


birds_mt_inter <- jags_model_multinomial(marr = marray(birds_inter[,4:ncol(birds_inter)]),
                                          hyper_par = list(alpha_phi = 4,
                                                           beta_phi = 16,
                                                           alpha_p = 4,
                                                           beta_p = 16),
                                          #inits = survival_inits_mt,
                                          model_type="mt",
                                          params = c("phi", "p"))

MCMCtrace(birds_mt_inter,
          params = c("phi[1]", "phi[2]", "phi[3]", "phi[4]", "phi[5]", "phi[6]",
                     "phi[7]", "phi[8]", "phi[9]", "phi[10]"),
          ISB = FALSE,
          exact = TRUE,
          priors = rbeta(15000,4,16),
          pdf = FALSE,
          Rhat = TRUE)

#### model fixed group and time effects
survival_inits_mth <-function(){list(z=create_z_init_matrix(blackcap_inter,
                                                            create_f1_values(blackcap_inter)))}


# Convert mcmc output to data frames
blackcap_df <- data.frame(as.matrix(blackcap_samples_mt$samples))
chiffchaff_df <- data.frame(as.matrix(chiffchaff_samples_mt$samples))
robin_df <- data.frame(as.matrix(robin_samples_mt$samples))
birds_df <- data.frame(as.matrix(birds_samples_mt$samples)) # new data set

# Reshape data to long format
blackcap_long <- blackcap_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(species = "Blackcap",
         parameter_type = str_extract(parameter, "^[a-z]+"),
         parameter_number = as.numeric(str_extract(parameter, "(?<=\\.)\\d+(?=\\.)")))

chiffchaff_long <- chiffchaff_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(species = "Chiffchaff",
         parameter_type = str_extract(parameter, "^[a-z]+"),
         parameter_number = as.numeric(str_extract(parameter, "(?<=\\.)\\d+(?=\\.)")))

robin_long <- robin_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(species = "Robin",
         parameter_type = str_extract(parameter, "^[a-z]+"),
         parameter_number = as.numeric(str_extract(parameter, "(?<=\\.)\\d+(?=\\.)")))

birds_long <- birds_df %>% # new data set
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(species = "Birds",
         parameter_type = str_extract(parameter, "^[a-z]+"),
         parameter_number = as.numeric(str_extract(parameter, "(?<=\\.)\\d+(?=\\.)")))

# Combine all data into a single data frame
all_species <- bind_rows(blackcap_long, chiffchaff_long, robin_long, birds_long) # new data set added

# Arrange by parameter type and number
all_species <- all_species %>%
  arrange(parameter_type, parameter_number) # order by descending parameter_type to get phi first

parameters <- unique(all_species$parameter[order(all_species$parameter_type,
                                                 all_species$parameter_number)])

# Split the list into "phi" and "p" parameters
phi_parameters <- parameters[str_detect(parameters, "^phi")]
p_parameters <- parameters[str_detect(parameters, "^p\\.")]

# Concatenate the lists back together with "phi" first
sorted_parameters <- c(phi_parameters, p_parameters)

# convert parameter to a factor and specify the order of the levels
all_species$parameter <- factor(all_species$parameter,
                                levels = sorted_parameters)

# Define the color for each species
color_map <- c("Blackcap" = "#F8766D", "Chiffchaff" = "#00BA38",
               "Robin" = "#619CFF", "Birds" = "black")

# Plot densities
ggplot(all_species, aes(x = value, color = species)) +
  geom_density() +
  facet_wrap(~parameter, scales = "free", dir = "h") +
  scale_color_manual(values = color_map) +
  theme_minimal() +
  labs(x = "Parameter Value", y = "Density", title = "Comparative Density Plots of MCMC Output for Different Bird Species") +
  theme(legend.position = "bottom")
ggsave("report/images/mt_models.png")

####### adding in age

birds_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((birds_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(birds_inter)]),
  marr.adult = marray((birds_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(birds_inter)]),
  inits = NULL,
  model_type = "m0",
  hyper_par = list(alpha_phi_j = 4, beta_phi_j= 16, alpha_phi_ad= 4, beta_phi_ad= 16, alpha_p = 4, beta_p= 16),
  params = c("phi.juv0", "phi.ad0", "p0"))

birds_m0_age


blackcap_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((blackcap_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(blackcap_inter)]),
  marr.adult = marray((blackcap_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(blackcap_inter)]),
  inits = NULL,
  model_type = "m0",
  hyper_par = list(alpha = 4, beta= 16),
  params = c("phi.juv0", "phi.ad0", "p0"))

blackcap_m0_age


chiffchaff_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((chiffchaff_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(chiffchaff_inter)]),
  marr.adult = marray((chiaffchaff_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(chiffchaff_inter)]),
  inits = NULL,
  model_type = "m0",
  hyper_par = list(alpha = 4, beta= 16),
  params = c("phi.juv0", "phi.ad0", "p0"))

chiffchaff_m0_age

robin_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((robin_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(robin_inter)]),
  marr.adult = marray((robin_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(robin_inter)]),
  inits = NULL,
  model_type = "m0",
  hyper_par = list(alpha_phi_j = 4, beta_phi_j= 16, alpha_phi_ad= 4, beta_phi_ad= 16, alpha_p = 4, beta_p= 16),
  params = c("phi.juv0", "phi.ad0", "p0"))

robin_m0_age


blackcap_mt_age <- jags_model_multinomial_age(
  marr.juv = marray((blackcap_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(blackcap_inter)]),
  marr.adult = marray((blackcap_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(blackcap_inter)]),
                                      inits = NULL,
  model_type ="mt",
  hyper_par = list(alpha_phi_j = 4, beta_phi_j= 16, alpha_phi_ad= 4, beta_phi_ad= 16, alpha_p = 4, beta_p= 16),
                                      params = c("phi.juv", "phi.ad", "p"))

blackcap_mt_age

chiffchaff_mt_age <- jags_model_multinomial_age(
  marr.juv = marray((chiffchaff_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(chiffchaff_inter)]),
  marr.adult = marray((chiffchaff_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(chiffchaff_inter)]),
                                              inits = NULL,
                                              hyper_par = list(alpha = 4, beta= 16),
                                              params = c("phi.juv", "phi.ad", "p"))

chiffchaff_mt_age

robin_mt_age <- jags_model_multinomial_age(
  marr.juv = marray((robin_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(robin_inter)]),
  marr.adult = marray((robin_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(robin_inter)]),
                                                inits = NULL,
                                                hyper_par = list(alpha = 4, beta= 16),
                                                params = c("phi.juv", "phi.ad", "p"))

robin_mt_age

# Convert mcmc output to data frame
data_df <- data.frame(as.matrix(birds_mt_age$samples))

# Reshape data to long format
data_long <- data_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(
    parameter_type = str_extract(parameter, "^[a-z]+"),
    age_type = ifelse(str_detect(parameter, "\\.ad\\."), "ad",
                      ifelse(str_detect(parameter, "\\.juv\\."), "juv", NA)),
    index = as.numeric(str_extract(parameter, "(?<=\\.)\\d+(?=\\.)")),
    plot_order = ifelse(parameter_type == "phi", index, index + 10)
  )

# Define custom labeller function
custom_labeller <- function(variable,value){
  return(paste(unique(data_long$parameter_type[data_long$plot_order == value]),
               unique(data_long$index[data_long$plot_order == value]),
               sep="["))
}

# Plot densities
ggplot(data_long, aes(x = value, color = age_type)) +
  geom_density() +
  facet_wrap(~ plot_order, labeller = custom_labeller, scales = "free_y", dir = "h", ncol = 5) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  labs(x = "Parameter Value", y = "Density", title = "Density Plots of MCMC Output for  MT Age") +
  theme(legend.position = "bottom")
ggsave("report/images/blackcap_mt_age.png")




CH.J <- chiffchaff_inter %>% filter(age_at_ringing == "juvenile") %>% select(c(-id, -age_at_ringing, -species))
CH.A <- chiffchaff_inter %>% filter(age_at_ringing == "adult") %>% select(c(-id, -age_at_ringing, -species))

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
CH.A.marray <- marray(CH.A.m)
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
CH.J.marray <- CH.J.R.marray + CH.J.N.marray


params <- c("phi.juv", "phi.ad", "p")
releases.juv <- rowSums(CH.J.marray)
releases.adult <- rowSums(CH.A.marray)
jdata <- list(marr.j = CH.J.marray, marr.a = CH.A.marray, releases.juv = releases.juv, releases.adult = releases.adult,
              n.occasions = dim(CH.J.marray)[2], hyper_par = list(alpha_phi_j = 4, beta_phi_j= 16, alpha_phi_ad= 4, beta_phi_ad= 16, alpha_p = 4, beta_p= 16))


samples <- jags(jdata, NULL, c("mean.phijuv", "mean.phiad", "mean.p", "phi.diff"),
     paste0("models/multinomial_m0_age_2.jags"), DIC=TRUE,
     n.chains=3, n.adapt=1000, n.iter=10000, parallel = TRUE,
     seed = 42)
MCMCplot(samples, params = c("phi.diff"))

samples <- jags(jdata, NULL, c("phi.ad", "phi.juv", "mean.p", "phi.diff"),
     paste0("models/multinomial_mt_age_2.jags"), DIC=TRUE,
     n.chains=3, n.adapt=1000, n.iter=10000, parallel = TRUE,
     seed = 42)

MCMCplot(samples, params = c("phi.diff"))


