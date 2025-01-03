library(gridExtra)



ni <-10000
nt <-3
nb <-2000
nc <-3



birds_mt_intra <- jags_model_multinomial(
  marr = marray(robin_intra[,5:ncol(robin_intra)]),
  inits = survival_inits_m0(),
  model_type="mt",
  inter_intra = 'inter',
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi", "p"))

MCMCplot(birds_m0_intra, horiz = FALSE, params = "phi")

blackcap_mt_intra <- jags_model_multinomial(
  marr = marray(blackcap_intra[,5:ncol(blackcap_intra)]),
  inits = survival_inits_m0(),
  model_type="mt",
  inter_intra = 'inter',
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi", "p"))

MCMCplot(blackcap_mt_intra, horiz = FALSE, params = "phi")

chiffchaff_mt_intra <- jags_model_multinomial(
  marr = marray(chiffchaff_intra[,5:ncol(chiffchaff_intra)]),
  inits = survival_inits_m0(),
  model_type="mt",
  inter_intra = 'inter',
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi", "p"))

MCMCplot(chiffchaff_mt_intra, horiz = FALSE, params = "phi")

robin_mt_intra <- jags_model_multinomial(
  marr = marray(robin_intra[,5:ncol(robin_intra)]),
  inits = survival_inits_m0(),
  model_type="mt",
  inter_intra = 'inter',
  hyper_par = list(alpha_phi = 4, beta_phi = 16, alpha_p = 4, beta_p = 16),
  params = c("phi", "p"))

MCMCplot(robin_mt_intra, horiz = FALSE, params = "phi")

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
  inter_intra = 'inter',
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
  inter_intra = 'inter',
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
  inter_intra = 'inter',
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
  model_type="mt",
  hyper_par = list(alpha_phi = 2,
                   beta_phi = 4,
                   alpha_p = 2,
                   beta_p = 4),
  inter_intra = 'inter',
  params = c("phi", "p"))

blackcap_mt_inter
plot(blackcap_mt_inter)
MCMCsummary(blackcap_mt_inter)
mcmc_trace(blackcap_mt_inter)
mcmc_areas(blackcap_mt_inter, pars = c("p[1]", "p[2]", "p[3]", "p[4]", "p[5]",
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


chiffchaff_mt_inter <- jags_model_multinomial(marr = marray(chiffchaff_inter[,4:ncol(chiffchaff_inter)]),
                                                hyper_par = list(alpha_phi = 2,
                                                                 beta_phi = 4,
                                                                 alpha_p = 2,
                                                                 beta_p = 4),
                                                #inits = survival_inits_mt,
                                                model_type="mt",
                                                inter_intra = 'inter',
                                                params = c("phi", "p"))
chiffchaff_mt_inter
plot(chiffchaff_samples_mt)

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
                                                inter_intra = 'inter',
                                                params = c("phi", "p"))

robin_mt_inter
plot(robin_mt_inter)

MCMCplot(blackcap_samples_mt, robin_samples_mt,
         params = c("phi"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')


birds_mt <- jags_model_multinomial(marr = marray(birds_inter[,4:ncol(birds_inter)]),
                                          hyper_par = list(alpha_phi = 4,
                                                           beta_phi = 16,
                                                           alpha_p = 4,
                                                           beta_p = 16),
                                          #inits = survival_inits_mt,
                                          model_type="mt",
                                          inter_intra = 'inter',
                                          params = c("phi", "p"))

MCMCtrace(birds_mt,
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
  inter_intra = 'inter',
  hyper_par = list(alpha = 4, beta= 16),
  params = c("phi.juv0", "phi.adult0", "p0"))

birds_m0_age


blackcap_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((blackcap_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(blackcap_inter)]),
  marr.adult = marray((blackcap_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(blackcap_inter)]),
  inits = NULL,
  model_type = "m0",
  inter_intra = 'inter',
  hyper_par = list(alpha = 4, beta= 16),
  params = c("phi.juv0", "phi.adult0", "p0"))

blackcap_m0_age


chiffchaff_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((chiffchaff_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(chiffchaff_inter)]),
  marr.adult = marray((chiaffchaff_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(chiffchaff_inter)]),
  inits = NULL,
  model_type = "m0",
  inter_intra = 'inter',
  hyper_par = list(alpha = 4, beta= 16),
  params = c("phi.juv0", "phi.adult0", "p0"))

chiffchaff_m0_age

robin_m0_age <- jags_model_multinomial_age(
  marr.juv = marray((robin_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(robin_inter)]),
  marr.adult = marray((robin_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(robin_inter)]),
  inits = NULL,
  model_type = "m0",
  inter_intra = 'inter',
  hyper_par = list(alpha = 4, beta= 16),
  params = c("phi.juv0", "phi.adult0", "p0"))

robin_m0_age


blackcap_mt_age <- jags_model_multinomial_age(
  marr.juv = marray((blackcap_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(blackcap_inter)]),
  marr.adult = marray((blackcap_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(blackcap_inter)]),
                                      inits = NULL,
                                      inter_intra = 'inter',
                                      hyper_par = list(alpha = 4, beta= 16),
                                      params = c("phi.juv", "phi.adult", "p"))

blackcap_mt_age

chiffchaff_mt_age <- jags_model_multinomial_age(
  marr.juv = marray((chiffchaff_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(chiffchaff_inter)]),
  marr.adult = marray((chiffchaff_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(chiffchaff_inter)]),
                                              inits = NULL,
                                              inter_intra = 'inter',
                                              hyper_par = list(alpha = 4, beta= 16),
                                              params = c("phi.juv", "phi.adult", "p"))

chiffchaff_mt_age

robin_mt_age <- jags_model_multinomial_age(
  marr.juv = marray((robin_inter %>%
                       filter(age_at_ringing == "juvenile"))[,4:ncol(robin_inter)]),
  marr.adult = marray((robin_inter %>%
                      filter(age_at_ringing == "adult"))[,4:ncol(robin_inter)]),
                                                inits = NULL,
                                                inter_intra = 'inter',
                                                hyper_par = list(alpha = 4, beta= 16),
                                                params = c("phi.juv", "phi.adult", "p"))

robin_mt_age

# Convert mcmc output to data frame
data_df <- data.frame(as.matrix(birds_mt_age$samples))

# Reshape data to long format
data_long <- data_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(
    parameter_type = str_extract(parameter, "^[a-z]+"),
    age_type = ifelse(str_detect(parameter, "\\.ad\\."), "ad", ifelse(str_detect(parameter, "\\.juv\\."), "juv", NA)),
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
