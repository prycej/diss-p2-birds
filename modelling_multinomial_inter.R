##############################
# INTER-WINTER MODELLING
##############################
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
    phi.adult0 = runif(1),
    p0 = runif(1))
}

inits_mt_age_inter <- function(chain){
  list(
    phi.juv = runif(10, 0, 1),
    phi.adult = runif(10, 0, 1),
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
                                       params =  c("phi.juv", "phi.adult", "p")) {

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
  releases.adult <- rowSums(marr.a)
  jdata <- list(marr.j = marr.j, marr.a = marr.a,
                  releases.juv = releases.juv, releases.adult = releases.adult,
                  n.occasions=ncol(marr.j),
                  alpha_phi_j = hyper_par$alpha_phi_j,
                  beta_phi_j = hyper_par$beta_phi_j,
                  alpha_phi_ad = hyper_par$alpha_phi_j,
                  beta_phi_ad = hyper_par$beta_phi_ad,
                  alpha_p = hyper_par$alpha_p,
                  beta_p = hyper_par$beta_p)


  out <- jags(jdata, inits, params,
              paste0("models/multinomial_",model_type,"_age_2.jags"), DIC=TRUE,
              n.chains=nc, n.adapt=na, n.iter=ni, parallel = TRUE,  seed = 42)
  out
}

jags_model_multinomial_intra <- function(df, hyper_par, inits,
                                             model_type = "mt",
                                             params =  c("phi.juv", "phi.adult",
                                                         "p")) {

  marr <-  marray(df %>% select(c(-id, -age_at_ringing, -species, -winter)))
  releases <- rowSums(marr)
  jdata <- list(marr=marr,
                R = releases,
                n.occasions=ncol(marr),
                alpha_phi= hyper_par$alpha_phi,
                beta_phi = hyper_par$beta_phi,
                alpha_p = hyper_par$alpha_p,
                beta_p = hyper_par$beta_p)


  out <- jags(jdata, inits, params,
              paste0("models/multinomial_",model_type,"_intra.jags"),
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
                                                         "phi.adult0",
                                                         "p0", "mean.phi.diff",
                                                       "fit.juv", "fit.new.juv",
                                                       "fit.ad", "fit.new.ad"))

blackcap_m0_age_inter
MCMCplot(blackcap_m0_age_inter, params = c("phi.juv0", "phi.adult0", "p0",
                                           "mean.phi.diff"))
print(xtable(as.data.frame(blackcap_m0_age_inter$summary[,-11]),
    caption = "Blackcap $M_0$ inter-winter model diagnostic and summary output",
    label = "tab:blackcap_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/blackcap_m0_inter.tex")

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
                                              params = c("phi.juv0", "phi.adult0",
                                                         "p0", "mean.phi.diff",
                                                           "fit.juv",
                                                           "fit.new.juv",
                                                           "fit.ad",
                                                           "fit.new.ad"))

chiffchaff_m0_age_inter
MCMCplot(chiffchaff_m0_age_inter, params = c("phi.juv0", "phi.adult0", "p0",
                                       "mean.phi.diff"))
print(xtable(as.data.frame(chiffchaff_m0_age_inter$summary[,-11]),
  caption = "Chiffchaff $M_0$ inter-winter model diagnostic and summary output",
             label = "tab:chiffchaff_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/chiffchaff_m0_inter.tex")
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
                                                      "phi.adult0",
                                                      "p0", "mean.phi.diff",
                                                      "fit.juv", "fit.new.juv",
                                                      "fit.ad", "fit.new.ad"))

robin_m0_age_inter
print(xtable(as.data.frame(robin_m0_age_inter$summary[,-11]),
       caption = "Robin $M_0$ inter-winter model diagnostic and summary output",
        label = "tab:robin_m0_inter_summary_output", type = "latex"),
      file = "report/model_outputs/robin_m0_inter.tex")
MCMCplot(robin_m0_age_inter, params = c("phi.juv0", "phi.adult0", "p0",
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
                                              params = c("phi.juv", "phi.adult",
                                                       "p", "phi.diff",
                                                       "fit.juv", "fit.new.juv",
                                                       "fit.ad", "fit.new.ad"))
(model <- blackcap_mt_age_inter)
MCMCplot(blackcap_mt_age_inter, params = c("phi.diff"))
print(xtable(as.data.frame(blackcap_mt_age$summary[,-11]),
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
                                          params = c("phi.juv", "phi.adult", "p",
                                                     "phi.diff", "fit.juv",
                                                     "fit.new.juv", "fit.ad",
                                                     "fit.new.ad"))

(model <- chiffchaff_mt_age_inter)
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
                                           params = c("phi.juv", "phi.adult",
                                                      "p", "phi.diff",
                                                      "fit.juv", "fit.new.juv",
                                                      "fit.ad", "fit.new.ad"))

(model <- robin_mt_age_inter)
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
                                 params = c("phi",
                                            "p",
                                            "fit", "fit.new"))
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
                                 params = c("phi",
                                            "p",
                                            "fit", "fit.new"))
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
                         filter(age_at_ringing == "juvenile"))[, 4:11]),
            rep(1, nrow((birds_intra %>%
                           filter(age_at_ringing == "juvenile"))[, 4:11])))

overall_CJS(as.matrix(blackcap_intra[, 4:11]),
            rep(1, nrow(blackcap_intra[, 4:11])))
overall_CJS(as.matrix(chiffchaff_intra[, 4:11]),
            rep(1, nrow(chiffchaff_intra[, 4:11])))
overall_CJS(as.matrix(robin_intra[, 4:11]),
            rep(1, nrow(robin_intra[, 4:11])))



model_comparison <- data.frame(
  Species = c("Blackcap", "Chiffchaff", "Robin"),
  M_0 = c(blackcap_m0_age$DIC, chiffchaff_m0_age$DIC, robin_m0_age$DIC),
  M_t = c(blackcap_mt_age$DIC, chiffchaff_mt_age$DIC, robin_mt_age$DIC)
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
  juveniles_m0 = c(calc_pvalue(blackcap_m0_age, ".juv"),
                   calc_pvalue(chiffchaff_m0_age, ".juv"),
                   calc_pvalue(robin_m0_age, ".juv")),
  adults_m0 = c(calc_pvalue(blackcap_m0_age, ".ad"),
                calc_pvalue(chiffchaff_m0_age, ".ad"),
                calc_pvalue(robin_m0_age, ".ad")),
  juveniles_mt = c(calc_pvalue(blackcap_mt_age, ".juv"),
                   calc_pvalue(chiffchaff_mt_age, ".juv"),
                   calc_pvalue(robin_mt_age, ".juv")),
  adults_mt = c(calc_pvalue(blackcap_mt_age, ".ad"),
                calc_pvalue(chiffchaff_mt_age, ".ad"),
                calc_pvalue(robin_mt_age, ".ad"))
)


# Creating a LaTeX table using xtable and saving it to a file
latex_pvalues <- xtable(p_values,
      caption="Inter-winter model goodness of fit results (bayesian p-values)",
      label="tab:gof_results")
print(latex_pvalues, type="latex", file="report/model_outputs/gof_results.tex")


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
