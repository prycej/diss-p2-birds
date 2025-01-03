
ni <-50000
nt <-3
nb <-2000
nc <-3

create_f1_values <- function(df) {
  f1 <- apply(X = df %>% select(-c(id,age_at_ringing,species)),
              MARGIN = 1, FUN = function(x) min(which(x == 1)))

  f1
}

create_z_init_matrix <- function(df, f1) {
  z.init <- matrix(NA, nrow = nrow(df),
                   ncol = ncol(df %>% select(-c(id,age_at_ringing,species))))

  for(i in 1:dim(z.init)[1]){
    z.init[i, f1[i]:dim(z.init)[2]] <- 1
    z.init[i,f1[i]] <- NA
  }

  z.init
}

create_survival_probability_jags_model <- function(df, inits, species, f1,
                                                   model_type, group = NULL,
                                                   params =  c("phi", "p"),
                                                   inter_intra) {
   if ((model_type == "mt" | model_type == "m0") & missing(group)) {
    jags.data<-list(y=df %>%
                      select(-c(id,age_at_ringing,species)), f1 = f1,
                    nInd=dim(df)[1],
                    n.occasions= dim(df %>%
                                       select(-c(id,age_at_ringing,species)))[2])
    } else {
      jags.data<-list(y=df %>%
                      select(-c(id,age_at_ringing,species)), f1 = f1,
                      nInd=dim(df)[1],
                      group=group,
                      n.occasions= dim(df %>%
                                         select(-c(id,age_at_ringing,species)))[2])

    }



  model <- jags.model(paste0("models/survival_",model_type,".jags"), data = jags.data,
                      inits=inits, n.chains = 3)
  saveRDS(model, paste0("jag_runs/survival/",species,"_",model_type,"_",inter_intra,"_model"))
  update(model, nb)
  samples <- coda.samples(model, variable.names = params, n.iter = ni, thin =nt)
  saveRDS(samples, paste0("jag_runs/survival/",species,"_", model_type,"_",inter_intra,"_samples"))
  return(samples)
}

survival_inits_m0 <-function(){list(z=create_z_init_matrix(blackcap_inter,
                                                           create_f1_values(blackcap_inter)),
                                    phi=runif(1,0,1),
                                    p =runif(1,0,1))}

blackcap_samples <- create_survival_probability_jags_model(blackcap_inter,
                                       'blackcap',
                                       create_f1_values(blackcap_inter),
                                       inits = survival_inits_m0(),
                                       model_type="m0", inter_intra = 'inter')

survival_inits_m0 <-function(){list(z=create_z_init_matrix(chiffchaff_inter,
                                                           create_f1_values(chiffchaff_inter)),
                                    phi=runif(1,0,1),
                                    p =runif(1,0,1))}

chiffchaff_samples <- create_survival_probability_jags_model(chiffchaff_inter,
                                                             survival_inits_m0(),
                                       'chiffchaff',
                                       create_f1_values(chiffchaff_inter),
                                       model_type="m0", inter_intra = 'inter')

survival_inits_m0 <-function(){list(z=create_z_init_matrix(robin_inter,
                                                           create_f1_values(robin_inter)),
                                    phi=runif(1,0,1),
                                    p =runif(1,0,1))}

robin_samples <- create_survival_probability_jags_model(robin_inter,
                                                        survival_inits_m0(),
                                       'robin',
                                       create_f1_values(robin_inter),
                                       model_type="m0", inter_intra = 'inter')


survival_inits_mt <-function(){list(z=create_z_init_matrix(blackcap_inter,
                                                           create_f1_values(blackcap_inter)))}

#,
#                                    alpha=runif(ncol(blackcap_inter %>%
#                                                       select(-c(id,
#                                                                 age_at_ringing,
#                                                                 species))),0,1),
#                                    beta=runif(ncol(blackcap_inter %>%
#                                                      select(-c(id,
#                                                                age_at_ringing,
#                                                                species)))),0,1
#))}
#
blackcap_samples_mt <- create_survival_probability_jags_model(blackcap_inter,
                                                           survival_inits_mt(),
                                                           'blackcap',
                                                           create_f1_values(blackcap_inter),
                                                           model_type = "mt",
                                                           inter_intra = 'inter')

plot(blackcap_samples_mt)
MCMCsummary(blackcap_samples_mt)
mcmc_trace(blackcap_samples_mt)
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


survival_inits_mt <-function(){list(z=create_z_init_matrix((chiffchaff_inter), create_f1_values(chiffchaff_inter)))}

chiffchaff_samples_mt <- create_survival_probability_jags_model(chiffchaff_inter,
                                                           survival_inits_mt(),
                                                           'blackcap',
                                                           create_f1_values(chiffchaff_inter),
                                                           model_type = "mt",
                                                           inter_intra = 'inter')

MCMCplot(blackcap_samples_mt, chiffchaff_samples_mt,
         params = c("phi"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')

survival_inits_mt <-function(){list(z=create_z_init_matrix(robin_inter, create_f1_values(robin_inter)))}

robin_samples_mt <- create_survival_probability_jags_model(robin_inter,
                                                                survival_inits_mt(),
                                                                'blackcap',
                                                                create_f1_values(robin_inter),
                                                                model_type = "mt",
                                                           inter_intra = 'inter')
MCMCplot(robin_samples_mt,
         params = c("p"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')


MCMCplot(robin_samples_mt,
         params = c("phi"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')

MCMCtrace(robin_samples_mt,
          params = c("phi[1]", "phi[2]", "phi[3]", "phi[4]", "phi[5]", "phi[6]",
                     "phi[7]", "phi[8]", "phi[9]", "phi[10]"),
          ISB = FALSE,
          exact = TRUE,
          priors = runif(15000,0,1),
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

#### model fixed group and time effects
survival_inits_mth <-function(){list(z=create_z_init_matrix(blackcap_inter, create_f1_values(blackcap_inter)))}

robin_samples_mth <- create_survival_probability_jags_model(blackcap_inter,
                                                           survival_inits_mth(),
                                                           'robin',
                                                           f1 = create_f1_values(blackcap_inter),
                                                           model_type = "mth",
                                                           group = as.factor(blackcap_inter$age_at_ringing),
                                                           params = c("phi.g1", "phi.g2", "phi.g3"),
                                                           inter_intra = 'inter')

MCMCplot(as.matrix(robin_samples_mth),
         params = c("phi.g2"),
         horiz = FALSE,
         ylab = 'PARAMETER ESTIMATE')



