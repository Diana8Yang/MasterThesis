## Initialise

rm(list=ls()) # Remove every previous functions and variables
gc() # Free up some space
{
  setwd("~/Desktop/MThesis")

  install.packages("rjags")
  # install.packages("runjags")
  install.packages("MASS")
  install.packages("fitdistrplus")
  library("coda")
  library("rjags")
  # library("runjags")
  library("MASS")
  library("fitdistrplus")
  library("mc2d")
}

##   data and functions
{
  source("DataCollectionChicken.R")
  # source("ExtraData.R")
}

##################################################################
############# INDIVIDUAL MODELS WITH ORIGINAL PRIORS #############
##################################################################
{
##### MIKKELA'S PREVALENCE AND CONCENTRATION MODEL #########
if(FALSE) 
  {
  ## Runs the model
  {
    Mikkela_jags <- jags.model(file = "MikkelaModelC.R",
                               data = Mikkela_all_data_c,
                               inits = initdata_c,
                               n.chains = 4,
                               n.adapt = 1000) # for adapting or tuning.
    
    update(Mikkela_jags, 10^4) # burn-in
    
    Mikkela_parametertosave <- c("pb", "pw", "pf", "pba","pfa", "pvarb", "pvarw",
                                 "vartot",  'taue', 'logitpb', 'e', 'mu_zero',
                                 'tauw','sdb')
    
    Mikkela_coda <- coda.samples(Mikkela_jags,
                                Mikkela_parametertosave,
                                10^4)
    
    Mikkela_mat <- as.matrix(Mikkela_coda)
    
    post_Mikkela_orig_prior <- list(
    logitpb_1 = Mikkela_mat[,grep("logitpb",colnames(Mikkela_mat))][, 1], # logitpb is the first month
    taue = Mikkela_mat[,grep("taue",colnames(Mikkela_mat))],
    e_1 = Mikkela_mat[,grep("\\be\\b",colnames(Mikkela_mat))][, 1], # e_1 is the first month
    pw = Mikkela_mat[,grep("pw",colnames(Mikkela_mat))],
    tauw = Mikkela_mat[,grep("tauw",colnames(Mikkela_mat))],
    mu_zero = Mikkela_mat[,grep("mu_zero",colnames(Mikkela_mat))],
    sdb = Mikkela_mat[,grep("sdb",colnames(Mikkela_mat))]
    )
  }
  
  ## Fit distributions to marginal posterior to get informed priors for Mikkela's model
  {
    # pw: Prior (not informed) for pw is U(0, 1) (Mikkela et al.)
    {
      iprior.pw <- fitdist(post_Mikkela_orig_prior$pw, "beta") #, list(shape1=0.5, shape2=0.5))
      iprior.pw$convergence # 0, i.e. successful convergence
    }
    
    # tau_e: Prior (not informed) for tau_e is gamma(0.001, 0.001)
    {
      iprior.taue <- fitdist(post_Mikkela_orig_prior$taue, "lnorm" , probs=c(0.05, 0.95), method="qme")
      iprior.taue$convergence # 0
    }
    
    # logit_pb: Prior (not informed) for logit_pb is normal(0, 0.001)
    {
      iprior.logitpb_1 <- fitdist(post_Mikkela_orig_prior$logitpb_1, "norm", method = "mle")
      iprior.logitpb_1$convergence # 0
    }
    
    # e_1: Prior (not informed) for e_1 is N(0, 0.001)
    {
      iprior.e_1 <- fitdist(post_Mikkela_orig_prior$e_1, "norm")
      iprior.e_1$convergence # 0
    }
    
    # mu_zero: Prior (not informed) for mu_zero is N(0, 0.001)
    {
      iprior.mu_zero <- fitdist(post_Mikkela_orig_prior$mu_zero, "norm")
      iprior.mu_zero$convergence # 0
    }
    
    # tau_w: Prior (not informed) for tau_w is Gamma(0.001, 0.001)
    {
      iprior.tauw <- fitdist(post_Mikkela_orig_prior$tauw, "gamma", probs=c(0.05, 0.95),
                             start=list(shape=8, scale=1), method = "qme")
      iprior.tauw$convergence # 0
    }
    
    # sdb_t: Prior (not informed) for sdb is U(0, 100)
      # use a gamma
      {
        iprior.sdb <- fitdist(post_Mikkela_orig_prior$sdb, "gamma")
        iprior.sdb$convergence # 0
      }
    {
      save_prior.pw <- rbind(c(iprior.pw$estimate[1], iprior.pw$estimate[2]),
                             c(iprior.pw$estimate[1], iprior.pw$estimate[2]))
      save_prior.taue <- rbind(c(iprior.taue$estimate[1], iprior.taue$estimate[2]),
                               c(iprior.taue$estimate[1], iprior.taue$estimate[2]))
      save_prior.logitpb_1 <- rbind(c(iprior.logitpb_1$estimate[1], iprior.logitpb_1$estimate[2]),
                                  c(iprior.logitpb_1$estimate[1], iprior.logitpb_1$estimate[2]))
      save_prior.e_1 <- rbind(c(iprior.e_1$estimate[1], iprior.e_1$estimate[2]),
                            c(iprior.e_1$estimate[1], iprior.e_1$estimate[2]))
      save_prior.mu_zero <- rbind(c(iprior.mu_zero$estimate[1], iprior.mu_zero$estimate[2]),
                                  c(iprior.mu_zero$estimate[1], iprior.mu_zero$estimate[2]))
      save_prior.tauw <- rbind(c(iprior.tauw$estimate[1], iprior.tauw$estimate[2]),
                               c(iprior.tauw$estimate[1], iprior.tauw$estimate[2]))
      save_prior.sdb <- rbind(c(iprior.sdb$estimate[1], iprior.sdb$estimate[2]),
                              c(iprior.sdb$estimate[1], iprior.sdb$estimate[2]))
      
      informed_prior_Mikkela <- list("iprior.logitpb_1" = save_prior.logitpb_1,
                                     "iprior.taue" = save_prior.taue,
                                     "iprior.e_1" = save_prior.e_1,
                                     "iprior.pw" = save_prior.pw,
                                     "iprior.tauw" = save_prior.tauw,
                                     "iprior.mu_zero" = save_prior.mu_zero,
                                     "iprior.sdb" = save_prior.sdb)
    }
  }
    # Save the lists to a RData file
    save(Mikkela_coda, file="Mikkela_coda.RData")
    save(informed_prior_Mikkela, file="Mikkela_iprior.RData")
    save(post_Mikkela_orig_prior, file="Mikkela_posterior_orig_prior.RData")
}
  load("Mikkela_coda.RData")
  load("Mikkela_iprior.RData") #informed_prior_Mikkela
  load("Mikkela_posterior_orig_prior.RData") #post_Mikkela_orig_prior

##### TEUNIS DOSE-RESPONSE MODEL ###########################
if(FALSE)
  {
## Runs the model 
{
    Teunis_jags <- jags.model(file = "TeunisModel_z.R",
                       data = drdata,
                       n.chains = 4,
                       n.adapt = 1000)
    
    update(Teunis_jags, 10^4) # burn-in
    
    parametertosaveDR <- c('mu.w', 'prior.z', 'z1', 'z2',
                           'logconc', 'ltp0',
                           'prior.conc', 'prior.p0',
                           'a', 'b', 'eta', 'r', # 'w1', 'w2',
                           'prill','prinf')
                          
    Teunis_coda <- coda.samples(Teunis_jags,
                           parametertosaveDR,
                           10^4)
    
    sammatDR <- as.matrix(Teunis_coda)
    
    post_Teunis_orig_prior <- list(mu.w = sammatDR[,grep('mu.w',colnames(sammatDR))],
              prior.z = sammatDR[,grep('prior.z', colnames(sammatDR))],
              z1 = sammatDR[,grep('z1',colnames(sammatDR))],
              z2 = sammatDR[,grep('z2',colnames(sammatDR))],
              logconc = sammatDR[,grep('logconc',colnames(sammatDR))],
              ltp0 = sammatDR[,grep('ltp0',colnames(sammatDR))],
              prior.conc = sammatDR[,grep('prior.conc',colnames(sammatDR))],
              prior.p0 = sammatDR[,grep('prior.p0',colnames(sammatDR))],
              # w1 = sammatDR[,grep('w1',colnames(sammatDR))], # 12, too high
              # w2 = sammatDR[,grep('w2',colnames(sammatDR))],
              a = sammatDR[,grep('\\ba\\b',colnames(sammatDR))],
              b = sammatDR[,grep('\\bb\\b',colnames(sammatDR))],
              eta = sammatDR[,grep('eta',colnames(sammatDR))], # e-16, too low.
              r = sammatDR[,grep('\\br\\b',colnames(sammatDR))], # 0.9, too high))
              prill = sammatDR[,grep('prill',colnames(sammatDR))],
              prinf = sammatDR[,grep('prinf',colnames(sammatDR))]) 
 
}
  
## Fit distributions to marginal posteriors to get informed priors
{
    # mu.w ~ normal
    {
      iprior.mu.w11 <- fitdist(post_Teunis_orig_prior$mu.w[,1], "norm")  
      # plot(iprior.mu.w11)
      
      iprior.mu.w12 <- fitdist(post_Teunis_orig_prior$mu.w[,4], "norm")  
      # plot(iprior.mu.w12)
      
      iprior.mu.w21 <- fitdist(post_Teunis_orig_prior$mu.w[,2], "norm")  
      # plot(iprior.mu.w21)
      
      iprior.mu.w22 <- fitdist(post_Teunis_orig_prior$mu.w[,5], "norm")  
      # plot(iprior.mu.w22)
      
      iprior.mu.w31 <- fitdist(post_Teunis_orig_prior$mu.w[,3], "norm")  
      # plot(iprior.mu.w31)
      
      iprior.mu.w32 <- fitdist(post_Teunis_orig_prior$mu.w[,6], "norm")  
      # plot(iprior.mu.w32)
      
      iprior.w <- rbind(c(iprior.mu.w11$estimate[1], iprior.mu.w11$estimate[2],
                          iprior.mu.w12$estimate[1], iprior.mu.w12$estimate[2]),
                        c(iprior.mu.w21$estimate[1], iprior.mu.w21$estimate[2],
                          iprior.mu.w22$estimate[1], iprior.mu.w22$estimate[2]),
                        c(iprior.mu.w31$estimate[1], iprior.mu.w31$estimate[2],
                          iprior.mu.w32$estimate[1], iprior.mu.w32$estimate[2]))
    }
    
    # z ~ normal, need to be fixed
    {
      iprior.z <- rbind(c(mean(post_Teunis_orig_prior$prior.z[,1]), mean(post_Teunis_orig_prior$prior.z[,4])),
                        c(mean(post_Teunis_orig_prior$prior.z[,2]), mean(post_Teunis_orig_prior$prior.z[,5])),
                        c(mean(post_Teunis_orig_prior$prior.z[,3]), mean(post_Teunis_orig_prior$prior.z[,6])))

    }
    
    # logconc ~ normal
  {
    iprior.conc <- c(mean(post_Teunis_orig_prior$prior.conc[,1]), mean(post_Teunis_orig_prior$prior.conc[,2]))
    
  }
    
    
    # ltp0 ~ normal
  {
    iprior.ltp0 <- c(mean(post_Teunis_orig_prior$prior.p0[,1]), mean(post_Teunis_orig_prior$prior.p0[,2])) 
  }
  
    # Save the informed priors                    
    informed_prior_Teunis <- list("iprior.w" = iprior.w,
                                  "iprior.z" = iprior.z,
                                  "iprior.conc" = iprior.conc,
                                  "iprior.ltp0" = iprior.ltp0)
}
  # Save the lists to a RData file
  save(Teunis_coda, file="Teunis_coda.RData")
  save(informed_prior_Teunis, file="Teunis_iprior.RData")
  save(post_Teunis_orig_prior, file="Teunis_posterior_orig_prior.RData")
  
}
load("Teunis_coda.RData")  
load("Teunis_iprior.RData")# informed_prior_Teunis
load("Teunis_posterior_orig_prior.RData") # post_Teunis_orig_prior

##### CHRISTENSEN CPM MODEL ################################
if(FALSE)
  {
  # Run the model
  {
    CCPM_Christensen <- jags.model(file = "ChristensenCPM.R", # Christensen CPM
                                   data = Christensen_data_and_prior,
                                   n.chains = 4,
                                   n.adapt = 1000)
    
    update(CCPM_Christensen, 10^4) # burn-in
    
    parametertosaveCCPM <- c("prior.t.t_ce", "prior.t.t_es",
                             "prior.s.t_ce", "prior.s.t_es",
                             "t_cesim", "t_essim")
    
    CCPM_coda <- coda.samples(CCPM_Christensen,
                            parametertosaveCCPM,
                            10^4)
    sammatCCPM <- as.matrix(CCPM_coda)
  }
  
  # Save posterior samples
  {
    post_ChristensenCPM_orig_prior <- list(t_ce = sammatCCPM[,grep("t_cesim",colnames(sammatCCPM))],
      t_es = sammatCCPM[,grep("t_essim",colnames(sammatCCPM))],
      prior.t.t_ce = sammatCCPM[,grep("prior.t.t_ce",colnames(sammatCCPM))],
      prior.s.t_ce = sammatCCPM[,grep("prior.s.t_ce",colnames(sammatCCPM))],
      prior.t.t_es = sammatCCPM[,grep("prior.t.t_es",colnames(sammatCCPM))],
      prior.s.t_es = sammatCCPM[,grep("prior.s.t_es",colnames(sammatCCPM))])
  }
  
  # Fit distributions to marginal posterior toget informed priors
  {
    iprior.t.t_ce <- fitdist(post_ChristensenCPM_orig_prior$prior.t.t_ce, distr="beta")
    iprior.s.t_ce <- fitdist(post_ChristensenCPM_orig_prior$prior.s.t_ce, distr="gamma")
    iprior.t.t_es <- fitdist(post_ChristensenCPM_orig_prior$prior.t.t_es, distr="beta")
    iprior.s.t_es <- fitdist(post_ChristensenCPM_orig_prior$prior.s.t_es, distr="gamma")
  
    save_prior.t_ce <- rbind(c(iprior.t.t_ce$estimate[1], iprior.t.t_ce$estimate[2]),
                               c(iprior.s.t_ce$estimate[1], iprior.s.t_ce$estimate[2]))
    save_prior.t_es <- rbind(c(iprior.t.t_es$estimate[1], iprior.t.t_es$estimate[2]),
                               c(iprior.s.t_es$estimate[1], iprior.s.t_es$estimate[2]))
      
    informed_prior_ChristensenCPM <- list("prior.t_ce"=save_prior.t_ce,
                                            "prior.t_es"=save_prior.t_es)
  }
  # Save the lists to a RData file
  save(CCPM_coda, file="CCPM_coda.RData")
    save(informed_prior_ChristensenCPM, file="ChristensenCPM_iprior.RData")
    save(post_ChristensenCPM_orig_prior, file="ChristensenCPM_posterior_orig_prior.RData")
}
load("CCPM_coda.RData")
load("ChristensenCPM_iprior.RData") # informed_prior_ChristensenCPM
load("ChristensenCPM_posterior_orig_prior.RData") # post_ChristensenCPM_orig_prior

##### MYLIUS CPM ###########################################
if(FALSE)
  {
  # Run the model
  {
    CPM_Mylius <- jags.model(file = "MyliusCPM.R", # Mylius CPM
                             data = Mylius_data_and_prior,
                             n.chains = 4,
                             n.adapt = 1000)
    
    update(CPM_Mylius, 10^4) # burn-in
    
    parametertosaveMCPM <- c("prior.t.t_ch", "prior.s.t_ch",
                             "prior.t.t_cb", "prior.s.t_cb",
                             "prior.mu.t_bs","prior.precision.t_bs",
                             "t_ch_sim", "t_cb_sim",
                             "neg_log_t_bs_sim")
    
    MCPM_coda <- coda.samples(CPM_Mylius,
                            parametertosaveMCPM,
                            10^4)
    sammatMCPM <- as.matrix(MCPM_coda)
  }
  
  # Get the posterior samples
  {
    post_MyliusCPM_orig_prior <- list("t.t_ch" = sammatMCPM[,grep("prior.t.t_ch",colnames(sammatMCPM))],
                                           "s.t_ch" = sammatMCPM[,grep("prior.s.t_ch",colnames(sammatMCPM))],
                                           "t.t_cb"=sammatMCPM[,grep("prior.t.t_cb",colnames(sammatMCPM))],
                                           "s.t_cb"=sammatMCPM[,grep("prior.s.t_cb",colnames(sammatMCPM))],
                                           "mu.t_bs"=sammatMCPM[,grep("prior.mu.t_bs",colnames(sammatMCPM))],
                                           "precision.t_bs"=sammatMCPM[,grep("prior.precision.t_bs",colnames(sammatMCPM))],
                                            t_chmat = sammatMCPM[,grep("t_ch_sim",colnames(sammatMCPM))],
                                            t_cbmat = sammatMCPM[,grep("t_cb_sim",colnames(sammatMCPM))],
                                            neg_log_t_bs_mat = sammatMCPM[,grep("neg_log_t_bs_sim",colnames(sammatMCPM))])
    
  }
  
  # Fit distributions to marginal posteriors to get informed priors
  {
    iprior.t.t_ch <- fitdist(post_MyliusCPM_orig_prior$t.t_ch, distr="beta")
    iprior.s.t_ch <- fitdist(post_MyliusCPM_orig_prior$s.t_ch, distr="gamma")
    iprior.t.t_cb <- fitdist(post_MyliusCPM_orig_prior$t.t_cb, distr="beta")
    iprior.s.t_cb <- fitdist(post_MyliusCPM_orig_prior$s.t_cb, distr="gamma")
    iprior.mu.t_bs <- fitdist(post_MyliusCPM_orig_prior$mu.t_bs, distr="norm")
    iprior.precision.t_bs <- fitdist(post_MyliusCPM_orig_prior$precision.t_bs, distr="gamma")

  }

  # Save the data and informed hyperparameters
  {
    {
      save_prior.t_ch <- rbind(c(iprior.t.t_ch$estimate[1], iprior.t.t_ch$estimate[2]),
                               c(iprior.s.t_ch$estimate[1], iprior.s.t_ch$estimate[2]))
      save_prior.t_cb <- rbind(c(iprior.t.t_cb$estimate[1], iprior.t.t_cb$estimate[2]),
                               c(iprior.s.t_cb$estimate[1], iprior.s.t_cb$estimate[2]))
      save_prior.t_bs <- rbind(c(iprior.mu.t_bs$estimate[1], iprior.mu.t_bs$estimate[2]),
                               c(iprior.precision.t_bs$estimate[1], iprior.precision.t_bs$estimate[2]))
      
      informed_prior_MyliusCPM <- list("prior.t_ch"=save_prior.t_ch,
                                       "prior.t_cb"=save_prior.t_cb,
                                       "prior.neg_log_t_bs"=save_prior.t_bs)
    }
    
    # Save the lists to a RData file
    save(MCPM_coda, file="MCPM_coda.RData")
    save(informed_prior_MyliusCPM, file="MyliusCPM_iprior.RData")
    save(post_MyliusCPM_orig_prior, file="MyliusCPM_posterior_orig_prior.RData")
  }     
}
load("MCPM_coda.RData")
load("MyliusCPM_iprior.RData") # informed_prior_MyliusCPM
load("MyliusCPM_posterior_orig_prior.RData") # post_MyliusCPM_orig_prior
}

##################################################################
############# INDIVIDUAL MODELS WITH INFORMED PRIORS #############
##################################################################
{
source("FunctionsToUse.R")
# Different data sets along with the informed prior
alldataT <- get_data_for_ind_info_model_and_combined_model(prior.t_ec, prior.t_nodata)
# save(alldataT, file='alldataT.Rdata')
load("alldataT.Rdata")


##### MIKKELA'S PREVALENCE AND CONCENTRATION MODEL WITH INFORMED PRIOR ####
if(FALSE)
{
  ## Runs the model
  {
    Mikkela_jags_info <- jags.model(file = "MikkelaModelC_info.R",
                               data = alldata$Mikkela_data_and_informed_priors,
                               inits = alldata$initdatanew,
                               n.chains = 4,
                               n.adapt = 1000) # for adapting or tuning.
    
    update(Mikkela_jags_info, 10^4) # burn-in
    
    Mikkela_parametertosave <- c("pb", "pw", "pf", "pba","pfa", "pvarb", "pvarw",
                                 "vartot",  'taue', 'logitpb', 'e', 'mu_zero',
                                 'tauw','sdb')
    
    Mikkela_coda_info <- coda.samples(Mikkela_jags_info,
                                Mikkela_parametertosave,
                                10^4)
    
    Mikkela_mat_info <- as.matrix(Mikkela_coda_info)
    
    post_Mikkela_info_prior <- list(
      logitpb_1 = Mikkela_mat_info[,grep("logitpb",colnames(Mikkela_mat_info))][, 1], # logitpb is the first month
      taue = Mikkela_mat_info[,grep("taue",colnames(Mikkela_mat_info))],
      e_1 = Mikkela_mat_info[,grep("\\be\\b",colnames(Mikkela_mat_info))][, 1], # e_1 is the first month
      pw = Mikkela_mat_info[,grep("pw",colnames(Mikkela_mat_info))],
      tauw = Mikkela_mat_info[,grep("tauw",colnames(Mikkela_mat_info))],
      mu_zero = Mikkela_mat_info[,grep("mu_zero",colnames(Mikkela_mat_info))],
      sdb = Mikkela_mat_info[,grep("sdb",colnames(Mikkela_mat_info))]
    )
  }
  
  # Save the lists to a RData file
  save(Mikkela_coda_info, file="Mikkela_coda_info.RData")
  save(post_Mikkela_info_prior, file="Mikkela_posterior_info_prior.RData")
}
load("Mikkela_coda_info.RData")
load("Mikkela_posterior_info_prior.RData") #post_Mikkela_info_prior

##### TEUNIS DOSE-RESPONSE MODEL WITH INFORMED PRIOR #######
if(FALSE)
{
  ## Runs the model 
  {
    Teunis_jags_info <- jags.model(file = "TeunisModel_info.R",
                              data = alldata$Teunis_data_and_informed_priors,
                              n.chains = 4,
                              n.adapt = 1000)
    
    update(Teunis_jags_info, 10^4) # burn-in
    
    parametertosaveDR <- c('mu.w', 'w1', 'z1', 'w2', 'z2', 'logconc',
                           'ltp0', 
                           'a', 'b', 'eta', 'r',
                           'prill','prinf'
                           )
    
    Teunis_coda_info <- coda.samples(Teunis_jags_info,
                           parametertosaveDR,
                           10^4)
    
    sammatDR_info <- as.matrix(Teunis_coda_info)
    
    post_Teunis_info_prior <- list(mu.w = sammatDR_info[,grep('mu.w',colnames(sammatDR_info))],
                                  z1 = sammatDR_info[,grep('z1',colnames(sammatDR_info))],
                                  z2 = sammatDR_info[,grep('z2',colnames(sammatDR_info))],
                                  logconc = sammatDR_info[,grep('logconc',colnames(sammatDR_info))],
                                  ltp0 = sammatDR_info[,grep('ltp0',colnames(sammatDR_info))],
                                  a = sammatDR_info[,grep('\\ba\\b',colnames(sammatDR_info))],
                                  b = sammatDR_info[,grep('\\bb\\b',colnames(sammatDR_info))],
                                  eta = sammatDR_info[,grep('eta',colnames(sammatDR_info))], # e-14, too low. 
                                  r = sammatDR_info[,grep('\\br\\b',colnames(sammatDR_info))],
                                  prill = sammatDR_info[,grep('prill',colnames(sammatDR_info))],
                                  prinf = sammatDR_info[,grep('prinf',colnames(sammatDR_info))]
                                  ) # 58, too high)
    
  }
    # Save the list to a RData file
    save(Teunis_coda_info, file="Teunis_coda_info.RData")
    save(post_Teunis_info_prior, file="Teunis_posterior_info_prior.RData")
}
load("Teunis_coda_info.RData")
load("Teunis_posterior_info_prior.RData") # post_Teunis_info_prior

##### CHRISTENSEN CPM MODEL WITH INFORMED PRIOR ############
if(FALSE)
{
  # Run the model
  {
    CCPM_Christensen_info <- jags.model(file = "ChristensenCPM.R", # Christensen CPM
                                   data = alldata$Christensen_data_and_informed_priors,
                                   n.chains = 4,
                                   n.adapt = 1000)
    
    update(CCPM_Christensen_info, 10^4) # burn-in
    
    parametertosaveCCPM <- c("prior.t.t_ce", "prior.t.t_es",
                             "prior.s.t_ce", "prior.s.t_es",
                             "t_cesim", "t_essim")
    # parametertosaveCCPM <- c("t_ce", "t_es", "transneglogt_ce", "transneglogt_es")
    
    CCPM_coda_info <- coda.samples(CCPM_Christensen_info,
                            parametertosaveCCPM,
                            10^4)
    sammatCCPM_info <- as.matrix(CCPM_coda_info)
  }
  
  # Save posterior samples
  {
    # transneglogt_cemat <- sammatCCPM_info[,grep("transneglogt_ce",colnames(sammatCCPM_info))]
    # transneglogt_esmat <- sammatCCPM_info[,grep("transneglogt_es",colnames(sammatCCPM_info))]
    post_ChristensenCPM_info_prior <- list(t_ce = sammatCCPM_info[,grep("t_cesim",colnames(sammatCCPM_info))],
                                           t_es = sammatCCPM_info[,grep("t_essim",colnames(sammatCCPM_info))],
                                           prior.t.t_ce = sammatCCPM_info[,grep("prior.t.t_ce",colnames(sammatCCPM_info))],
                                           prior.s.t_ce = sammatCCPM_info[,grep("prior.s.t_ce",colnames(sammatCCPM_info))],
                                           prior.t.t_es = sammatCCPM_info[,grep("prior.t.t_es",colnames(sammatCCPM_info))],
                                           prior.s.t_es = sammatCCPM_info[,grep("prior.s.t_es",colnames(sammatCCPM_info))])
  }
  
  # Save the lists to a RData file
  save(CCPM_coda_info, file="CCPM_coda_info.RData")
  save(post_ChristensenCPM_info_prior, file="ChristensenCPM_posterior_info_prior.RData")
}    
load("CCPM_coda_info.RData")
load("ChristensenCPM_posterior_info_prior.RData") # post_ChristensenCPM_info_prior

##### MYLIUS CPM WITH INFORMED PRIOR #######################
if(FALSE)
{
  # Run the model
  {
    CPM_Mylius_info <- jags.model(file = "MyliusCPM.R", # Mylius CPM
                             data = alldata$Mylius_data_and_informed_priors,
                             n.chains = 4,
                             n.adapt = 1000)
    
    update(CPM_Mylius_info, 10^4) # burn-in
    
    parametertosaveMCPM <- c("prior.t.t_ch", "prior.s.t_ch",
                             "prior.t.t_cb", "prior.s.t_cb",
                             "prior.mu.t_bs","prior.precision.t_bs",
                             "t_ch_sim", "t_cb_sim",
                             "neg_log_t_bs_sim")
    
    MCPM_coda_info <- coda.samples(CPM_Mylius_info,
                            parametertosaveMCPM,
                            10^4)
    sammatMCPM_info <- as.matrix(MCPM_coda_info)
  }
  
  # Get the posterior samples
  {
    post_MyliusCPM_info_prior <- list("t.t_ch" = sammatMCPM_info[,grep("prior.t.t_ch",colnames(sammatMCPM_info))],
                                      "s.t_ch" = sammatMCPM_info[,grep("prior.s.t_ch",colnames(sammatMCPM_info))],
                                      "t.t_cb"=sammatMCPM_info[,grep("prior.t.t_cb",colnames(sammatMCPM_info))],
                                      "s.t_cb"=sammatMCPM_info[,grep("prior.s.t_cb",colnames(sammatMCPM_info))],
                                      "mu.t_bs"=sammatMCPM_info[,grep("prior.mu.t_bs",colnames(sammatMCPM_info))],
                                      "precision.t_bs"=sammatMCPM_info[,grep("prior.precision.t_bs",colnames(sammatMCPM_info))],
                                      t_chmat = sammatMCPM_info[,grep("t_ch_sim",colnames(sammatMCPM_info))],
                                      t_cbmat = sammatMCPM_info[,grep("t_cb_sim",colnames(sammatMCPM_info))],
                                      neg_log_t_bs_mat = sammatMCPM_info[,grep("neg_log_t_bs_sim",colnames(sammatMCPM_info))])
    
  }
    
  # Save the lists to a RData file
  save(MCPM_coda_info, file = "MCPM_coda_info.RData")
  save(post_MyliusCPM_info_prior, file="MyliusCPM_posterior_info_prior.RData")
}
load("MCPM_coda_info.RData")
load("MyliusCPM_posterior_info_prior.RData") # post_MyliusCPM_info_prior
}

##################################################################
############# COMBINED MODEL WITH INFORMED PRIOR #################
##################################################################
{
 # alldata <- get_data_for_ind_info_model_and_combined_model()
#### Combined model with Christensen CPM ###################
if(FALSE)
  {    
    # Run the model  
    { 
      # alldata$combi_Christensen$iprior.z <- rbind(c(0, 1), c(0,1), c(0,1))
      # C for Christensen
      combi_C <- jags.model(file = "PrevCCPMILLChicken.R", 
                          data = alldataT$combi_Christensen,
                          inits = alldataT$initdatanew,
                          n.chains = 4,
                          n.adapt = 1000)

      # update(combi_C, 10^4) # burn-in  
      update(combi_C, 10^4) # burn-in  
      
      resultattosaveC <- c('pw', 'taue', 'logitpb', 'e', 'mu_zero',
                           'tauw','sdb',
                           'prior.t.t_ce', 'prior.t.t_es',
                           'prior.s.t_ce', 'prior.s.t_es',
                           'mu.w', 'z1', 'z2', 'logconc', 'ltp0',
                           'a', 'b', 'eta', 'r', 
                           'dose_pred_mean',
                           'prill_pred_outbreak_tot', 'prinf_pred_human_tot',
                           'prill_pred_human_tot', 'prinf_pred_outbreak_tot',
                           'pf')
      
      Ccombo_coda <- coda.samples(combi_C,
                               resultattosaveC,
                               10^4)
      sammatC <- as.matrix(Ccombo_coda)
      # range(sammatC[,grep("dose_pred",colnames(sammatC))])
    }
    
    # Get the posterior 
    {
      post_combi_Christensen_info_prior <- list(
      logitpb_1 = sammatC[,grep("logitpb",colnames(sammatC))][, 1],
      taue = sammatC[,grep("taue",colnames(sammatC))],
      e_1 = sammatC[,grep("\\be\\b",colnames(sammatC))][, 1],
      pw = sammatC[,grep("pw",colnames(sammatC))],
      tauw = sammatC[,grep("tauw",colnames(sammatC))],
      mu_zero = sammatC[,grep("mu_zero",colnames(sammatC))],
      sdb = sammatC[,grep("sdb",colnames(sammatC))],
      prior.t.t_cemat = sammatC[,grep("prior.t.t_ce",colnames(sammatC))],
      prior.s.t_cemat = sammatC[,grep("prior.s.t_ce",colnames(sammatC))],
      prior.t.t_esmat = sammatC[,grep("prior.t.t_es",colnames(sammatC))],
      prior.s.t_esmat = sammatC[,grep("prior.s.t_es",colnames(sammatC))],
      mu.w = sammatC[,grep('mu.w',colnames(sammatC))],
      z1 = sammatC[,grep('z1',colnames(sammatC))],
      z2 = sammatC[,grep('z2',colnames(sammatC))],
      logconc = sammatC[,grep('logconc',colnames(sammatC))],
      ltp0 = sammatC[,grep('ltp0',colnames(sammatC))],
      
      a = sammatC[,grep('\\ba\\b',colnames(sammatC))], # 58, too high
      b = sammatC[,grep('\\bb\\b',colnames(sammatC))], # e-14, too low 
      eta = sammatC[,grep('eta',colnames(sammatC))], # e-14, too low. 
      r = sammatC[,grep('\\br\\b',colnames(sammatC))], # 58, too high
      
      prill_pred_human_tot = sammatC[,grep('prill_pred_human_tot',colnames(sammatC))],
      prill_pred_outbreak_tot = sammatC[,grep('prill_pred_outbreak_tot',colnames(sammatC))],
      prinf_pred_human_tot = sammatC[,grep('prinf_pred_human_tot',colnames(sammatC))],
      prinf_pred_outbreak_tot = sammatC[,grep('prinf_pred_outbreak_tot',colnames(sammatC))],
      # dose_if_cont_pred = sammatC[,grep('dose_if_cont_pred',colnames(sammatC))],
      dose_pred_mean = sammatC[,grep('dose_pred_mean',colnames(sammatC))],
      pf = sammatC[,grep('\\bpf\\b',colnames(sammatC))]
      # c_pred = sammatC[,grep('c_pred',colnames(sammatC))]
      ) 
    
    }
  
    # Save objects to RData files  
    save(Ccombo_coda, file="Ccombo_coda2.RData")
    save(sammatC, file="sammatC2.Rdata")
    save(post_combi_Christensen_info_prior, file="combi_Christensen_posterior_info_prior2.RData")
}
load("Ccombo_coda2.RData")
load("sammatC2.Rdata")  
load("combi_Christensen_posterior_info_prior2.RData")  
  
######### Combined model with Mylius CPM ####################
if(FALSE)
  {  
    # Run the model
    {
      # M for Mylius
      combi_M <- jags.model(file = "PrevMyliusILLChicken.R", 
                                      data = alldata$combi_Mylius,
                                      inits = alldata$initdatanew,
                                      n.chains = 4,
                                      n.adapt = 1000)
      
      update(combi_M, 10^4) # burn-in  
      
      resultattosaveM <- c('pw', 'taue', 'logitpb', 'e', 'mu_zero',
                           'tauw','sdb',
                           'prior.t.t_ch', 'prior.t.t_cb',
                           'prior.s.t_ch', 'prior.s.t_cb',
                           'prior.mu.t_bs', 'prior.precision.t_bs',
                           # "t_hs", "t_hh", "t_bb", "t_ss",
                           'mu.w', 'z1', 'z2', 'logconc', 'ltp0',
                           'a', 'b', 'eta', 'r',
                           'dose_pred_mean',
                           'prill_pred_outbreak_tot', 'prinf_pred_human_tot',
                           'prill_pred_human_tot', 'prinf_pred_outbreak_tot',
                           'pf'
                           )
     
      Mcombo_coda <- coda.samples(combi_M,
                               resultattosaveM,
                               10^4)
      
      sammatM <- as.matrix(Mcombo_coda)
    }
    
    # Get the posterior
    {
    post_combi_Mylius_info_prior <- list(
      logitpb_1 = sammatM[,grep("logitpb",colnames(sammatM))][, 1],
      taue = sammatM[,grep("taue",colnames(sammatM))],
      e_1 = sammatM[,grep("\\be\\b",colnames(sammatM))][, 1],
      pw = sammatM[,grep("pw",colnames(sammatM))],
      tauw = sammatM[,grep("tauw",colnames(sammatM))],
      mu_zero = sammatM[,grep("mu_zero",colnames(sammatM))],
      sdb = sammatM[,grep("sdb",colnames(sammatM))],
      prior.t.t_ch = sammatM[,grep("prior.t.t_ch",colnames(sammatM))],
      prior.s.t_ch = sammatM[,grep("prior.s.t_ch",colnames(sammatM))],
      prior.t.t_cb = sammatM[,grep("prior.t.t_cb",colnames(sammatM))],
      prior.s.t_cb = sammatM[,grep("prior.s.t_cb",colnames(sammatM))],
      prior.mu.t_bs = sammatM[,grep("prior.mu.t_bs",colnames(sammatM))],
      prior.precision.t_bs = sammatM[,grep("prior.precision.t_bs",colnames(sammatM))],
      mu.w = sammatM[,grep('mu.w',colnames(sammatM))],
      z1 = sammatM[,grep('z1',colnames(sammatM))],
      z2 = sammatM[,grep('z2',colnames(sammatM))],
      logconc =  sammatM[,grep('logconc',colnames(sammatM))],
      ltp0 = sammatM[,grep('ltp0',colnames(sammatM))],
     
      a = sammatM[,grep('\\ba\\b',colnames(sammatM))], # 58, too high
      b = sammatM[,grep('\\bb\\b',colnames(sammatM))], # e-14, too low 
      eta = sammatM[,grep('eta',colnames(sammatM))], # e-14, too low. 
      r = sammatM[,grep('\\br\\b',colnames(sammatM))], # 58, too high
     
      prill_pred_human_tot = sammatM[,grep('prill_pred_human_tot',colnames(sammatM))],
      prill_pred_outbreak_tot = sammatM[,grep('prill_pred_outbreak_tot',colnames(sammatM))],
      prinf_pred_human_tot = sammatM[,grep('prinf_pred_human_tot',colnames(sammatM))],
      prinf_pred_outbreak_tot = sammatM[,grep('prinf_pred_outbreak_tot',colnames(sammatM))],
      dose_pred_mean = sammatM[,grep('dose_pred_mean',colnames(sammatM))],
      pf = sammatM[,grep('\\bpf\\b',colnames(sammatM))]
    )
    }
  
    
    # Save objects to RData files
    save(Mcombo_coda, file="Mcombo_coda2.RData")
    save(sammatM, file="sammatM2.Rdata")
    save(post_combi_Mylius_info_prior, file="combi_Mylius_posterior_info_prior2.RData")
    
}
load("Mcombo_coda.RData")
load("sammatM.Rdata")
load("combi_Mylius_posterior_info_prior.RData")
}

##################################################################
############# GENERATE DIFFERENT FIGURES AND PLOTS ###############
##################################################################

source("FunctionsToUse.R")  

# Generate the Figure 4 and 5, Mikkela's article
Recreate_Mikkela_article_Fig_4_and_5(Mikkela_coda)

# Generate the dose-response curve in Figure 4, Teunis' article
Recreate_Teunis_article_Fig_4()

# Generate the corresponding Table A4 in the Appendix, Teunis' article
Generate_Teunis_table_A4()


# Generate the tables with informed priors
{
get_table_with_informed_prior_Mikkela()
get_table_with_informed_prior_Christensen()
get_table_with_informed_prior_Mylius()
get_table_with_informed_prior_Teunis()
}

# Plot all the posteriors for comparison
{
plot_the_posterior_Mikkela(post_Mikkela_orig_prior, 
                           post_Mikkela_info_prior,
                           post_combi_Christensen_info_prior, 
                           post_combi_Mylius_info_prior)

plot_the_posterior_Christensen(post_ChristensenCPM_info_prior, 
                               post_ChristensenCPM_info_prior,
                               post_combi_Christensen_info_prior)

plot_the_posterior_Mylius(post_MyliusCPM_orig_prior, 
                          post_MyliusCPM_info_prior,
                          post_combi_Mylius_info_prior)


plot_the_posterior_Teunis(post_Teunis_orig_prior,
                          post_Teunis_info_prior,
                          post_combi_Christensen_info_prior,
                          post_combi_Mylius_info_prior)
}

# Plot the pf, mean dose and prill
{
plot_pf_meandose_prill(post_combi_Christensen_info_prior, 
                       maxy = 1, marginalline =3, "Christensen")

plot_pf_meandose_prill(post_combi_Mylius_info_prior, 
                       maxy=0.1, marginalline=3, "Mylius")
}


##################################################################
####### ESTIMATION OF HOW MANY TIMES PEOPLE EAT CHICKEN ##########
##################################################################
{
##### Here the model is run under a scenario made up from Finnish health data  

# The model to estimate the number of times people eat chicken
{
  # for(m in 1:12){
ms ='model{
 consumers_reporters ~ dunif(minc,maxc)
 consumers_reporters1 <- round(consumers_reporters)
 for(i in 1:n.mcmc){
    for(j in 1:n.obs){
     y1[j,i] ~ dbin(prill_CPM1[i,month[j]],consumers_reporters1)
     y2[j,i] ~ dbin(prill_CPM2[i,month[j]],consumers_reporters1)
    }}
  
}'
}

# Input data
{
number_of_ill <- finnishdata$val #*10 # 7-11%
monthindex <- finnishdata$Month

p1 <- post_combi_Christensen_info_prior$prill_pred_human_tot[
  sample.int(nrow(post_combi_Christensen_info_prior$prill_pred_human_tot), 1000),]

p2 <- post_combi_Mylius_info_prior$prill_pred_human_tot[
  sample.int(nrow(post_combi_Mylius_info_prior$prill_pred_human_tot), 1000),]

p1[p1==0] <- 0.00001 # p can't be 0 in dbin
p2[p2==0] <- 0.00001

finska.data <- list(y1=matrix(rep(number_of_ill,nrow(p1)),ncol=nrow(p1)),
                    y2=matrix(rep(number_of_ill,nrow(p2)),ncol=nrow(p2)),
                    month=monthindex, n.obs = length(number_of_ill),
                    n.mcmc = nrow(p1),
                    prill_CPM1=p1,prill_CPM2=p2,
                    minc=100000, maxc=7000000)
}

# Run the model
{
jags <- jags.model(textConnection(ms),
                   data = finska.data,
                   n.chains = 2,
                   n.adapt = 1000)

update(jags, 10^4) # burn-in  

parametertosaveX <- c("consumers_reporters1") #, "y1_pred", "y2_pred")

ThedataX <- coda.samples(jags,
                         parametertosaveX,
                         10^4)

matrixX <- as.matrix(ThedataX)

boxplot(p1)
boxplot(p2)
save(matrixX, file="matrixX.RData")
}
  
n0 = 314700  
plot(colMeans(p1*n0),colMeans(p2*n0))
abline(0,1)
head(finska.data$y1)

# Make predictions about number of campylobacteriosis cases
{  
  ## new code
  {
    p1mean <- colMeans(p1)
    p2mean <- colMeans(p2)
    pred1 <- p1mean[finnishdata$Month]
    pred2 <- p2mean[finnishdata$Month]
    n1 <- mean(finnishdata$val/pred1) # crude estimate of the number eating chickens calibrated on the estimats from the first model, assuming its the same number of chicken eaters every month
    n2 <- mean(finnishdata$val/pred2) # same as above, but caliibrated baed on the second model
    
  }
  
load("matrixX.RData")
z1 <- matrix(, 1000, 12)
z2 <- matrix(, 1000, 12)
for(i in 1:1000){
  for(j in 1:12){
    z1[i,j] <- rnorm(1, n1*p1mean[j], sqrt(n1*p1mean[j]*(1-p1mean[j])))
    z2[i,j] <- rnorm(1, n2*p2mean[j], sqrt(n2*p2mean[j]*(1-p2mean[j])))
  
#    z1[i,j] <- rbinom(1, n1, p1mean[j])
#    z2[i,j] <- rbinom(1, n2, p2mean[j])
  }}
}
  
save(z1,z2,file='tillullrika.Rdata')
load("tillullrika.Rdata")
# Plot the estimate of the number of times people eat chicken
{
  plot(density(matrixX), xlab = "Number of times people eat chicken", 
       col="darkblue")
  mtext("", outer = TRUE, cex = 1.5)
} 
  
# Plot the predictions
plot_pred_case_study(z1, z2, finnishdata)    
    
}

##################################################################
####### SENSITIVITY ANALYSIS ########################## ##########
##################################################################
{
  # install.packages("BASS")
  library('BASS')
  library('ggplot2')
  
  # The factors and output of interest
  {
  Christensen_factor <- c('pw', 'taue', 'logitpb\\[1,1\\]', '\\be[1,1]\\b', 'mu_zero',
    'tauw','sdb',
    'prior.t.t_ce', 'prior.t.t_es',
    'prior.s.t_ce', 'prior.s.t_es',
    'mu.w', 'z1\\[9,2\\]', 'z2\\[9,2\\]', 'z1\\[9,3\\]', 'z2\\[9,3\\]', 'logconc\\[5\\]', 
  'ltp0\\[5\\]')
  
  Mylius_factor <- c('pw', 'taue', 'logitpb\\[1,1\\]', '\\be[1,1]\\b', 'mu_zero',
    'tauw','sdb',
    'prior.t.t_ch', 'prior.t.t_cb',
    'prior.s.t_ch', 'prior.s.t_cb',
    'prior.mu.t_bs', 'prior.precision.t_bs',
    'mu.w', 'z1\\[9,2\\]', 'z2\\[9,2\\]', 'z1\\[9,3\\]', 'z2\\[9,3\\]', 'logconc\\[5\\]',
    'ltp0\\[5\\]')
  
  
  SA_output <- c('dose_pred_mean',
                 'prinf_pred_outbreak_tot',
                 'prinf_pred_human_tot',
                 'prill_pred_outbreak_tot', 
                 'prill_pred_human_tot', 
                 'pf')
  
  # List with the factors and the corresponding data (mcmc)
  xxC <- sammatC[,sort(unique(unlist(lapply(1:length(Christensen_factor),function(i){
    grep(Christensen_factor[i],colnames(sammatC))}))))]
  # colnames(xxC)
  
  xxM <- sammatM[,sort(unique(unlist(lapply(1:length(Mylius_factor),function(i){
    grep(Mylius_factor[i],colnames(sammatM))}))))]
  
  y <- sammatC[,grep(SA_output[4],colnames(sammatC))][,7] # 7 for the july month
  }
  
  # Analysis for the combined model with Christensen CPM
  {
  k <- ncol(xxC)
  
  ## Fit a response surface (polynom regression) using y as response and X as predictor 
  mod<-bass(xxC,y)
  # plot(mod) ## Plot the residuals and see if everything looks ok
  
  ## Calculate the Sobols sensitivity indices
  sens<-sobol(mod)
  plot(sens,cex.axis=.5) 
  
  colnames(mod$xx.des)
  
  ## Plot the first picture of all sensitivities (interested in only the 
  ## first k main effects and the second is the total sensitivity. 
  
  ## Code to make it a more nice picture 
 
  ## Create data frame with all sensitivity values 
  ## Pick out the values and name for parameters
  Tind <- as.numeric(colnames(sens$T))
  dat <- do.call('rbind',lapply(1:ncol(sens$T),function(j){
   data.frame(first_order_sensitivity=sens$S[,j],total_sensitivity = sens$T[,j],
               param = rep(colnames(mod$xx.des)[Tind[j]],nrow(sens$S)))
   }))
  
  #dat$param <- as.factor(dat$param)
  
  ## Choose the parameters that have biggest influence/impact 
  first_mean <- aggregate(dat$first_order_sensitivity,by=list(dat$param),mean)
  total_mean <- aggregate(dat$total_sensitivity,by=list(dat$param),mean)
  bet1 <- order(first_mean$x,decreasing=TRUE)
  bet2 <- order(total_mean$x,decreasing=TRUE)

  # Numbers of interest
  NumbersOfInterest=10
  
  #dat_urval <- dat
  dat_urval1 <- dat[dat$param %in%  first_mean$Group.1[bet1[1:NumbersOfInterest]],]
  p1 <- ggplot(dat_urval1,aes(x=param,y=first_order_sensitivity,fill = param)) +
    geom_boxplot(width=0.15) +
    geom_jitter(shape=16, position=position_dodge(0.8), cex = 0.5) +
    labs(title='Sobol\'s First order sensitivity, Christensen',
         x="parameter", y = "proportion variance")
  
  print(p1)
  
  dat_urval2 <- dat[dat$param %in%  total_mean$Group.1[bet2[1:NumbersOfInterest]],]
  p2 <- ggplot(dat_urval2,aes(x=param,y=total_sensitivity,fill = param)) +
    geom_boxplot(width=0.15) +
    geom_jitter(shape=16, position=position_dodge(0.8), cex = 0.5) +
    labs(title='Sobol\'s total sensitivity, Christensen',x="parameter", y = "proportion variance")
  
  print(p2) # pdf
  }  
  
  

hist(first_mean$x)
hist(total_mean$x)

}
