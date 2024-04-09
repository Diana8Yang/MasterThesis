######################## Data collection ########################

#### Mikkelä et al. data ############
{
  # Data from the article 
  {
  #B = number of sampled batches (chicken, turkey)
  #x = numer of contaminated samples, n = number of samples
  #indm = month index (1, ...,12), indt = meat type index (1,2)
  #logconcentration = logarithm(base-10) of the concentration in a contaminated sample
  #indj = batch index for the contaminated batches (1, ...,48), 
  #lower = lower limit for the concentration (log10(1/ws), ws = 25)
  #upper = upper limit for the concentration (log10(0.5))
  
  logconcentration = c(-0.3010300,-0.3010300,-0.3010300,1.0000000,
                       -0.3010300,0.00,0.2787536,0.3010300,-0.3010300,-0.3010300,-0.3010300,0.3010300,0.5440680,
                       0.8450980,
                       -0.3010300,0.1760913,-0.3010300,0,1,0,-0.3010300,0,0,-0.3010300,-0.3010300,-0.3010300,0,
                       0.8129134,0.9030900,0.5440680,0.6334685,1.3944517,1.209515,1.580925,
                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.3010300,-0.3010300,
                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  
  is.censored <- ifelse(is.na(logconcentration), 1, 2) # t[which(c)] = NA
  # logconcentration_noNA <- ifelse(is.na(logconcentration),(log10(0.5)-log10(0.04))/2, logconcentration) # 0.548455
  
  
  Mikkeladata <- list(B = c(226,185),
                      x = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,1,0,0,0,0,0,0,0,
                            0,0,0,0,5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,
                            2,1,5,1,1,1,2,3,5,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,1,1,1,2,2,2,4,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,
                            2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0),
                      n = c(2,3,2,6,5,2,1,1,5,1,2,1,5,1,1,1,2,2,1,2,1,1,2,5,1,4,4,1,2,1,1,10,3,2,3,2,5,2,3,2,2,7,
                            1,1,1,1,2,2,1,1,6,1,2,1,2,2,1,2,2,1,1,1,1,1,3,1,3,1,2,2,4,1,1,3,1,2,2,4,1,2,1,2,5,1,1,1,3,
                            1,2,1,2,4,2,2,5,5,1,1,1,1,1,3,3,1,5,5,1,1,7,5,6,5,5,5,2,3,1,4,5,5,3,1,7,2,2,3,5,7,2,1,4,3,
                            5,8,6,3,2,1,1,3,3,7,3,5,2,3,5,3,5,6,3,6,9,1,5,2,8,2,2,2,4,4,5,2,3,5,2,2,1,1,1,5,4,1,2,3,4,
                            1,3,1,1,1,2,1,2,1,1,5,5,2,7,1,7,1,1,1,5,3,1,4,5,4,8,2,5,5,2,1,1,2,1,1,1,1,1,1,4,2,2,2,3,1,
                            1,1,1,1,2,1,5,3,1,4,1,3,1,8,1,2,1,4,1,2,5,4,2,5,4,1,1,3,2,4,1,2,4,2,9,1,2,1,1,1,2,1,1,3,1,
                            1,5,5,1,3,1,6,3,1,1,8,2,3,3,3,1,2,2,2,3,2,5,7,3,16,2,2,1,1,4,3,2,3,3,1,2,1,1,4,5,2,5,3,10,
                            1,4,7,10,5,2,3,3,5,6,2,2,2,11,5,2,4,1,2,3,2,6,2,1,1,1,3,3,3,1,7,5,2,3,8,3,4,4,1,5,3,9,4,2,
                            4,5,1,3,3,2,1,1,3,7,3,2,2,5,1,3,2,5,1,7,2,5,2,3,2,5,2,1,1,5,1,1,1,1,5,5,4,2,2,1,1,2,3,1,1,
                            3,4,3,2,1,1,1,3,2,7,1),
                      indm = c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,
                               4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,
                               6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                               7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,
                               9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,
                               11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,1,1,1,1,1,1,1,1,1,2,2,2,2,
                               2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,
                               6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,
                               8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
                               10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                               11,11,11,11,12,12,12,12,12,12)
                      ,indt = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                                2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                      logconcentration = c(-0.3010300,-0.3010300,-0.3010300,1.0000000,
                                           -0.3010300,0.00,0.2787536,0.3010300,-0.3010300,-0.3010300,-0.3010300,0.3010300,0.5440680,
                                           0.8450980,
                                           -0.3010300,0.1760913,-0.3010300,0,1,0,-0.3010300,0,0,-0.3010300,-0.3010300,-0.3010300,0,
                                           0.8129134,0.9030900,0.5440680,0.6334685,1.3944517,1.209515,1.580925,
                                           NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                           NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.3010300,-0.3010300,
                                           NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                      indj = c(1,1,1,1,3,3,3,3,7,9,9,9,9,9,10,10,12,13,17,19,21,21,21,23,23,24,24,27,
                               27,29,29,29,31,31,1,2,3,4,5,6,6,7,8,10,11,14,14,14,15,15,15,15,16,16,18,18,19,19,19,19,20,
                               20,22,23,23,23,23,24,24,24,25,26,26,28,29,30,40,40,32,33,33,34,35,35,35,35,36,37,38,38,38,
                               38,38,39,41,41,41,42,43,44,44,45,45,46,46,46,47,48),lower = -1.39794,upper = -0.30103,
                      is.censored = is.censored)
                      # is.censored=c(rep(0, 34), rep(1, 42), rep(0, 2), rep(1, 30)))

  # Inits for both chicken and turkey
  # initdata <- list(taue = c(0.5,0.5),
  #                  e = structure(.Data = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),.Dim = c(2,12)),
  #                  logitpb = structure(.Data = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),.Dim = c(2,12)),
  #                  tauw = c(0.01,0.01),mu_zero = c(0.1,0.1),sdb = c(1,1))
  # 
  initdata.func <- function(){list(taue = c(0.5,0.5),
                                   e = structure(.Data = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),.Dim = c(2,12)),
                                   logitpb = structure(.Data = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),.Dim = c(2,12)),
                                   tauw = c(0.01,0.01),mu_zero = c(0.1,0.1),sdb = c(1,1))}
  }

  # Values of the priors 
  {
  prior.pw <- rbind(c(0, 1), c(0, 1))
  prior.taue <- rbind(c(0.001, 0.001), c(0.001, 0.001))
  prior.logitpb <- rbind(c(0, 0.001), c(0, 0.001))
  prior.e <- rbind(c(0, 0.001), c(0, 0.001))
  
  prior.mu_zero <- rbind(c(0, 0.001), c(0, 0.001))
  prior.tauw <- rbind(c(0.001, 0.001), c(0.001, 0.001))
  prior.sdb <- rbind(c(0, 100), c(0, 100))
  
  priordata <- list("prior.pw" = prior.pw, "prior.taue" = prior.taue,
                    "prior.logitpb" = prior.logitpb, "prior.e" = prior.e,
                    "prior.mu_zero" = prior.mu_zero, "prior.tauw" = prior.tauw,
                    "prior.sdb" = prior.sdb)
  }
  
  # Data for chicken only 
  {
  B_c = 226
  x_c <- Mikkeladata$x[1:226]
  n_c <- Mikkeladata$n[1:226]
  indm_c <- Mikkeladata$indm[1:226]
  indt_c <- Mikkeladata$indt[1:226]
  logconcentration_c <- Mikkeladata$logconcentration[1:76]
  indj_c <- Mikkeladata$indj[1:76]
  is.censored_c <- ifelse(is.na(logconcentration_c), 1, 2)
  
  Mikkela_data_c <- list("B"=B_c, "x"=x_c, "n"=n_c, "indm"=indm_c, "indt"=indt_c, 
                         "logconcentration" = logconcentration_c, "indj"=indj_c, 
                         "lower" = -1.39794,"upper" = -0.30103, "is.censored"=is.censored_c)
  
  # Inits for chicken only
  initdata_c <- function(){list(taue = c(0.5),
                     e = structure(.Data = c(2,2,2,2,2,2,2,2,2,2,2,2),.Dim = c(1,12)),
                     logitpb = structure(.Data = c(1,1,1,1,1,1,1,1,1,1,1,1),.Dim = c(1,12)),
                     tauw = c(0.01),mu_zero = c(0.1),sdb = c(1))}

  
  # All data for chicken only
  Mikkela_all_data_c <- c(Mikkela_data_c, priordata)
  }
}

#### CPM and Luber et al. data ################
{
  # Christensen CPM, Luber's data
  {
  t_ce = c(0.9, 1.1, 1.6, 1.9, 0.2, 2.3, 0.9, 0.4, 1.2, 0.2, 1.4)/100
  t_es = c(0.0001, 9.1, 5.0, 10.5, 14.3, 7.7, 33.3, 20.0, 6.7, 0.0001, 6.7)/100
  # mlogt_ce = -log10(c(0.9, 1.1, 1.6, 1.9, 0.2, 2.3, 0.9, 0.4, 1.2, 0.2, 1.4)/100),
  # mlogt_es = -log10(c(0.0, 9.1, 5.0, 10.5, 14.3, 7.7, 33.3, 20.0, 6.7, 0.0, 6.7)/100),

  # Transfrom from Pert to Beta distribution, see wiki
  # prior.shape1 <- c(1.8, 1.8) 
  # prior.shape2 <- c(4.2, 4.2) 

  t1 = mean(t_ce) # 0.011
  t2 = mean(t_es) # 0.1030002
  s = length(t_ce) # 11
  
  hyper.t_ce <- rbind(c(1,1), c(1,1)) # beta(1,1), gamma(1,1)
  hyper.t_es <- rbind(c(1,1), c(1,1))
  
  Christensen_data <- list(t_ce = t_ce, t_es = t_es, ntr = s)
  Christensen_prior <- list(prior.t_ce=hyper.t_ce, prior.t_es=hyper.t_es)
  
  Christensen_data_and_prior <- c(Christensen_data, Christensen_prior)
  }
  # Additional data to be used in combined model
  {
    prior.t_ec <- list(prior.t_ec=rbind(c(1,1), # for t, mean
                                        c(1,1))) # for s, sample size
  }
  # Mylius CPM, Luber's data
  {
  legs_to_hands <- c(19.7, 3.1, 0.5, 2.1, 4.5, 0.2, 0.3, 0.5, 0.2, 0.6, 1.5)
  filets_to_hands <- c(1.2, 1.0, 1.6, 0.2, 0.3, 3.6, 0.6, 7.8, 2.5, 2.6, 0.9)
  both_to_hands <- c(23.9, 2.5, 5.3, 2.7, 1.5)
  
  t_ch <- c(legs_to_hands, filets_to_hands, both_to_hands)/100
  t_cb <- t_ce
  t_bs <- t_es
  
  ntr2 <- length(t_ch)
  
  prior.t_ch <- rbind(c(1, 1), c(1, 1))
  prior.t_cb <- rbind(c(1, 1), c(1, 1))
  prior.neg_log_t_bs <- rbind(c(0, 1), c(1, 1))
  
  Mylius_data <- list(t_ch=t_ch, t_cb=t_cb, t_bs=t_bs, ntr=s, ntr2=ntr2)
  Mylius_prior <- list(prior.t_ch=prior.t_ch, prior.t_cb=prior.t_cb, 
                       prior.neg_log_t_bs=prior.neg_log_t_bs)
  
  Mylius_data_and_prior <- c(Mylius_data, Mylius_prior)
  
  # # for the combined model with both Christensen and Mylius. Double ntr o/w.
  # Mylius_data_CM <- list(t_ch=t_ch, t_cb=t_cb, t_bs=t_bs, ntr2=ntr2)
  }
  # Additional data to be used in combined model
  {
    prior.t_nodata <- list(prior.neg_log_t_hs = c(1.90, 1/0.606^2),
                           prior.t_hh = c(0.24, 6.67),
                           prior.t_bb = c(0.25, 400),
                           prior.t_ss = c(3.25, 4.7))
  }
}

#### Teunis et al. data ###############
{
  # R code to collect literature data and prepare R objects for transfer to JAGS.
  
  ## Human volunteer data
  {
    # Black et al. (1988) Experimental Campylobacter jejuni infection in humans
    # Journal of Infectious Diseases 157(3):472-279
    human.1 <- list(strain="A3249",host="human",
                    dose=c(8e2,8e3,9e4,8e5,1e6,1e8,1e8), # final dose in think bicarbonate
                    expos=c(10,10,13,11,19,5,4),
                    infec=c(5,6,11,8,15,5,4),
                    sympt=c(1,1,6,1,2,0,2));
    
    human.2 <- list(strain="81-176",host="human",
                    dose=c(1e6,1e8,2e9),
                    expos=c(7,10,22),
                    infec=c(7,10,22),
                    sympt=c(3,6,9));
    
    # Tribble et al. (2009) Campylobacter jejuni strain CG8421: a refined model
    # for the study of campylobacteriosis and evaluation of Campylobacter
    # vaccines in human subjects
    # Clinical Infectious Diseases 49(10):1512-1519
    human.3 <- list(strain="CG8421",host="human",
                    dose=c(0.97e6,0.84e5,0.54e5),
                    expos=c(8,7,8),
                    infec=c(8,7,8),
                    sympt=c(8,6,8));
    
    # Tribble et al. (2010) Assessment of the duration of protection in
    # Campylobacter jejuni experimental infection in humans
    # Infection and Immunity 78(4):1750-1759
    human.4 <- list(strain="81-176",host="human",
                    dose=c(1e5,1e7,1e9,1e9,1e9), # in bicarbonate buffer
                    expos=c(5,5,36,8,7), # last two rechallenges (1-2months;1year)
                    infec=c(5,5,36,6,7),
                    sympt=c(3,2,33,0,4));
    
    # Kirkpatrick et al. (2013) Lack of homologous protection against
    # Campylobacter jejuni CG8421 in a human challenge model
    # Clinical Infectious Diseases 57(8):1106-1113
    ## same as Tribble et al. 2009
    human.5 <- list(strain="CG8421",host="human",
                    dose=c(5e5,5e5), # second dose rechallenge
                    expos=c(15,8),
                    infec=c(15,8),
                    sympt=c(14,8));
    
    # Unpublished data supplied by David Tribble (April 2017)
    human.6 <- list(strain="CG8421",host="human",
                    dose=c(1.9e4,2.0e4,1.4e5),
                    expos=c(6,7,4),
                    infec=c(6,7,4),
                    sympt=c(2,5,2));
    
    # Unpublished data supplied by David Tribble (April 2017)
    human.7 <- list(strain="CG8421",host="human",
                    dose=c(1.8e5),
                    expos=c(13),
                    infec=c(13),
                    sympt=c(11));
  }
  ## Human outbreak data
  {
    # Teunis et al. (2005) A reconsideration of the Campylobacter dose-response
    # relation
    # Epidemiology and Infection 133(4):583-592
    outbreak.1 <- list(strain=NA,host="human",
                       dose=c(0,1/6,0.5,1,2)*0.188, # contaminated milk (liters)
                       expos=c(35,12,18,21,6),
                       infec=c(NA,NA,NA,NA,NA),
                       sympt=c(2,2,7,13,6));
    
    # Evans et al. (1996) A milk-borne Campylobacter outbreak following 
    # an educational farm visit
    # Epidemiol. Infect. 117: (457-462)
    outbreak.2 <- list(strain=NA,host="human",
                       dose=c(0,0.5,1,2)*0.188, # contaminated milk (liters)
                       expos=c(17,7,21,5),
                       infec=c(NA,NA,NA,NA),
                       sympt=c(2,3,14,4));
    
    # Korlath et al. (1985) A point-source outbreak of campylobacteriosis
    # associated with consumption of raw milk
    # Journal of Infectious Diseases 152(3):592-596
    outbreak.3 <- list(strain="81-176",host="human",
                       dose=c(0,1,2)*0.188, # contaminated milk (liters); 2 or more
                       expos=c(20,20,30),
                       infec=c(NA,NA,NA),
                       sympt=c(0,7,18));
    
    # Blaser et al. (1987) The influence of immunity on raw milk–associated
    # Campylobacter infection
    # Journal of the American Medical Association 257(1):43-46
    outbreak.4 <- list(strain=NA,host="human",
                       dose=c(0,1,2,3)*0.188, # contaminated milk (liters); 3 or more
                       expos=c(2,11,6,8),
                       infec=c(NA,8,6,8), # 0 - < NA (pinf*pill|inf = p0)
                       sympt=c(0,6,5,8));
  }
  ## Primate challenge data
  {
    # Russell et al. (1989) Experimental Campylobacter jejuni infection in
    # Macaca nemestrina
    # Infection and Immunity 57(5):1438-1444
    primate.1 <- list(strain="81-176",host="Macaca nemestrina",
                      dose=c(3e11),
                      expos=c(4),
                      infec=c(4),
                      sympt=c(4));
    
    # Russell et al. (1993) Early colonic damage and invasion of Campylobacter
    # jejuni in experimentally challenge infant Macaca mulatta
    # Journal of Infectious Diseases 168(1):210-215
    primate.2 <- list(strain="78-37",host="Macaca mulatta",
                      dose=c(2.7e10),
                      expos=c(1),
                      infec=c(1),
                      sympt=c(1));
    
    # Fitzgeorge et al. (1981) Experimental infection of Rhesus monkeys with
    # a human strain of Campylobacter jejuni
    # Journal of Hygiene 86(3):343-351
    primate.3 <- list(strain="v212x",host="Rhesus",
                      dose=c(1e9,1e10,1e9), # last dose is a rechallenge
                      expos=c(6,2,2), # of the animals who got 1e10
                      infec=c(6,2,2),
                      sympt=c(5,1,0));
    
    # Jones et al. (2006) New world monkey Aotus nancymae as a model for
    # Campylobacter jejuni infection and immunity
    # Infection and Immunity 74(1):790-793
    primate.4 <- list(strain="81-176",host="Aotus nancymae",
                      dose=c(8e8,6e10,5e12),
                      expos=c(6,6,6),
                      infec=c(6,6,6),
                      sympt=c(1,4,5));
    
    # Islam et al. (2006) Establishment of a non-human primate Campylobacter
    # disease model for the pre-clinical evaluation of Campylobacter vaccine
    # formulations
    # Vaccine 24(18):3762-3771
    primate.5 <- list(strain="81-176",host="Macaca mulatta",
                      dose=c(1e7,1e9,1e11),
                      expos=c(10,10,10),
                      infec=c(10,10,10),
                      sympt=c(2,8,7));
    
  }
  # Organize the data  
  {  
    hstnam <- c("primate","human","outbreak");
    stnum <- c(1,2,3,4,5, 1,2,3,4,5,6,7, 1,2,3,4); # different studies
    obnum <- c(rep(NA,5), rep(NA,7), 1,2,3,4); # ?
    # number of [] in each study:
    strains <- c(1,2,3,1,1, 4,1,5,1,5,5,5, 6,7,1,8); 
    hosts <- c(1,1,1,1,1, 2,2,2,2,2,2,2, 3,3,3,3); # hstnam coded
    n.dose <- c(1,1,2,3,3, 6,3,3,3,1,3,1, 5,4,3,4); 
    
    strain <- c(); host <- c(); dose <- c();
    expos <- c(); infec <- c(); sympt <- c();
    stnam <- c(); obn <- c();
    
    # Create a big collected obs. data with all the data above
    for(k in 1:length(n.dose)){ # for each observation
      strain <- c(strain,rep(strains[k],n.dose[k])); # what kind of strain in each obs
      stn <- eval(parse(text=paste(hstnam[hosts[k]],".",stnum[k],sep="")))$strain;
      obn <- c(obn,rep(obnum[k],n.dose[k])); # ?
      host <- c(host,rep(hosts[k],n.dose[k])); # what kind of host in each obs
      
      # collect all data, dose, expos, intec, sympt
      dse <- eval(parse(text=paste(hstnam[hosts[k]],".",stnum[k],sep="")))$dose;
      exs <- eval(parse(text=paste(hstnam[hosts[k]],".",stnum[k],sep="")))$expos;
      inf <- eval(parse(text=paste(hstnam[hosts[k]],".",stnum[k],sep="")))$infec;
      smp <- eval(parse(text=paste(hstnam[hosts[k]],".",stnum[k],sep="")))$sympt;
      
      stnam <- c(stnam,stn); # strain in each studies
      dose <- c(dose,dse[1:n.dose[k]]); # how much dose in each obs
      expos <- c(expos,exs[1:n.dose[k]]); # how many that got exposed in each obs
      infec <- c(infec,inf[1:n.dose[k]]); # how many that got infected in each obs
      sympt <- c(sympt,smp[1:n.dose[k]]); # how many that got symptoms in each obs
    }
    
    # Add a single subject without observed data (for prediction)
    # to human challenge studies (2) and to human outbreaks (3)
    strains <- c(strains,c(max(strains)+1,max(strains)+1));
    strain <- c(strain, c(max(strains),max(strains)));
    hosts <- c(hosts,c(2,3)); host <- c(host,c(2,3));
    n.dose <- c(n.dose,c(1,1)); dose <- c(dose,c(1e1,1e1));
    expos <- c(expos,c(1,1)); infec <- c(infec,c(NA,NA));
    sympt <- c(sympt,c(NA,NA)); obn <- c(obn,c(NA,max(obn,na.rm=TRUE)+1));
    # put the 2 new obs into right places
    byhost <- order(host);
    strain <- strain[byhost]; host <- host[byhost];
    dose <- dose[byhost]; expos <- expos[byhost];
    infec <- infec[byhost]; sympt <- sympt[byhost];
    obn <- obn[byhost]; n.ob <- max(obn,na.rm=TRUE); # number of outbreaks
    n.ch <- length(which(host <=2)); # 31, only primates and human plus one prediction
    last <- length(obn); # 48
  }
  
  # Prior
  {
  prior.w <- rbind(c(0,0.01,0,0.1),c(0,0.01,0,0.1),c(0,0.01,0,0.1));
  prior.z <- rbind(c(0.0001,1),c(0.0001,1),c(0.0001,1));
  tau.w <- c(1,1,1);
  
  prior.conc <- c(log(152),1/(4^2));
  prior.p0 <- c(-2,0.1);
  }
  
  drdata <- list("prior.w"=prior.w,"tau.w"=tau.w, # "prior.z"=prior.z,
                 # "prior.conc"=prior.conc,"prior.p0"=prior.p0,
                 "n.ch"=n.ch,"n.ob"=n.ob,"obn"=obn,"last"=last,
                 "n.strain"=max(strain),"strain"=strain,
                 "n.host"=max(host),"host"=host,
                 "dose"=dose,"expos"=expos,"infec"=infec,"sympt"=sympt);
  
  drdata_graph <- list("prior.w"=prior.w,"tau.w"=tau.w, "prior.z"=prior.z,
                      "prior.conc"=prior.conc,"prior.p0"=prior.p0,
                      "n.ch"=n.ch,"n.ob"=n.ob,"obn"=obn,"last"=last,
                      "n.strain"=max(strain),"strain"=strain,
                      "n.host"=max(host),"host"=host,
                      "dose"=dose,"expos"=expos,"infec"=infec,"sympt"=sympt);
}





#### Additional data to use when running the whole model, to be fixed ##########

#### Illness data and additional priors for the whole model #######

# # Illness data in Finland https://sampo.thl.fi/pivot/prod/sv/ttr/shp/fact_shp?row=area-12260.&column=time-12048.12060.12207.&column=time-12235.12348.12383.12185.12430.12038.12029.12136.12045.12130.12140.12058.12063.12198.12366.12392.&filter=reportgroup-12162#

finnishdata <- read.csv("fact_shp.csv", header=T, sep=",") 

# # fix the data 
# 
# # https://population.un.org/wpp/Download/Standard/Population/
# N_exposed <- c(5412980,5436616,5459717) # 2012, 2013, 2014 
# N_exposed_mean <- mean(N_exposed) # 5436438
# N_exposed_sd <- sd(N_exposed) # 23369.01

# OBS use the factor of 10-11, see Sundstrom

# additional_data <- list(symp_pred = F_ill$
#   N_exposed_mean = N_exposed_mean, N_exposed_sd = N_exposed_sd)
# 
