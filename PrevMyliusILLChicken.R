# Prevalence, Christensen CPM and dose-response (ill) model 
data{
  for(i in 1:ntr){
    neg_log_t_bs[i] <- -log(t_bs[i])/log(10); # log10
  }
}
model{
  
  ########### Prevalence and concentration model, Mikkelä #################
  
  # Prevalence model
  # i = batch, m = month, t = meat type (1 = chicken, 2 = turkey)
  for(i in 1:B[1]){ #+B[2]){
    x[i] ~ dbin(pp[i],n[i]); # model for the positive samples in a batch i
    pp[i] <- I[i]*pw[indt[i]]; # probability of sampling positive sample
    I[i] ~ dbern(pb[indt[i],indm[i]]); # model for the "true contamination status" of a batch i
  }
  
  # Markovian time series for the between-batch prevalence
  for(t in 1){ #:2){
    for(m in 2:12){
      logitpb[t,m] ~ dnorm(logitpb[t,m-1],e[t,m]);
      e[t,m]~ dnorm(0, taue[t]); # model for the difference between consecutive months
    }}
  
  for(t in 1){ #:2){
    for(m in 1:12){
      pb[t,m] <- exp(logitpb[t,m])/(1+exp(logitpb[t,m])); # between-batch prevalence
    }}
  
  for(t in 1){ #:2){
    # pw[t] ~ dunif(0,1); # prior for the within-batch prevalence
    # taue[t] ~ dgamma(0.001,0.001); # prior for the precision (1/variance)
    # logitpb[t,1] ~ dnorm(0,0.001);
    # e[t,1] ~ dnorm(0,0.001);
    
    # pw[t] ~ dunif(iprior.pw[t, 1], iprior.pw[t, 2]); # prior for the within-batch prevalence
    pw[t] ~ dbeta(iprior.pw[t, 1], iprior.pw[t, 2]);
    # taue[t] ~ dgamma(iprior.taue[t, 1], iprior.taue[t, 2]); # prior for the precision (1/variance)
    taue[t] ~ dlnorm(iprior.taue[t, 1], 1/iprior.taue[t, 2]^2);
    logitpb[t,1] ~ dnorm(iprior.logitpb_1[t, 1], 1/iprior.logitpb_1[t, 2]^2);
    e[t,1] ~ dnorm(iprior.e_1[t, 1], 1/iprior.e_1[t, 2]^2);
    
  }
  
  # Predictions
  for(t in 1){ #:2){
    for(m in 1:12){
      pf[t,m] <- pw[t]*pb[t,m]; # monthly retail prevalence
    }
    pba[t] <- sum(pb[t,])/12; # average monthly between-batch prevalence
    pfa[t] <- sum(pf[t,])/12; # average monthly retail prevalence
  }
  
  # Concentration model
  
  # j = contaminated batch, k = contaminated sample, t = meat type (1 = chicken, 2 = turkey)
  # hierarchical log-normal model for the positive concentrations
  for(k in 1:34){# exact observations chicken
    logconcentration[k] ~ dnorm(mu[indj[k]],tauw[1]);
  }
  
  for(k in 35:76){# censored observations chicken
    is.censored[k] ~ dinterval(logconcentration[k], c(lower,upper));
    logconcentration[k] ~ dnorm(mu[indj[k]],tauw[1]);
    # logconcentration[k] ~ dnorm(mu[indj[k]],tauw[1])C(lower,upper); #in bugs
  }
  
  ## Batch
  for(j in 1:31){
    mu[j] ~ dnorm(mu_zero[1],taub[1]); # mean concentration in chicken batch j
  }
  
  for(t in 1){ #:2){
    # mu_zero[t] ~ dnorm(0,0.001); # prior for the mean concentration in all contaminated batches
    # tauw[t] ~ dgamma(0.001,0.001); # prior for the within-batch precision (1/variance)
    # sdb[t] ~ dunif(0,100); # prior for the between-batch standard deviation
    # taub[t] <- pow(sdb[t],-2); # between-batch precision (1/variance)
    mu_zero[t] ~ dnorm(iprior.mu_zero[t, 1], 1/iprior.mu_zero[t, 2]^2); # prior for the mean concentration in all contaminated batches
    tauw[t] ~ dgamma(iprior.tauw[t, 1], 1/iprior.tauw[t, 2]); # prior for the within-batch precision (1/variance)
    sdb[t] ~ dgamma(iprior.sdb[t, 1], iprior.sdb[t, 2]); # prior for the between-batch standard deviation
    taub[t] <- pow(sdb[t],-2); # between-batch precision (1/variance)
  }
  
  # Predictions
  for(t in 1){ #:2){
    murep[t] ~ dnorm(mu_zero[t],taub[t]); # predicted mean concentration in a random 
    # contaminated batch
    log10crep[t] ~ dnorm(murep[t],tauw[t]); # predicted log10-concentration in a random 
    # contaminated retail meat unit
    crep[t] <- pow(10,log10crep[t]); # predicted concentration in a random contaminated 
    # retail meat unit
    varb[t] <- 1/taub[t]; # between-batch variance
    varw[t]<- 1/tauw[t]; # within-batch variance
    vartot[t] <- varb[t]+varw[t]; # total variance
    pvarb[t] <- varb[t]/vartot[t]; # proportion of the total variance explained by 
    # between-batch variance
    pvarw[t]<- varw[t]/vartot[t]; # proportion of the total variance explained by 
    # within-batch variance
  }
  
  
  ########### CPM, Mylius ############################################
  # {
  # Parameters in the CPM model
  
  
  ## Assumed transfer rates
  # from hand to salad
  neg_log_t_hs ~ dnorm(prior.neg_log_t_hs[1], prior.neg_log_t_hs[2])T(0,);
  t_hs <- 10^(-neg_log_t_hs); 
  
  # hand washing
  t_hh_distr1 <- 1; # no washing
  t_hh_distr2 ~ dbeta(prior.t_hh[1], prior.t_hh[2]); # washing
  t_hh_prob ~ dunif(0, 1); 
  indicator_hh <- ifelse(t_hh_prob <=0.2, 1, 0);
  t_hh <- (t_hh_distr1 * indicator_hh) + (t_hh_distr2 * (1 - indicator_hh));
  
  # board washing
  t_bb_distr1 <- 0;
  t_bb_distr2 <- 1;
  t_bb_distr3 ~ dbeta(prior.t_bb[1], prior.t_bb[2]);
  t_bb_prob ~ dunif(0, 1);
  indicator_bb1 <- ifelse(t_bb_prob <= 0.33, 1, 0);
  indicator_bb2 <- ifelse(0.33 < t_bb_prob && t_bb_prob <= (0.33+0.05), 1, 0);
  indicator_bb3 <- ifelse(t_bb_prob >(0.33+0.05) , 1, 0);
  t_bb <- (t_bb_distr1 * indicator_bb1) + (t_bb_distr2 * indicator_bb2) + (t_bb_distr3*indicator_bb3);
  
  # salad washing
  t_ss_distr1 <- 1;
  neg_log_t_ss_distr2 ~ dbeta(prior.t_ss[1], prior.t_ss[2]);
  t_ss_distr2 <- 10^(-neg_log_t_ss_distr2);
  t_ss_prob ~ dunif(0, 1);
  indicator_ss <- ifelse(t_ss_prob <=0.4, 1, 0);
  t_ss <- (t_ss_distr1 * indicator_ss) + (t_ss_distr2 * (1 - indicator_ss));
  
  # prior on parameters for t_ch, t_cb and t_bs
  prior.t.t_ch ~ dbeta(prior.t_ch[1, 1], prior.t_ch[1,2]);
  prior.s.t_ch ~ dgamma(prior.t_ch[2, 1], prior.t_ch[2,2]);
  prior.t.t_cb ~ dbeta(prior.t_cb[1, 1], prior.t_cb[1,2]);
  prior.s.t_cb ~ dgamma(prior.t_cb[2, 1], prior.t_cb[2,2]);
  
  prior.mu.t_bs ~ dnorm(prior.neg_log_t_bs[1,1], 1/prior.neg_log_t_bs[1,2]^2);
  prior.precision.t_bs ~ dgamma(prior.neg_log_t_bs[2,1], prior.neg_log_t_bs[2,2]);
  
  for(i in 1:ntr2){
    # from chicken to hand
    t_ch[i] ~ dbeta(prior.t.t_ch*prior.s.t_ch, (1-prior.t.t_ch)*prior.s.t_ch);
  }
  for(i in 1:ntr){ # ntr= 11
    # from chicken to board
    t_cb[i] ~ dbeta(prior.t.t_cb*prior.s.t_cb, (1-prior.t.t_cb)*prior.s.t_cb);
    # from board to salad
    neg_log_t_bs[i] ~ dnorm(prior.mu.t_bs, prior.precision.t_bs)T(0,); # truncated at t_bs = 1
    # t_bs[i] <- 10^(-log_t_bs[i]);
  }
  
  ## Draw values for transfer rates
  t_ch_pred ~ dbeta(prior.t.t_ch*prior.s.t_ch, (1-prior.t.t_ch)*prior.s.t_ch);
  t_cb_pred ~ dbeta(prior.t.t_cb*prior.s.t_cb, (1-prior.t.t_cb)*prior.s.t_cb);
  neg_log_t_bs_pred ~ dnorm(prior.mu.t_bs, prior.precision.t_bs)T(0,);
  t_bs_pred <- 10^(-neg_log_t_bs_pred);

  # Probability of at least one colony-forming unit in a portion of meet
  p_tr_pred <- (t_ch_pred*t_hh*t_hs+t_cb_pred*t_bb*t_bs_pred)*t_ss;
  
  # #  Prediction for each of month where we have prevalence and illness statistics
  # {
  # exposure for NN persons chicken eaters
  # chicken_eaters <- 100;
  # assumed number of chicken eaters per month
  for(m in 1:12){
    for(p in 1:N_person){
      ## Assumed variability in portion sizes
      w_c[p, m] ~ dlnorm(189, 1/127^2)T(0, 1000); #(5, 1/0.8^2)T(0, 1000);  # Log portion sizes, in gram.
      #w_c[p, m] <- exp(logw_c[p, m]); # Portion sizes, in gram
      #w_c[p, m] = w_c[p, m][w_c[p, m]<= 1000];
      mu_pred[p, m] ~ dnorm(mu_zero[1],taub[1]); # predicted mean concentration in a random if there is something ther
      # contaminated batch
      log10c_pred[p, m] ~ dnorm(mu_pred[p, m],tauw[1]); # predicted log10-concentration in a random
      # contaminated retail meat unit
      c_pred[p, m] <- pow(10,log10c_pred[p, m]); # predicted concentration in a random contaminated
      # the portion to end up in the dose
      lambda_pred[p, m] <- c_pred[p, m]*w_c[p, m] # lambda from Poisson model for N_portion
      Nportion[p, m] ~ dpois(lambda_pred[p, m]); # Number of Campylobacters (cfu) on one portion of consumed chicken meat
      dose_if_cont_pred[p, m] ~ dbin(p_tr_pred, Nportion[p, m]) # dose in Teunis code, cV in article (in cfu).
      cont[p, m] ~ dbin(pf[1,m],1) ## if the person eat chicken from a cont batch
      dose_pred[p, m] <- dose_if_cont_pred[p, m]*cont[p, m]
      ## ger dos för en slumpmässig kycklingätare person i finland per månad
    }
    dose_pred_mean[m] <- mean(dose_pred[,m]);  
  }

  # ########### Dose-response model, Teunis #################################
  # { 
  # parameters/variables which doesnt depend on dose
  
  # parent nodes: hyperparameters
  mu.w[1,1] ~ dnorm(iprior.w[1,1],1/iprior.w[1,2]^2);
  mu.w[1,2] ~ dnorm(iprior.w[1,3],1/iprior.w[1,4]^2);
  mu.w[2,1] ~ dnorm(iprior.w[2,1],1/iprior.w[2,2]^2);
  mu.w[2,2] ~ dnorm(iprior.w[2,3],1/iprior.w[2,4]^2);
  mu.w[3,1] ~ dnorm(iprior.w[3,1],1/iprior.w[3,2]^2);
  mu.w[3,2] ~ dnorm(iprior.w[3,3],1/iprior.w[3,4]^2);
  for(k.st in 1:n.strain){
    # infection
    for(k.hs in 1:n.host){
      w1[k.st,k.hs] ~ dnorm(mu.w[1,1] + mu.w[1,2],tau.w[1]);
      # z1[k.st,k.hs] ~ dnorm(iprior.z1[k.st, k.hs*2-1],1/iprior.z1[k.st,k.hs*2]^2);
      z1[k.st,k.hs] ~ dnorm(iprior.z[1,1], iprior.z[1,2]);
      u1[k.st,k.hs] <- exp(w1[k.st,k.hs])/(1+exp(w1[k.st,k.hs]));
      v1[k.st,k.hs] <- exp(z1[k.st,k.hs]);
      a[k.st,k.hs] <- u1[k.st,k.hs]*v1[k.st,k.hs];
      b[k.st,k.hs] <- (1-u1[k.st,k.hs])*v1[k.st,k.hs];
    }
    # symptoms
    for(k.hs in 1:(n.host-1)){
      w2[k.st,k.hs] ~ dnorm(mu.w[2,1] + mu.w[2,2],tau.w[2]);
      # z2[k.st,k.hs] ~ dnorm(iprior.z2[k.st,k.hs*2-1],1/iprior.z2[k.st,k.hs*2]^2);
      z2[k.st,k.hs] ~ dnorm(iprior.z[2,1], iprior.z[2,2]);
      u2[k.st,k.hs] <- exp(w2[k.st,k.hs])/(1+exp(w2[k.st,k.hs]));
      v2[k.st,k.hs] <- exp(z2[k.st,k.hs]);
      eta[k.st,k.hs] <- u2[k.st,k.hs]*v2[k.st,k.hs];
      r[k.st,k.hs] <- (1-u2[k.st,k.hs])*v2[k.st,k.hs];
    }
    # For outbreak
    w2[k.st,n.host] ~ dnorm(mu.w[3,1] + mu.w[3,2],tau.w[3]);
    # z2[k.st,n.host] ~ dnorm(iprior.z2[k.st,n.host*2-1],1/iprior.z2[k.st,n.host*2]^2);
    z2[k.st,n.host] ~ dnorm(iprior.z[3,1], iprior.z[3,2]);
    u2[k.st,n.host] <- exp(w2[k.st,n.host])/(1+exp(w2[k.st,n.host]));
    v2[k.st,n.host] <- exp(z2[k.st,n.host]);
    eta[k.st,n.host] <- u2[k.st,n.host]*v2[k.st,n.host];
    r[k.st,n.host] <- (1-u2[k.st,n.host])*v2[k.st,n.host];
  }
  
  # Likelihood for Teunis model
  
  for(k in 1:n.ch){ # primates and human + one prediction (k = 31)
    # infection
    num[k] ~ dpois(dose[k]);
    gamma[k] <- loggam(a[strain[k],host[k]]+b[strain[k],host[k]]) -
      loggam(a[strain[k],host[k]]+b[strain[k],host[k]]+num[k]) +
      loggam(b[strain[k],host[k]]+num[k]) -
      loggam(b[strain[k],host[k]]);
    prinf[k] <- (1-exp(gamma[k]));
    infec[k] ~ dbin(prinf[k],expos[k]);
    # symptoms
    prill[k] <- (1-pow(1+(dose[k]/eta[strain[k],host[k]]),
                       -r[strain[k],host[k]]));
    sympt[k] ~ dbin(prill[k],infec[k]);
  }
  for(k.ob in 1:n.ob){ # to be replaced with C_ret from prevalence model?
    logconc[k.ob] ~ dnorm(iprior.conc[1],1/iprior.conc[2]^2);
    conc[k.ob] <- exp(logconc[k.ob]);
    ltp0[k.ob] ~ dnorm(iprior.ltp0[1],1/iprior.ltp0[2]^2);
    p0[k.ob] <- exp(ltp0[k.ob])/(1+exp(ltp0[k.ob]));
  }
  for(k in (n.ch+1):last){ # outbreaks only plus one prediction (k=48)
    num[k] ~ dpois(conc[obn[k]]*dose[k]);
    gamma[k] <- loggam(a[strain[k],host[k]]+b[strain[k],host[k]]) -
      loggam(a[strain[k],host[k]]+b[strain[k],host[k]]+num[k]) +
      loggam(b[strain[k],host[k]]+num[k]) -
      loggam(b[strain[k],host[k]]);
    prinf[k] <- ifelse(dose[k]==0,1,1-exp(gamma[k]));
    infec[k] ~ dbin(prinf[k],expos[k]);
    # symptoms
    prill[k] <- ifelse(dose[k]==0,p0[obn[k]],1-(1-p0[obn[k]])*(1-prinf[k]*
                      (1-pow(1+(conc[obn[k]]*dose[k]/eta[strain[k],host[k]]),
                      -r[strain[k],host[k]]))));
    # sympt[k] ~ dbin(1-(1-p0[obn[k]])*(1-prinf[k]*prill[k]),expos[k]);
    sympt[k] ~ dbin(prill[k],expos[k]);
  }
  
  # # Prediction using the Teunis model for predicted dose, one person prediction only
  # {
  k_human <- n.ch; # Prediction for human challenge, k=n.ch= 31, or prediction for outbreak, k = 48
  k_outbreak <- last; # prediction for outbreak, k=last= 48
  
  for(m in 1:12){
    for(p in 1:N_person){ # simulate 100 persons
    
      # infections
      num_pred_human[p, m] ~ dpois(dose_pred[p, m]);
      gamma_pred_human[p, m] <- loggam(a[strain[k_human],host[k_human]]+b[strain[k_human],host[k_human]]) -
        loggam(a[strain[k_human],host[k_human]]+b[strain[k_human],host[k_human]]+num_pred_human[p, m]) +
        loggam(b[strain[k_human],host[k_human]]+num_pred_human[p, m]) -
        loggam(b[strain[k_human],host[k_human]]);
      prinf_pred_human[p, m] <- ifelse((1-exp(gamma_pred_human[p, m]))<0,0,(1-exp(gamma_pred_human[p, m])));
      #prinf_pred_human[p, m] <-1-exp(gamma_pred_human[p, m]);
      infec_pred_human[p, m] ~ dbin(prinf_pred_human[p, m],expos[k_human]); # expos[k_human] = 1
      # symptoms
      prill_pred_human[p, m] <- (1-pow(1+(dose_pred[p, m]/eta[strain[k_human],host[k_human]]),
                                       -r[strain[k_human],host[k_human]]));
      sympt_pred_human[p, m] ~ dbin(prill_pred_human[p, m],infec_pred_human[p, m]);
      
      # k=n.ob=5, for prediction of the outbreak 
      ltp0_pred_outbreak[p, m] ~ dnorm(iprior.ltp0[1],1/iprior.ltp0[2]^2);
      p0_pred_outbreak[p, m] <- exp(ltp0_pred_outbreak[p, m])/(1+exp(ltp0_pred_outbreak[p, m]));
      
      num_pred_outbreak[p, m] ~ dpois(dose_pred[p, m]);
      gamma_pred_outbreak[p, m] <- loggam(a[strain[k_outbreak],host[k_outbreak]]+b[strain[k_outbreak],host[k_outbreak]]) -
        loggam(a[strain[k_outbreak],host[k_outbreak]]+b[strain[k_outbreak],host[k_outbreak]]+num_pred_outbreak[p, m]) +
        loggam(b[strain[k_outbreak],host[k_outbreak]]+num_pred_outbreak[p, m]) -
        loggam(b[strain[k_outbreak],host[k_outbreak]]);
      
      # prinf_pred_outbreak[p, m] <- 1-exp(gamma_pred_outbreak[p, m]);
      prinf_pred_outbreak[p, m] <- ifelse(1-exp(gamma_pred_outbreak[p, m])<0, 0, 1-exp(gamma_pred_outbreak[p, m]));
      infec_pred_outbreak[p, m] ~ dbin(prinf_pred_outbreak[p, m],expos[k_outbreak]);
      # symptoms
      prill_pred_outbreak[p, m] <- ifelse(dose_pred[p, m]==0,p0_pred_outbreak[p, m], 1-(1-p0_pred_outbreak[p, m])*(1-prinf_pred_outbreak[p, m]*
                                          (1-pow(1+(dose_pred[p, m]/eta[strain[k_outbreak],
                                          host[k_outbreak]]), -r[strain[k_outbreak],host[k_outbreak]]))));
      # prill_pred_outbreak[p, m] <- 1-(1-prinf_pred_outbreak[p, m]*
      #                                   (1-pow(1+(dose_pred[p, m]/eta[strain[k_outbreak],
      #                                   host[k_outbreak]]), -r[strain[k_outbreak],host[k_outbreak]])));
      sympt_pred_outbreak[p, m] ~ dbin(prill_pred_outbreak[p, m],expos[k_outbreak]);
      
    }
    ## For a every month in a year
    prinf_pred_human_tot[m] <- mean(infec_pred_human[,m]);
    prill_pred_human_tot[m] <- mean(sympt_pred_human[,m]);
    
    prinf_pred_outbreak_tot[m] <- mean(infec_pred_outbreak[,m]);
    prill_pred_outbreak_tot[m] <- mean(sympt_pred_outbreak[,m]);
    
  }

}    