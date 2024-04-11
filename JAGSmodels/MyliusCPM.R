##### Mylius CPM

data{
  for(i in 1:ntr){
  neg_log_t_bs[i] <- -log(t_bs[i])/log(10); # log10
  }
}
model{
  
  prior.t.t_ch ~ dbeta(prior.t_ch[1, 1], prior.t_ch[1,2]);
  prior.s.t_ch ~ dgamma(prior.t_ch[2, 1], prior.t_ch[2,2]);
  prior.t.t_cb ~ dbeta(prior.t_cb[1, 1], prior.t_cb[1,2]);
  prior.s.t_cb ~ dgamma(prior.t_cb[2, 1], prior.t_cb[2,2]);
  
  prior.mu.t_bs ~ dnorm(prior.neg_log_t_bs[1,1], prior.neg_log_t_bs[1,2]);
  prior.precision.t_bs ~ dgamma(prior.neg_log_t_bs[2,1], prior.neg_log_t_bs[2,2]);
      
  for(i in 1:ntr2){
    t_ch[i] ~ dbeta(prior.t.t_ch*prior.s.t_ch, (1-prior.t.t_ch)*prior.s.t_ch);
  }
  for(i in 1:ntr){ # ntr= 11
    t_cb[i] ~ dbeta(prior.t.t_cb*prior.s.t_cb, (1-prior.t.t_cb)*prior.s.t_cb);
    neg_log_t_bs[i] ~ dnorm(prior.mu.t_bs, prior.precision.t_bs)T(0,); # truncated at t_bs = 1
    # t_bs[i] <- 10^(-neg_log_t_bs[i]);
  }
  
  
  t_ch_sim ~ dbeta(prior.t.t_ch*prior.s.t_ch, (1-prior.t.t_ch)*prior.s.t_ch);
  t_cb_sim ~ dbeta(prior.t.t_cb*prior.s.t_cb, (1-prior.t.t_cb)*prior.s.t_cb);
  neg_log_t_bs_sim ~ dnorm(prior.mu.t_bs, prior.precision.t_bs)T(0,);
  # t_bs_sim <- 10^(-neg_log_t_bs_sim);
  
}

# hist(testpert)
# fitdist(testtrans, "beta")
# test3 <- rbeta(1000, 0.253, 340.37)
# hist(test3)
# ks.test(testtrans, test3)
