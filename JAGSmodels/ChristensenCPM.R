##################### CPM, Christensen ############################


model{
  ## Data in Luber can be used to update the parameter t_ce and t_es of the CPM model
  prior.t.t_ce ~ dbeta(prior.t_ce[1, 1], prior.t_ce[1,2]);
  prior.t.t_es ~ dbeta(prior.t_es[1, 1], prior.t_es[1,2]);
  prior.s.t_ce ~ dgamma(prior.t_ce[2, 1], prior.t_ce[2,2]);
  prior.s.t_es ~ dgamma(prior.t_es[2, 1], prior.t_es[2,2]);
  
  for(i in 1:ntr){ # ntr= 11
      # t_ce[i] ~ dbeta(prior.shape1[1], prior.shape2[1]); # beta instead of 10^(-pert(..))
      # t_es[i] ~ dbeta(prior.shape1[2], prior.shape2[2]);
      t_ce[i] ~ dbeta(prior.t.t_ce*prior.s.t_ce, (1-prior.t.t_ce)*prior.s.t_ce); # beta instead of 10^(-pert(..))
      t_es[i] ~ dbeta(prior.t.t_es*prior.s.t_es, (1-prior.t.t_es)*prior.s.t_es);
      
      # transneglogt_ce[i] ~ dbeta(prior.shape1[1], prior.shape2[1]); # beta instead of pert
      # transneglogt_es[i] ~ dbeta(prior.shape1[2], prior.shape2[2]);
    }
  t_cesim ~ dbeta(prior.t.t_ce*prior.s.t_ce, (1-prior.t.t_ce)*prior.s.t_ce); # beta instead of 10^(-pert(..))
  t_essim ~ dbeta(prior.t.t_es*prior.s.t_es, (1-prior.t.t_es)*prior.s.t_es);
  
  # Need to find good prior distribution for beta parameters
  # for(i in 1:2){
  # prior.shape1[i] ~ dgamma(1, 3);
  # prior.shape2[i] ~ dgamma();
  # }
}

