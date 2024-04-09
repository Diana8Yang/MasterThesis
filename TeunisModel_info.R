
model{
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
  # ? for prediction?
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
}