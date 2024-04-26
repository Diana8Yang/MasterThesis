

## Functions to load data for the models with informed priors
get_data_for_ind_info_model_and_combined_model <- function(C_additional_data, M_additional_data){
  load("Mikkela_iprior.RData")
  load("ChristensenCPM_iprior.RData")
  load("MyliusCPM_iprior.RData")
  load("Teunis_iprior.RData")
  
  load("Mikkela_posterior_orig_prior.RData")
  load("ChristensenCPM_posterior_orig_prior.RData")
  load("MyliusCPM_posterior_orig_prior.RData")
  load("Teunis_posterior_orig_prior.RData")
  
  # Replace the priors values in the input data with informed priors
  
  # Mikkela  
  Mikkela_data_and_informed_priors <- c(Mikkeladata, informed_prior_Mikkela)
  
  initdatanew <- list(taue = c(0.5),
                      e = structure(.Data = c(2,2,2,2,2,2,2,2,2,2,2,2),.Dim = c(1,12)),
                      logitpb = structure(.Data = c(1,1,1,1,1,1,1,1,1,1,1,1),.Dim = c(1,12)),
                      tauw = c(0.01),mu_zero = c(0.1),sdb = c(1)) #c(1/100))
  
  Christensen_data_and_informed_priors <- c(Christensen_data, informed_prior_ChristensenCPM)
  
  Mylius_data_and_informed_priors <- c(Mylius_data, informed_prior_MyliusCPM)
  
  # Teunis    
  Teunis_data <- list("tau.w"=tau.w, "n.ch"=n.ch,"n.ob"=n.ob,"obn"=obn,"last"=last,
                      "n.strain"=max(strain),"strain"=strain,
                      "n.host"=max(host),"host"=host,
                      "dose"=dose,"expos"=expos,"infec"=infec,"sympt"=sympt);
  
  Teunis_data_and_informed_priors <- c(Teunis_data, informed_prior_Teunis)
  
  alldata_CPMChristensen1000000 <- c(Mikkela_data_and_informed_priors, Christensen_data_and_informed_priors,
                              C_additional_data, Teunis_data_and_informed_priors, N_person = 1000000)
  
  alldata_CPMMylius1000000 <- c(Mikkela_data_and_informed_priors, Mylius_data_and_informed_priors,
                         M_additional_data, Teunis_data_and_informed_priors, N_person = 1000000) #200)
  
  alldataT <- vector('list', 7)
  alldataT[[1]] <- Mikkela_data_and_informed_priors
  alldataT[[2]] <- Teunis_data_and_informed_priors
  alldataT[[3]] <- Christensen_data_and_informed_priors
  alldataT[[4]] <- Mylius_data_and_informed_priors
  alldataT[[5]] <- alldata_CPMChristensen1000000 #alldata_CPMChristensen
  alldataT[[6]] <- alldata_CPMMylius1000000 #alldata_CPMMylius
  alldataT[[7]] <- initdatanew
  
  names(alldataT) <- c('Mikkela_data_and_informed_priors', 
                      'Teunis_data_and_informed_priors',
                      'Christensen_data_and_informed_priors',
                      'Mylius_data_and_informed_priors',
                      'combi_Christensen','combi_Mylius', 'initdatanew')
  
  return(alldataT)    
}    

# alldata <- get_data_for_ind_info_model_and_combined_model(prior.t_ec, prior.t_nodata)


############ Generate result for report #################

#### Recreation of result from articles ####
  
# Recreate graphs from Mikkela article, Fig 4 and 5 (chicken only)
Recreate_Mikkela_article_Fig_4_and_5 <- function(Mikkela_coda){
  # summary(Mikkela_sam)[1,] # to double-check with the result from the article
  Mikkela_mat <- as.matrix(Mikkela_coda)
  pfa <- Mikkela_mat[,grep("pfa",colnames(Mikkela_mat))] # chicken, has to be changed
  plot(density(pfa))
  
  pb_m <- Mikkela_mat[,grep("pb",colnames(Mikkela_mat))][, 13:24]
  colnames(pb_m) <- month.abb
  #plot(density(pf_season))
  boxplot(pb_m*100, outline = F,
          ylab = "Between-batch prevalence (%)", ylim=c(0, 100))
  
  pf_m <- Mikkela_mat[,grep("pf",colnames(Mikkela_mat))][, 1:12]
  colnames(pf_m) <- month.abb
  #plot(density(pf_season))
  boxplot(pf_m*100, outline = F,
          ylab = "Prevalence (%)", ylim=c(0, 100))
  
  
}
  
# To double-check the parameters with the Teunis article
# Run the Teunis model (altered code) in order to recreate the dose-response curve
Recreate_Teunis_article_Fig_4 <- function(){
  {
  Topredict_dose = c(10^(-2), 10^(-1), 10^(-0), 10^(1), 10^(2),
                     10^(3), 10^(4), 10^(5), 10^(6), 10^(7), 
                     10^(8), 10^(9), 10^(10), 10^(11), 10^(12))
  
  Predict_dose <-list(Topredict_dose = Topredict_dose, 
                      NpredictDose = length(Topredict_dose))
  Test_Teunis_data <- c(drdata_graph, Predict_dose)
  
  jags <- jags.model(file = "TeunisModelGraphs.R",
                     data = Test_Teunis_data,
                     n.chains = 4,
                     n.adapt = 10^4)
  
  update(jags, 1000) # burn-in
  
  parametertosaveDRgraph <- c( 'eta', 'r', 'a', 'b',
                               "Topredict_dose", "prill_pred_human", "prill_pred_outbreak",
                               "prinf_pred_human", "prinf_pred_outbreak")
  
  sam2DRgraph <- coda.samples(jags,
                              parametertosaveDRgraph,
                              10^4)
  
  sammatDRG <- as.matrix(sam2DRgraph)
  
  post_Teunis_orig_prior_graph <- list(a = sammatDRG[,grep('\\ba\\b', colnames(sammatDRG))],
                                       b = sammatDRG[,grep('\\bb\\b',colnames(sammatDRG))],
                                       # v2 = sammatDRG[,grep('v2',colnames(sammatDRG))], # 54, too high
                                       # u2 = sammatDRG[,grep('u2',colnames(sammatDRG))], # e-12, too low 
                                       eta = sammatDRG[,grep('eta',colnames(sammatDRG))], # e-12, too low. 
                                       r = sammatDRG[,grep('\\br\\b',colnames(sammatDRG))], # 54, too high))
                                       prinf_pred_human = sammatDRG[,grep('prinf_pred_human',colnames(sammatDRG))],
                                       prinf_pred_outbreak = sammatDRG[,grep('prinf_pred_outbreak',colnames(sammatDRG))],
                                       prill_pred_human = sammatDRG[,grep('prill_pred_human',colnames(sammatDRG))],
                                       prill_pred_outbreak = sammatDRG[,grep('prill_pred_outbreak',colnames(sammatDRG))])
  
  
   save(post_Teunis_orig_prior_graph, file="TableA4.Rdata")
   
  TheSgraph <-summary(sam2DRgraph)  # Error in ar.yw.default(x, aic = aic, order.max = order.max, na.action = na.action,  :
  TheDosegraph <- TheSgraph$statistics[grep('Topredict_dose',rownames(TheSgraph$statistics)),1]
}
  
  # Recreate Teunis dose-response curve for infection from the article, Fig 2a
  {
    par(mfrow=c(1,1))
    plot(TheDosegraph, TheSgraph$statistics[grep('prinf_pred_human',rownames(TheSgraph$statistics)),1], type='l', log = 'x', col="blue", ylab="P(inf)", xlab='Dose', ylim=c(0,1))
    lines(TheDosegraph, TheSgraph$quantiles[grep('prinf_pred_human',rownames(TheSgraph$quantiles)),'2.5%'], col="blue", lty=2)
    lines(TheDosegraph, TheSgraph$quantiles[grep('prinf_pred_human',rownames(TheSgraph$quantiles)),'97.5%'], col="blue", lty=2)
    
    # plot(TheDosegraph, TheSgraph$statistics[grep('prinf_pred_outbreak',rownames(TheSgraph$statistics)),1], type='l', log = 'x', col="blue", ylab="P(inf)", xlab='Dose', ylim=c(0,1))
    # lines(TheDosegraph, TheSgraph$quantiles[grep('prinf_pred_outbreak',rownames(TheSgraph$quantiles)),'2.5%'], col="blue", lty=2)
    # lines(TheDosegraph, TheSgraph$quantiles[grep('prinf_pred_outbreak',rownames(TheSgraph$quantiles)),'97.5%'], col="blue", lty=2)
  }
  
  # Recreate Teunis dose-response curve for illness from the article, Fig 4
  {   
    # par(mfrow=c(1,2))
    plot(TheDosegraph, TheSgraph$statistics[grep('prill_pred_human',rownames(TheSgraph$statistics)),1], type='l', log = 'x', col="blue", ylab="P(ill)", xlab='Dose', main="Human",ylim=c(0,1))
    lines(TheDosegraph, TheSgraph$quantiles[grep('prill_pred_human',rownames(TheSgraph$quantiles)),'2.5%'], col="blue", lty=2)
    lines(TheDosegraph, TheSgraph$quantiles[grep('prill_pred_human',rownames(TheSgraph$quantiles)),'97.5%'], col="blue", lty=2)
    
    plot(TheDosegraph, TheSgraph$statistics[grep('prill_pred_outbreak',rownames(TheSgraph$statistics)),1], type='l', log = 'x', col="blue", ylab="P(ill)", xlab='Dose', main="Outbreak",ylim=c(0,1))
    lines(TheDosegraph, TheSgraph$quantiles[grep('prill_pred_outbreak',rownames(TheSgraph$quantiles)),'2.5%'], type='l', col="blue", lty=2)
    lines(TheDosegraph, TheSgraph$quantiles[grep('prill_pred_outbreak',rownames(TheSgraph$quantiles)),'97.5%'], type='l', col="blue", lty=2)

}  
  
}

Generate_Teunis_table_A4 <- function(){
  load("TableA4.Rdata")
  library(xtable)
  
  abrn <- post_Teunis_orig_prior_graph[c("a", "b", "r", "eta")]
  abrn_matrix <- cbind(abrn$a, abrn$b, abrn$r, abrn$eta)
  abrn_matrixtest <- abrn_matrix[, -grep(paste0("\\[4,1]|\\[5,1]|\\[6,1]|\\[7,1]|\\[8,1]|\\[9,1]|",
                                                "\\[2,2]|\\[3,2]|\\[6,2]|\\[7,2]|\\[8,2]|",
                                                "\\[2,3]|\\[3,3]|\\[4,3]|\\[5,3]"),
                                         colnames(abrn_matrix))]
  
  abrn_summary <- apply(abrn_matrixtest, 2, quantile, c(0.5, 0.025, 0.975))
  abrn_t <- t(abrn_summary)
  
  abrn_row <- cbind(abrn_t[grep("\\ba\\b", rownames(abrn_t)), ], 
                    abrn_t[grep("\\bb\\b", rownames(abrn_t)), ],
                    abrn_t[grep("\\br\\b", rownames(abrn_t)), ],
                    abrn_t[grep("\\beta\\b", rownames(abrn_t)), ])
  
  ab_final <- abrn_row[c(1:6, 9, 10, 8, 11, 7, 12), c(1:6)]
  rn_final <- abrn_row[c(1:6, 9, 10, 8, 11, 7, 12), c(7:12)]
  # options(digits = 7, scipen = 0) # default: 7, 0
  
  # xtable(ab_final, type = "latex", digits = 2, 
  #        display=c("s","fg","fg","f","fg","fg","f"))
  # xtable(rn_final, type = "latex", digits = 2, 
  #        display=c("s","fg","G","f","fg","G","f"))
  
  print(xtable(ab_final, type = "latex", digits = 2, 
               display=c("s","fg","fg","f","fg","fg","f")), file = "TeunisTableA4ab.tex")
  print(xtable(rn_final, type = "latex", digits = 2, 
               display=c("s","fg","g","f","fg","g","f")), file = "TeunisTableA4rn.tex")
  
}


#### Check the convergence
if(FALSE){
  {
# Mikkela 
{
  # traceplot(Mikkela_sam)
  gelman.diag(Mikkela_coda) # Error in chol.default(W) : the leading minor of order 51 is not positive definite
  gelman.plot(Mikkela_coda) # Not convergence in the parameter e, above 2
}
# Teunis
{
  # traceplot(sam2DR)
  gelman.diag(Teunis_coda) # Error in chol.default(W) : the leading minor of order 106 is not positive definite
  gelman.plot(Teunis_coda)
}
# Christensen
{
  # traceplot(samCCPM)
  gelman.diag(CCPM_coda)
  gelman.plot(CCPM_coda)
}
# Mylius
{
  # traceplot(samMCPM)
  gelman.diag(MCPM_coda)
  gelman.plot(MCPM_coda)
}
    
# Mikkela with informed prior
{  
  # traceplot(Mikkela_coda_info)
  gelman.diag(Mikkela_coda_info)
  gelman.plot(Mikkela_coda_info)
}
# Teunis with informed prior
{
  # traceplot(Teunis_coda_info)
  gelman.diag(Teunis_coda_info)
  gelman.plot(Teunis_coda_info)
}

# Christensen
{
  # traceplot(CCPM_coda_info)
  gelman.diag(CCPM_coda_info)
  gelman.plot(CCPM_coda_info)
}
# Mylius
{
  # traceplot(MCPM_coda_info)
  gelman.diag(MCPM_coda_info)
  gelman.plot(MCPM_coda_info)
}
    
# Combined with Christensen and informed prior
{
  gelman.diag(Ccombo_coda)
  gelman.plot(Ccombo_coda)
}

# Combined with Mylius and informed prior
{
  gelman.diag(Mcombo_coda)
  gelman.plot(Mcombo_coda)
}
}
}


# Tables with the parameters for informed priors
get_table_with_informed_prior_Mikkela <- function(){
  # i = 2
  informed_prior <- lapply(informed_prior_Mikkela, round, 3)
  
  parameter_name <- c(# Mikkela
    "logit($pb_1$)",  "$\\tau_e$", "$e_1$", "$p_w$",  
    "$\\tau_w$", "$\\mu_0$", "$\\sigma_{b}$")
  distr_name <- c("Normal", "LogNormal", "Normal", "Beta", "Gamma",
                  "Normal",  "Gamma")
  
  thestring <- list()
  for(i in 1:length(informed_prior)){
  
    thestring[i] <- paste(parameter_name[i], " & " , distr_name[i], " & ",
                     colnames(informed_prior[[i]])[1], " = ", informed_prior[[i]][1], 
                     " & ", colnames(informed_prior[[i]])[2], " = ", informed_prior[[i]][1,2],
                     "\\\\")
  
    # writeLines(thestring[i])
    whole_string <- do.call(paste, append(thestring, c(sep = "\n")))
    writeLines(whole_string, "info_table_Mikkela.txt")
  }
}
get_table_with_informed_prior_Christensen <- function(){
  # i = 2
  informed_prior <- lapply(informed_prior_ChristensenCPM, round, 3)
 
  parameter_name <- c("$t.t_{ce}$", "$s.t_{ce}$",
                       "$t.t_{es}$", "$s.t_{es}$")
  
  distr_name <- c("Beta", "Gamma", "Beta", "Gamma")
  shape_rate <- c("shape1", "shape2",
                  "shape", "rate")
  
  thestring <- list()
  # for(i in 1:length(parameter_name)){
  #   for(j in 1:2){
  #     thestring[i] <- paste(parameter_name[i], " & " , distr_name[i], " & ",
  #                           shape_rate[i], " = ", 
  #                           informed_prior[[i]][j,1], " & ", 
  #                           shape_rate[i], " = ", 
  #                           informed_prior[[i]][j,2], "\\\\")
      
      thestring[1] <- paste(parameter_name[1], " & " , distr_name[1], " & ",
                                 shape_rate[1], " = ", 
                                 informed_prior[[1]][1,1], " & ", 
                                 shape_rate[2], " = ", 
                                 informed_prior[[1]][1,2], "\\\\")
      
      thestring[2] <- paste(parameter_name[2], " & " , distr_name[2], " & ",
                            shape_rate[3], " = ", 
                            informed_prior[[1]][2,1], " & ", 
                            shape_rate[4], " = ", 
                            informed_prior[[1]][2,2], "\\\\")
      
      thestring[3] <- paste(parameter_name[3], " & " , distr_name[3], " & ",
                            shape_rate[1], " = ", 
                            informed_prior[[2]][1,1], " & ", 
                            shape_rate[2], " = ", 
                            informed_prior[[2]][1,2], "\\\\")
      
      thestring[4] <- paste(parameter_name[4], " & " , distr_name[4], " & ",
                            shape_rate[3], " = ", 
                            informed_prior[[2]][2,1], " & ", 
                            shape_rate[4], " = ", 
                            informed_prior[[2]][2,2], "\\\\")
      
      
  
      whole_string <- do.call(paste, append(thestring, c(sep = "\n")))
      writeLines(whole_string, "info_table_Christensen.txt")
    # }
  # }
}
get_table_with_informed_prior_Mylius <- function(){
  # i = 2
  informed_prior <- lapply(informed_prior_MyliusCPM, round, 3)
  
  parameter_name <- c("$t.t_{ch}$", "$s.t_{ch}$",
                      "$t.t_{cb}$", "$s.t_{cb}$",
                      "$mu.t_{bs}$", "$precision.t_{bs}$")
  
  distr_name <- c("Beta", "Gamma", "Beta", "Gamma", "Norm", "Gamma")
  shape_rate <- c("shape1", "shape2",
                  "shape", "rate",
                  "mean", "sd")
  
  thestring <- list()
  # for(i in 1:length(parameter_name)){
  #   for(j in 1:2){
  #     thestring[i] <- paste(parameter_name[i], " & " , distr_name[i], " & ",
  #                           shape_rate[i], " = ", 
  #                           informed_prior[[i]][j,1], " & ", 
  #                           shape_rate[i], " = ", 
  #                           informed_prior[[i]][j,2], "\\\\")
  
  thestring[1] <- paste(parameter_name[1], " & " , distr_name[1], " & ",
                        shape_rate[1], " = ", 
                        informed_prior[[1]][1,1], " & ", 
                        shape_rate[2], " = ", 
                        informed_prior[[1]][1,2], "\\\\")
  
  thestring[2] <- paste(parameter_name[2], " & " , distr_name[2], " & ",
                        shape_rate[3], " = ", 
                        informed_prior[[1]][2,1], " & ", 
                        shape_rate[4], " = ", 
                        informed_prior[[1]][2,2], "\\\\")
  
  thestring[3] <- paste(parameter_name[3], " & " , distr_name[3], " & ",
                        shape_rate[1], " = ", 
                        informed_prior[[2]][1,1], " & ", 
                        shape_rate[2], " = ", 
                        informed_prior[[2]][1,2], "\\\\")
  
  thestring[4] <- paste(parameter_name[4], " & " , distr_name[4], " & ",
                        shape_rate[3], " = ", 
                        informed_prior[[2]][2,1], " & ", 
                        shape_rate[4], " = ", 
                        informed_prior[[2]][2,2], "\\\\")
  
  thestring[5] <- paste(parameter_name[5], " & " , distr_name[5], " & ",
                        shape_rate[5], " = ", 
                        informed_prior[[3]][1,1], " & ", 
                        shape_rate[6], " = ", 
                        informed_prior[[3]][1,2], "\\\\")
  
  thestring[6] <- paste(parameter_name[6], " & " , distr_name[6], " & ",
                        shape_rate[3], " = ", 
                        informed_prior[[3]][2,1], " & ", 
                        shape_rate[4], " = ", 
                        informed_prior[[3]][2,2], "\\\\")
  
  whole_string <- do.call(paste, append(thestring, c(sep = "\n")))
  writeLines(whole_string, "info_table_Mylius.txt")
  # }
  # }
}
get_table_with_informed_prior_Teunis <- function(){
  # i = 2
  informed_prior <- lapply(informed_prior_Teunis, round, 3)
  
  parameter_name <- c("$\\mu_{1, strain}$", "$\\mu_{1, host}$",
                      "$\\mu_{2, strain}$", "$\\mu_{2, host}$", 
                      "$\\mu_{3, strain}$", "$\\mu_{3, host}$", 
                      "$z_1$","$z_2$", "$z_2$, outbreak only",
                      "log conc", "logit p_0")
  distr_name <- c("Normal")
  
  thestring <- list()
  k <- 0
  for(i in 1:3){
    for(j in c(1, 3)){
      k <- k+1
      thestring[k] <- paste(parameter_name[k], " & " , distr_name[1], 
                          " & ", informed_prior[[1]][i,j], 
                          " & ", informed_prior[[1]][i,j+1],
                          "\\\\")
      
    # writeLines(thestring[i])
    
  }
}

  # k = 7
  for(i in 1:3){
    k <- k+1
    thestring[k] <- paste(parameter_name[k], " & " , distr_name[1], 
                          " & ", informed_prior[[2]][i,1], 
                        " & ", informed_prior[[2]][i,2],
                        "\\\\")
  }
  
  k <- k+1
  thestring[k] <- paste(parameter_name[k], " & " , distr_name[1], 
                          " & ", informed_prior[[3]][1], 
                          " & ", informed_prior[[3]][2],
                          "\\\\")
  k <- k+1
  thestring[k] <- paste(parameter_name[k], " & " , distr_name[1], 
                        " & ", informed_prior[[4]][1], 
                        " & ", informed_prior[[4]][2],
                        "\\\\")
  
  whole_string2 <- do.call(paste, append(thestring, c(sep = "\n")))
  writeLines(whole_string2, "info_table_Teunis.txt")
  
  
}
  

# Plot of the posterior
plot_the_posterior_Mikkela <- function(post_orig, post_info, post_combi_C, post_combi_M){
  post_orig <- post_Mikkela_orig_prior
  post_info <- post_Mikkela_info_prior
  post_combi_C <- post_combi_Christensen_info_prior
  post_combi_M <- post_combi_Mylius_info_prior
  
  posterior_names <- c(# Mikkela
                        "logitpb_1", "taue","e_1", "pw", "tauw",
                        "mu_zero", "sdb")
  
  par(mfcol=c(1,1))
  # par(mfcol=c(4,2))
  # Posterior probability density
  for(i in 1:length(posterior_names)){
  # i = 1
    minhist_C <- min(post_info[[i]], post_combi_C[[i]], 
                     post_combi_M[[i]])
    maxhist_C <- max(post_info[[i]], post_combi_C[[i]], 
                     post_combi_M[[i]])
    
    plot(density(post_combi_C[[i]]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[i])
    lines(density(post_orig[[i]]),  col='black')#, lty=2)
    lines(density(post_info[[i]]),  col='blue')
    lines(density(post_combi_M[[i]]),  col='green')
    if(posterior_names[i] == "taue" || posterior_names[i] == "tauw" || posterior_names[i] == "sdb"){
      legend("topright", bg="transparent", legend=c("Original model",
                                  "Prior-informed model",
                                  "Combined model, Christensen",
                                  "Combined model, Mylius"),
             col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1,
             box.lty=0)
    } else{
      legend("topleft", bg="transparent", legend=c("Original model",
                                  "Prior-informed model",
                                  "Combined model, Christensen",
                                  "Combined model, Mylius"),
             col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1,
             box.lty=0)
    }
    plotpath<- file.path("PicForReport2",paste("Mikkelaplot_",posterior_names[i],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF 
    plot(ecdf(post_combi_C[[i]]),  col='darkred', #xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[i])
    lines(ecdf(post_orig[[i]]),  col='black')
    lines(ecdf(post_info[[i]]),  col='blue')
    lines(ecdf(post_combi_M[[i]]),  col='green')
    if(posterior_names[i] == "taue"){
      legend("bottomright", bg="transparent", legend=c("Original model",
                                                       "Prior-informed model",
                                                       "Combined model, Christensen",
                                                       "Combined model, Mylius"),
             col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1,
             box.lty=0)
    } else{
      legend("topleft", bg="transparent", legend=c("Original model",
                                                   "Prior-informed model",
                                                   "Combined model, Christensen",
                                                   "Combined model, Mylius"),
             col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1,
             box.lty=0)
    }
    plotpath<- file.path("PicForReport2",paste("Mikkelaplot_",posterior_names[i],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }

}
plot_the_posterior_Christensen <- function(post_orig, post_info, post_combi_C){
  post_orig <- post_ChristensenCPM_orig_prior
  post_info <- post_ChristensenCPM_info_prior
  post_combi_C <- post_combi_Christensen_info_prior
  
  posterior_names <- c(
    # Christensen
    "t_t_ce", "s_t_ce",
    "t_t_es", "s_t_es")
  
  par(mfcol=c(1,1))
  # par(mfcol=c(2,2))
  k = 7
  for(i in 1:length(posterior_names)){
    #i = 4
    minhist_C <- min(post_info[[2+i]], post_combi_C[[k+i]], post_orig[[2+i]])
    maxhist_C <- max(post_info[[2+i]], post_combi_C[[k+i]], post_orig[[2+i]])
    #maxhist_Cy <- max(post_info[[2+i]], post_combi_C[[k+i]], post_orig[[2+i]])
    plot(density(post_combi_C[[k+i]]),  col='darkred', xlim=c(minhist_C, maxhist_C),
          xlab = "", main = posterior_names[i])
    lines(density(post_orig[[2+i]]),  col='black') #, lty=2)
    lines(density(post_info[[2+i]]),  col='blue')
    legend("topright", bg="transparent", legend=c("Original model",
                               "Prior-informed model",
                               "Combined model, Christensen"),
           col=c("black", "blue","darkred"), lty=c(1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Christensenplot_",posterior_names[i],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_C[[k+i]]),  col='darkred',#, xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[i])
    lines(ecdf(post_orig[[2+i]]),  col='black') #), lty=2)
    lines(ecdf(post_info[[2+i]]),  col='blue')
    legend("bottomright", bg="transparent", legend=c("Original model",
                                                     "Prior-informed model",
                                                     "Combined model, Christensen"),
           col=c("black", "blue","darkred"), lty=c(1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Christensenplot_",posterior_names[i],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }
  
  
}
plot_the_posterior_Mylius <- function(post_orig, post_info, post_combi_M){
  # post_orig <- post_MyliusCPM_orig_prior
  # post_info <- post_MyliusCPM_info_prior
  # post_combi_M <- post_combi_Mylius_info_prior
  
  posterior_names <- c(
    # Mylius
    "t_t_ch", "s_t_ch",
    "t_t_cb", "s_t_cb",
    "mu_t_bs", "precision_t_bs")
  
  k = 7
  par(mfcol=c(1,1))
  # par(mfcol=c(3,2))
  for(i in 1:6){
    minhist_C <- min(post_info[[i]], post_combi_M[[k+i]])
    maxhist_C <- max(post_info[[i]], post_combi_M[[k+i]])
    plot(density(post_combi_M[[k+i]]),  col='green', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[i])
    lines(density(post_orig[[i]]),  col='black')#, lty=2)
    lines(density(post_info[[i]]),  col='blue')
    if(posterior_names[i] == "mu_t_bs"){
      legend("topleft", bg="transparent", legend=c("Original model",
                                 "Prior-informed model",
                                 "Combined model, Mylius"),
             col=c("black", "blue", "green"), lty=c(1,1,1), cex=1.2,
            box.lty=0)
    }else{
      legend("topright", bg="transparent", legend=c("Original model",
                                                   "Prior-informed model",
                                                   "Combined model, Mylius"),
             col=c("black", "blue", "green"), lty=c(1,1,1), cex=1.2,
             box.lty=0)
    }
    plotpath<- file.path("PicForReport2",paste("Myliusplot_",posterior_names[i],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_M[[k+i]]),  col='green', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[i])
    lines(ecdf(post_orig[[i]]),  col='black')#, lty=2)
    lines(ecdf(post_info[[i]]),  col='blue')
    legend("bottomright", bg="transparent", legend=c("Original model",
                                                     "Prior-informed model",
                                                     "Combined model, Mylius"),
           col=c("black", "blue", "green"), lty=c(1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Myliusplot_",posterior_names[i],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }

  
}
plot_the_posterior_Teunis <- function(post_orig, post_info, post_combi_C, post_combi_M){
  post_orig <- post_Teunis_orig_prior
  post_info <- post_Teunis_info_prior
  post_combi_C <- post_combi_Christensen_info_prior
  post_combi_M <- post_combi_Mylius_info_prior

  posterior_names <- c("mu_w11", "mu_w21", "mu_w31", "mu_w12", "mu_w22", "mu_w32",
                       "z1_92","z1_93", 
                       "z2_92","z2_93",
                       "logconc_5", "ltp0_5")
  
  # mu.w
  par(mfcol=c(1,1))
  # par(mfcol=c(3,2))
  for(j in 1:6){
    #j = 6
    minhist_C <- min(post_info$mu.w[,j], post_combi_C$mu.w[,j], 
                     post_combi_M$mu.w[,j])
    maxhist_C <- max(post_info$mu.w[,j], post_combi_C$mu.w[,j], 
                     post_combi_M$mu.w[,j])
    plot(density(post_combi_C$mu.w[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[j])
    lines(density(post_orig$mu.w[,j]),  col='black')#, lty=2)
    lines(density(post_info$mu.w[,j]),  col='blue')
    lines(density(post_combi_M$mu.w[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                               "Prior-informed model",
                               "Combined model, Christensen",
                               "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[j],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_C$mu.w[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[j])
    lines(ecdf(post_orig$mu.w[,j]),  col='black')#, lty=2)
    lines(ecdf(post_info$mu.w[,j]),  col='blue')
    lines(ecdf(post_combi_M$mu.w[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                                                     "Prior-informed model",
                                                     "Combined model, Christensen",
                                                     "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[j],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }
  
  # z1
  name <- 6
  for(j in c(18, 27)){
    # j = 18
    name <- name+1
    minhist_C <- min(post_info$z1[,j], post_combi_C$z1[,j], 
                     post_combi_M$z1[,j])
    maxhist_C <- max(post_info$z1[,j], post_combi_C$z1[,j], 
                     post_combi_M$z1[,j])
    plot(density(post_combi_C$z1[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[name])
    lines(density(post_orig$z1[,j]),  col='black') #, lty=2)
    lines(density(post_info$z1[,j]),  col='blue')
    lines(density(post_combi_M$z1[,j]),  col='green')
    legend("topleft",  bg="transparent", legend=c("Original model",
                               "Prior-informed model",
                               "Combined model, Christensen",
                               "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[name],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_C$z1[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[name])
    lines(ecdf(post_orig$z1[,j]),  col='skyblue')#, lty=2)
    lines(ecdf(post_info$z1[,j]),  col='skyblue')
    lines(ecdf(post_combi_M$z1[,j]),  col='darkgreen')
    legend("topleft",  bg="transparent", legend=c("Original model",
                                "Prior-informed model",
                                "Combined model, Christensen",
                                "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[name],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }
  
  # z2
  for(j in c(18, 27)){
    # j = 18
    name <- name+1
    minhist_C <- min(post_info$z2[,j], post_combi_C$z2[,j], 
                     post_combi_M$z2[,j])
    maxhist_C <- max(post_info$z2[,j], post_combi_C$z2[,j], 
                     post_combi_M$z2[,j])
    plot(density(post_combi_C$z2[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[name])
    lines(density(post_orig$z2[,j]),  col='black')#, lty=2)
    lines(density(post_info$z2[,j]),  col='blue')
    lines(density(post_combi_M$z2[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                                                  "Prior-informed model",
                                                  "Combined model, Christensen",
                                                  "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[name],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_C$z2[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[name])
    lines(ecdf(post_orig$z2[,j]),  col='black')#, lty=2)
    lines(ecdf(post_info$z2[,j]),  col='blue')
    lines(ecdf(post_combi_M$z2[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                                                  "Prior-informed model",
                                                  "Combined model, Christensen",
                                                  "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[name],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }
  
  # logconc
  {
    j = 5
    minhist_C <- min(post_info$logconc[,j], post_combi_C$logconc[,j], 
                     post_combi_M$logconc[,j])
    maxhist_C <- max(post_info$logconc[,j], post_combi_C$logconc[,j], 
                     post_combi_M$logconc[,j])
    plot(density(post_combi_C$logconc[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[11])
    lines(density(post_orig$logconc[,j]),  col='black')#, lty=2)
    lines(density(post_info$logconc[,j]),  col='blue')
    lines(density(post_combi_M$logconc[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                               "Prior-informed model",
                               "Combined model, Christensen",
                               "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[11],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_C$logconc[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[11])
    lines(ecdf(post_orig$logconc[,j]),  col='black')#, lty=2)
    lines(ecdf(post_info$logconc[,j]),  col='blue')
    lines(ecdf(post_combi_M$logconc[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                                                  "Prior-informed model",
                                                  "Combined model, Christensen",
                                                  "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[11],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }
  
  # ltp0
  {
    j = 5  
    minhist_C <- min(post_info$ltp0[,j], post_combi_C$ltp0[,j], 
                     post_combi_M$ltp0[,j])
    maxhist_C <- max(post_info$ltp0[,j], post_combi_C$ltp0[,j], 
                     post_combi_M$ltp0[,j])
    plot(density(post_combi_C$ltp0[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[12])
    lines(density(post_orig$ltp0[,j]),  col='black')#, lty=2)
    lines(density(post_info$ltp0[,j]),  col='blue')
    lines(density(post_combi_M$ltp0[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                               "Prior-informed model",
                               "Combined model, Christensen",
                               "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[12],".png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
    
    # CDF
    plot(ecdf(post_combi_C$ltp0[,j]),  col='darkred', xlim=c(minhist_C, maxhist_C),
         xlab = "", main = posterior_names[12])
    lines(ecdf(post_orig$ltp0[,j]),  col='black')#, lty=2)
    lines(ecdf(post_info$ltp0[,j]),  col='blue')
    lines(ecdf(post_combi_M$ltp0[,j]),  col='green')
    legend("topleft", bg="transparent", legend=c("Original model",
                                                  "Prior-informed model",
                                                  "Combined model, Christensen",
                                                  "Combined model, Mylius"),
           col=c("black", "blue","darkred", "green"), lty=c(1,1,1,1), cex=1.2,
           box.lty=0)
    plotpath<- file.path("PicForReport2",paste("Teunisplot_",posterior_names[12],"_CDF.png",sep=""))
    dev.print(png, plotpath, width = 600, height = 500)
  }

}


# Plot of pf, mean dose and prill
plot_pf_meandose_prill <- function(post, maxy, marginalline, CPM){
  
 # post <- post_combi_Mylius_info_prior
 # maxy=0.1
 # marginalline=3
 
  pf <- post$pf
  meandose <- post$dose_pred_mean
  prill_human <- post$prill_pred_human_tot 
  prill_outbreak <- post$prill_pred_outbreak_tot
  prinf <- post$prinf_pred_human_tot
  
  par(mar=c(5,4,2,4)+0.1)
  plot(c(1:12)-0.2, apply(pf, 2, 'quantile', 0.5) ,ylim=c(0,0.5),xlim=c(1, 12),
       ylab="Predicted probability", xlab="Month", axes=FALSE, main=CPM)
  segments(x0=c(1:12)-0.2, y0=apply(pf,2,'quantile',0.05),
           y1=apply(pf,2,'quantile',0.95))
  
  lines(c(1:12)+0.2, apply(prinf, 2, 'quantile', 0.5), type="p", col="blue" )
  segments(x0=c(1:12)+0.2, y0=apply(prinf,2,'quantile',0.05),
           y1=apply(prinf,2,'quantile',0.95),col="blue")
  
  lines(c(1:12)+0.4,apply(prill_human,2,'quantile',0.5), type="p", col="darkgreen")
  segments(x0=c(1:12)+0.4,y0=apply(prill_human,2,'quantile',0.05),
         y1=apply(prill_human,2,'quantile',0.95), col="darkgreen")
  
  lines(c(1:12)+0.3, apply(prill_outbreak, 2, 'quantile', 0.5), type="p", col="purple")
  segments(x0=c(1:12)+0.3, y0=apply(prill_outbreak,2,'quantile',0.05),
           y1=apply(prill_outbreak,2,'quantile',0.95), col="purple")
  
  axis(2, ylim=c(0,0.5),col="black",las=1)
  axis(1, xlim=c(1, 12), at=1:12, label=month.abb, xlab="Month")

  par(new=TRUE)
  plot(c(1:12)-0.1, apply(meandose, 2, 'quantile',0.5), type="p", col="red", 
       ylim=c(0,maxy),  axes = FALSE, ylab="Predicted probability", xlab="Month")
  segments(x0=c(1:12)-0.1,y0=apply(meandose,2,'quantile',0.05),
         y1=apply(meandose,2,'quantile',0.95), col="red")
  axis(4, ylim=c(0,maxy), col="red",col.axis="red",las=1, line=0.2)
  mtext("Predicted dose, cfu", side=4, line=marginalline, col="red")

  legend("topleft", legend=c("Prevalence, pf",
                            "Predicted mean dose",
                            "Predicted prob. of infect.", 
                             "Predicted prob. of ill., human",
                             "Predicted prob. of ill., outbreak"),
         col=c("black", "red","blue", "darkgreen", "purple"), lty=c(1,1,1), cex=0.8,
         box.lty=0)
  #plotpath<- file.path("PicForReport2",paste("pf_dose_prill_prinf_",CPM,".png",sep=""))
  #dev.print(png, plotpath, width = 700, height = 448)
  
}

  
# Functions to get the tables with informed priors


### Case study: Plot predictions of number of Campylobacteriosis cases vs the Finnish data
plot_pred_case_study <- function(z1, z2, finnishdata){#, finnishdata, CPM){
  
  # maxy=0.1
  # marginalline=3
  findata <- finnishdata$val
    
  Ymaxium <- 10000
  # Christensen
  par(mar=c(5,4,2,4)+0.1)
  plot(c(1:12), apply(z1, 2, 'quantile', 0.5), ylim=c(0,Ymaxium), xlim=c(1, 12),
       ylab="Predicted number of campylobacteriosis cases", xlab="Month",  axes=FALSE, 
       main="Predicted number of campylobacteriosis cases vs the Finnish data")
  segments(x0=c(1:12), y0=apply(z1,2,'quantile',0.05),
           y1=apply(z1,2,'quantile',0.95))
  
  axis(2, ylim=c(0,Ymaxium),col="black",las=1)
  axis(1, xlim=c(1, 12), at=1:12, label=month.abb, xlab="Month")
  
  par(new=TRUE)
  plot(c(1:12), c(rep(NA, 10), findata[1:2]), type="p", pch=15, col="blue",  axes = FALSE, 
       ylim=c(0, 1000), ylab="", xlab="Month")
  lines(c(1:12), c(findata[3:13], NA), type="p", pch=17, col="blue", axes = FALSE,
        ylim=c(0, 1000))
  lines(c(1:12), c(rep(NA, 5), findata[14:16], rep(NA, 4)), type="p", pch=19, col="blue", axes = FALSE, 
        ylim=c(0, 1000))
  axis(4, ylim=c(0,1000), col="blue",col.axis="blue",las=1, line=0.2)
  mtext("Number of campylobacteriosis cases, Finnish data", side=4, line=3, col="blue")
  
  legend("topleft", legend=c("Model with Christensen CPM",
                             "Finnish campylobacteriosis data, 2012", 
                             "Finnish campylobacteriosis data, 2013",
                             "Finnish campylobacteriosis data, 2014"),
         col=c("black", "blue", "blue", "blue"), pch=c(1, 15, 17, 19),cex=0.8,
         box.lty=0)
  
  plotpath<- file.path("PicForReport2",paste("Prediction_study_case_Christensen.png",sep=""))
  dev.print(png, plotpath, width = 700, height = 448)
  
  # Mylius
  par(mar=c(5,4,2,4)+0.1)
  plot(c(1:12), apply(z2, 2, 'quantile', 0.5),ylim=c(0,10),xlim=c(1, 12),
       ylab="Predicted number of campylobacteriosis cases", xlab="Month", axes=FALSE, 
       main="Predicted number of campylobacteriosis cases vs the Finnish data")
  segments(x0=c(1:12), y0=apply(z2,2,'quantile',0.05),
           y1=apply(z2,2,'quantile',0.95))
  
  axis(2, ylim=c(0,10000),col="black",las=1)
  axis(1, xlim=c(1, 12), at=1:12, label=month.abb, xlab="Month")
  
  par(new=TRUE)
  plot(c(1:12), c(rep(NA, 10), findata[1:2]), type="p", pch=15, col="blue",
       axes = FALSE, ylim=c(0, 1000), ylab="", xlab="Month")
  lines(c(1:12), c(findata[3:13], NA), type="p", pch=17, col="blue",
        axes = FALSE, ylim=c(0, 1000))
  lines(c(1:12), c(rep(NA, 5), findata[14:16], rep(NA, 4)), type="p", pch=19, col="blue",
        axes = FALSE, ylim=c(0, 1000))
  axis(4, ylim=c(0,1000), col="blue",col.axis="blue",las=1, line=0.2)
  mtext("Number of campylobacteriosis cases, Finnish data", side=4, line=3, col="blue")
  
  legend("topleft", legend=c("Model with Mylius CPM",
                             "Finnish campylobacteriosis data, 2012", 
                             "Finnish campylobacteriosis data, 2013",
                             "Finnish campylobacteriosis data, 2014"),
         col=c("black", "blue", "blue", "blue"), pch=c(1, 15, 17, 19),cex=0.8,
         box.lty=0)
  
  
  plotpath<- file.path("PicForReport2",paste("Prediction_study_case_Mylius.png",sep=""))
  dev.print(png, plotpath, width = 700, height = 448)
  
}

if(FALSE){

##########################################################
# Function that produce the SA results


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
load("sammatC.Rdata")  

xxC <- sammatC[,sort(unique(unlist(lapply(1:length(Christensen_factor),function(i){
  grep(Christensen_factor[i],colnames(sammatC))}))))]
# colnames(xxC)

xxM <- sammatM[,sort(unique(unlist(lapply(1:length(Mylius_factor),function(i){
  grep(Mylius_factor[i],colnames(sammatM))}))))]


# OBS: Rewrite this as a function
# generate_the_SA_result_Christensen <- function(xxC,){
  
  plot_names <- c('dose',
                  'prinf_outbreak',
                  'prinf_human',
                  'prill_outbreak', 
                  'prill_human', 
                  'pf')
  
  
  # par(mfrow = c(3,2))
  for(i in 1:6) {# length(allY)){

    y <- sammatC[,grep(SA_output[i],colnames(sammatC))][,7] # 7 for the july month
    k <- ncol(xxC)
    
    ## Först anpassar man en responsyta (polynom regression) med y som respons och X som prediktorer
    mod[[i]]<-bass(xxC,y)
    # plot(mod) ## ritar upp residualer och kollar att ser ok ut
    
    ## Beräknar Sobols sensitivity indices
    sens[[i]]<-sobol(mod[[i]])
    # plot(sens,cex.axis=.5) 
    
    # colnames(mod$xx.des)
    
    ## skapar data frame med alla känslighetsvärden 
    ## plocka ut värden och namn på parametrar
    Tind <- as.numeric(colnames(sens[[i]]$T))
    dat <- do.call('rbind',lapply(1:ncol(sens[[i]]$T),function(j){
      data.frame(first_order_sensitivity=sens[[i]]$S[,j],total_sensitivity = sens[[i]]$T[,j],
                 param = rep(colnames(mod$xx.des)[Tind[j]],nrow(sens[[i]]$S)))
    }))
    
    #dat$param <- as.factor(dat$param)
    
    ## even valja ut de param som har strorst inflytande
    first_mean <- aggregate(dat$first_order_sensitivity,by=list(dat$param),mean)
    total_mean <- aggregate(dat$total_sensitivity,by=list(dat$param),mean)
    bet1 <- order(first_mean$x,decreasing=TRUE)
    bet2 <- order(total_mean$x,decreasing=TRUE)
    # antalduvillha
    antalduvillha=10
    
    #dat_urval <- dat
    dat_urval1 <- dat[dat$param %in%  first_mean$Group.1[bet1[1:antalduvillha]],]
    p1[[i]] <- ggplot(dat_urval1,aes(x=param,y=first_order_sensitivity,fill = param)) +
      geom_boxplot(width=0.15) +
      geom_jitter(shape=16, position=position_dodge(0.8), cex = 0.5) +
      labs(title=paste('Sobol\'s First order sensitivity,', plot_names[i]),
           x="parameter", y = "proportion variance")

    
    dat_urval2 <- dat[dat$param %in%  total_mean$Group.1[bet2[1:antalduvillha]],]
    p2[[i]] <- ggplot(dat_urval2,aes(x=param,y=total_sensitivity,fill = param)) +
      geom_boxplot(width=0.15) +
      geom_jitter(shape=16, position=position_dodge(0.8), cex = 0.5) +
      labs(title=paste('Sobol\'s total sensitivity,', plot_names[i]),x="parameter", y = "proportion variance")
    
   

  }
  
  # list_p <- list(first=p1, total=p2)
  # return(list_p)
    
# }

# save(p1, p2,  file="SAresult.RData") #obs, Christensen only
load("SAresult.Rdata")

pdf("test11.pdf") #, width = 800, height = 700) 
for(i in 1:6) {# length(allY)){
  print(p1[[i]])
  # ggsave("plot.pdf", arrangeGrob(grobs = l))
  print(p2[[i]]) # pdf
  # ggsave("plot.pdf", arrangeGrob(grobs = l))
}
dev.off() 

