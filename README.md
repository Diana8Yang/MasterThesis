# Master's Thesis

In my [master's thesis](https://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=9100273&fileOId=9100277) at Lund University I examined the possibility of making a Quantitative Microbiological Risk Assessment (QMRA) of Campylobacter contamination in a Bayesian framework. The purpose of the work is to quantity uncertainty in parameters and model predictions using Bayesian calibration. 

This is the code repository of my master’s thesis. 

# JAGS models
To quickly understand the structure of my code, especially the different JAGS models, a brief overview of the models used in this thesis is presented here. The linked model consists by three “submodels” covering three different stages: 

1.	**Contamination of retail meat:** Prevalence and concentration model (prevalence model in short) which estimates the prevalence and concentration of Campylobacter in retail chicken meat. The model by [**Mikkela et al.**](https://onlinelibrary.wiley.com/doi/10.1111/risa.12572) is used as the prevalence model, The objective of the particular prevalence model was to estimate prevalence and concentration in retail meat using Bayesian methods while the data is censored and clustered. 

    The model in OpenBUGS was published in the article. In my thesis, I converted the model into JAGS. 

2.	**From retail meat to human exposure:** consumer phase model, CPM, which describes the contamination of food during preparation and handling of raw chicken meat in kitchen. Here we adopted two out of the CPMs used in the study of [Nauta et al.](https://onlinelibrary.wiley.com/doi/10.1111/j.1539-6924.2010.01481.x): **Christensen CPM** and **Mylius CPM**. 

    The prior distributions and hyperparameters of the CPMs was altered due to practical reasons. The data used to calibrate some of the parameters for the transfer rates in the two CPMs was collected by [Luber et al.](https://journals.asm.org/doi/10.1128/aem.72.1.66-70.2006?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed) The two models in JAGS were written by me. 

3.	**From human exposure to health effect:** dose-response model which estimates the probability of becoming infected or ill after consuming the contaminated meat, i.e. a given amount of Campylobacter. The dose-response model used by [**Teunis et al.**](https://www.sciencedirect.com/science/article/pii/S1755436517301366?via%3Dihub) was chosen in this thesis.

    The model in JAGS was published with the article and that was slightly adjusted in the thesis. 

After linking the submodels presented above, the result is two linked models with Christensen CPM respective Mylius CPM (**PrevCCPMILLChicken** respective **PrevMyliusCPMILLChicken**). 


