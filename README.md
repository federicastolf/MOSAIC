# Supporting information for "Bayesian inference on beta diversity via feature allocation models with imperfect detection"

Code and data for the analysis contained in Stolf, Rigon & Dunson (2026).

- #### functs
  - [mosaic_NegBin.stan](mosaic_NegBin.stan) and [rarefaction.R](rarefaction.R) contain the R and stan code to run the MCMC algorithm for MOSAIC.
  - [Cwrapper.cpp](Cwrapper.cpp) contains C++ utility functions used across the project.
  - [helper.R](helper.R) contains miscellaneous R utility functions used across the project.
  - [rarefaction.R](rarefaction.R) contains R functions to compute empirical and model-based site-specific rarefaction curves and rarefaction surfaces for shared species.

- #### fungi_analysis
  - [fungiGSSP.Rdata](fungiGSSP.Rdata) are the fungi data from the Global Spore Sampling Project analyzed in the paper.
  - [fitFungi.R](fitFungi.R) contains the code to replicate the fungal biodiversity analysis presented in Section 5 of the paper.
  - [fungi_oos_sm.R](fungi_oos_sm.R) contains the code to replicate the out-of-sample experiments for the fungi data presented in the supplementary material.
  
- #### simulations
  - [sim_Kq.R](sim_Kq.R) contains the code to replicate the experiments of the simulation study on the predictions of the site-specific rarefaction curves.
  - [sim_coverage.R](sim_coverage.R) and [sim_sm.R](sim_sm.R) contains the code to replicate the additional simulations for coverage and varying number of sites in the supplementary material.
  - [sim_shared.R](sim_shared.R) contains the code to replicate the experiments of the simulation study on the rarefaction surfaces for shared species.
