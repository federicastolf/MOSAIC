# Supporting information for "Bayesian inference on beta diversity via feature allocation models"

Code and data for the analysis contained in Stolf, Rigon & Dunson (2026).

- #### functs
  - [mosaic_NegBin.stan](functs/mosaic_NegBin.stan) and [gibbs_mosaic.R](functs/gibbs_mosaic.R) contain the R and stan code to run the MCMC algorithm for MOSAIC.
  - [Cwrapper.cpp](functs/Cwrapper.cpp) contains C++ utility functions used across the project.
  - [helper.R](functs/helper.R) contains miscellaneous R utility functions used across the project.
  - [rarefaction.R](rfunctsarefaction.R) contains R functions to compute empirical and model-based site-specific rarefaction curves and rarefaction surfaces for shared species.

- #### fungi_analysis
  - [fungiGSSP.Rdata](fungi_analysis/fungiGSSP.Rdata) are the fungi data from the Global Spore Sampling Project analyzed in the paper.
  - [fitFungi.R](fungi_analysis/fitFungi.R) contains the code to replicate the fungal biodiversity analysis presented in Section 5 of the paper.
  - [fungi_oos_sm.R](fungi_analysis/fungi_oos_sm.R) contains the code to replicate the out-of-sample experiments for the fungi data presented in the supplementary material.
  
- #### simulations
  - [sim_Kq.R](simulations/sim_Kq.R) contains the code to replicate the experiments of the simulation study on the predictions of the site-specific rarefaction curves.
  - [sim_coverage.R](simulations/sim_coverage.R) and [sim_sm.R](simulations/sim_sm.R) contains the code to replicate the additional simulations for coverage and varying number of sites in the supplementary material.
  - [sim_shared.R](simulations/sim_shared.R) contains the code to replicate the experiments of the simulation study on the rarefaction surfaces for shared species.

The script [example.R](example.R) shows an example of the use of the available routines.
