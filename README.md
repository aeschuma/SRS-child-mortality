# SRS-child-mortality
Code to reproduce run simulations and estimate child mortality from MCHSS data. 
- each "task" (i.e. simulations, analyzing the MCHSS data, making plots for results) has its own folder, and each of these has its own R project
- the R project denotes the home directory that all of the file paths for that task's specific code is based.
- all code is set up to save results into a Dropbox folder named 'SRS-child-mortality-output' via relative filepaths in relation to where the code is located
    - this is based on the code repo being cloned onto your Desktop, and that the root for Desktop is on the same level as the Dropbox---if this isn't the case for your operating system, you will have to change the relative filepath for the 'savedir' object in all code
    - subdirectories for SRS-child-mortality-output
        - SRS-child-mortality-output/results
        - SRS-child-mortality-output/graphs
            - /model_results
            - /he_compare
            - /model_testing
            - /fe_re_glms
            - /rw_analysis
    - these directories and subdirectories will get created in the code if you haven't created them already

NOTE: MCHSS data, as well as the comparison data from He et al (2017) are not publicly available. We have provided simulated test data sets to use with the code.

## Folder structure and file descriptions:

data # any test data for use with the code in this repository will be available in this folder

- mchss_test_data.RDS # simulated data set for use with MCHSS analysis code, and some simulation code
- he_2017_test_data # simulated data set of results similar to He et al (2017) for use with comparison plotting code

MCHSS_anlysis # code for the analysis of MCHSS data

- analyze_results  # code to produce graphs to analyze the results of the MCHSS modeling
  - compare_he_2017.R # produce graphs to compare our final estiamted cause fractions to those of He et al (2017) and the empirical estimates
  - csmf_plots.R # plot cause fractions over time
  - heatmaps.R # make heatmaps of the estimated change in mortality from 1996 to 2015
  - residual_plots.R # plot residuals of our final model across various combinations of region, age, time, and cause

- model_development # code for the analyses to choose fixed and random effect specification
  - fixed_effects_glms.R # fit GLMs to MCHSS data with all combinations of fixed effect interactions for region, age, and cause, then compare the AIC
  - rw_sd_analysis.R # fit separate random walk models for each strata for which we fit a random walk in the MCHSS data, estimate the standard deviation parameters, and plot them 
  - time_interaction_glms.R # fit many GLMs to the MCHSS data with different interactions between year and other dimensions, then plot the residuals

- model_fitting # code to fit the models to the MCHSS data and produce plots of results
  - model_fit_inla.R # fit our final model to the MCHSS data in INLA

simulations # code for various simulations

- compile_results_simulation_comparison_correlation.R # code to compile simulations from 'simulation_comparison_correlation.R' and plot results
- compile_results_simulation_comparison_no_correlation.R # code to compile simulations from 'simulation_comparison_no_correlation.R' and plot results
- simulation_comparison_correlation.R # code to simulate data with fixed effects and correlated random effects, and then fit a multistage model and a unified model
- simulation_comparison_no_correlation.R # code to simulate data with fixed effects and non-correlated random effects, and then fit a multistage model and a unified model
- simulation_correlation_exploration.R # code to simulate correlated data, fit a correctly specified model in STAN, and plot results. This is for the simulation in the supplementary materials that explores the amount of data needed to accurately estimate correlation parameters
- stan_models # .stan model files for fitting models in STAN
    - model_causeRE_timeRW2constr.stan # model with fixed effects, cause-correlated random effects, and random walks on time
    - model_causeRE_timeRW2constr_youngoldsplit.stan # model with fixed effects, cause-correlated random effects, and random walks on time, and allowing the data to have different causes for young age groups and old age groups
    - model_causeRE_only_all_age.stan # model with fixed effects and cause-correlated random effects
    - model_nocorrRE_timeRW2constr.stan # model with fixed effects, non-correlated random effects, and random walks on time
    - model_timeRW2constr.stan # model with fixed effects and random walk on time
