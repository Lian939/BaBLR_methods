# BaBLR_methods
This is the R code for "BAyesian Bent-Line Regression model for longitudinal data with an application to the study of cognitive performance trajectories in Wisconsin Registry for Alzheimer’s Prevention".

In the R script file entitled “R_example_fitting_BAyesian_Bent-Line_Regression _model_in_Stan.R”, we have provided R code which allows the user to create a single simulated dataset (similar to the WRAP dataset used in the main analysis of our paper) and then fit the BAyesian Bent-Line Regression model to this simulated dataset using Stan. Comments are also provided in the R script file which should explain the main steps required to fit the model, as well as looking at some diagnostics (e.g. trace plots, density plot, and potential scale reduction statistic values summary plot) and inference (e.g. model parameter estimates).

In the plain text file entitled “BAyesian_Bent-Line_Regression_model.stan”, we provide the Stan model code for fitting the BAyesian Bent-Line Regression model described in the main manuscript. 

In the RData file entitled “stanfitsim_n100_10000iteration_example” we provide a model fitting result on one random simulated dataset (n = 100).
