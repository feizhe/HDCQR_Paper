# HDCQR_Paper
R Code for Paper "Inference for High Dimensional Censored Quantile Regression"

This repo contains the R code to implement simulation example 3 in the paper "Inference for High Dimensional Censored Quantile Regression."

To implement the simulation example, run the main code "jobEg3.R", the rest are helper codes. I would recommend to start with smaller values of n, p, and B and test the code, for example, n=200, p=300, B=20.

To run Fused-HDCQR for your own data, use "JASAfunctions.R" and "JASAparameters.R", where you only need to run the latter and specify the data and model parameters.
