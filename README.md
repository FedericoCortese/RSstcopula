# RSstcopula
This repository containts the R functions for estimation and simulation of a regime switching Student-t copula model adopting the EM algorithm.
The model is described in the paper "Maximum likelihood estimation of multivariate regime switching Student-t copula model" by Federico P. Cortese, Francesco Bartolucci and Fulvia Pennoni.
The file contained in the repository are:
- RSest_C_RScop.R for simulating and estimating a RSStC model
- RSest_C.cpp which contains the functions implemented in C++. This file must be saved on the working directory.
- sim_github.R to performe the simulations studies explained in the paper
- boot_github.R to performe bootstrap for the estimation of standard errors
- decoding_github.R contains the R functions for the global decoding performed through the Viterbi algorithm

A simulated 5-dimensional dataset of log-returns, named simdata5.txt, is also available. It consists of 500 time observations.
