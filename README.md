# R code to estimate Regime switching Student-t copula models 
Code ralated to the paper titled:  **"Maximum likelihood estimation of multivariate regime switching Student-t copula model"**. Authors: Federico P. Cortese, Francesco Bartolucci and Fulvia Pennoni.

The repository contains the following files:

- *RSest_C_RScop.R* for simulating and estimating a RSStC model and for performing local and global decoding

- *RSest_C.cpp* with useful functions implemented in C++.
This file must be saved on the working directory.

- *sim_github.R* to perform the simulations proposed in the paper

- *boot_github.R* to perform non-parametric bootstrap for the estimation of standard errors of the proposed RSStC model

- *simdata5.txt*:  simulated 5-dimensional dataset of length 500 with a  similar structure of that employed for the empirical illustration in the paper related to the log-retuns of five cryptocurrencies.
