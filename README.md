# lmr  #

There are functions for robust linear modelling in various packages.
This small packages is largely wrappers to provide default settings
(Pena-Yohai initial values, MM-estimation with bisquare weight
functions and 85% Gaussian efficiency) that I want,
robust final prediction error and anova functions. It also has a
wrapper for Salibian-Barrera's fast bootstrap (github.com/msalibian/FRB).

The main purpose of this package is to get the defaults that I want,
have everything in one place, and to be able to compare robustbase::lmrob
with MASS::rlm for testing purposes.

