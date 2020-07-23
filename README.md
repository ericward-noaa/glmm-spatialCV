# glmm-spatialCV

This repo contains some simulations to ask whether spatiotemporal GLMM methods are affected by the need to do blocked spatial cross-validation. Several recent papers (Valavi et al. 2018)[https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13107] have developed new tools for blocked cross validation for SDMs, because not accounting for spatial autocorrelation may lead to overfitting the training data. 

There are 3 scenarios currently explored here in the simulations:

1. Example 1 contains a logistic model, without time dependence (single spatial field in a single year)

2. Example 2 is like Example 1, but with a normal response

3. Example 3 extends Example 2 to include a spatiotemporal component
