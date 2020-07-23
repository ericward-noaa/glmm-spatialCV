library(glmmfields)
library(blockCV)
library(tidyverse)
library(sdmTMB)
library(ggplot2)

# Parameters for simulation
n_sims = 50 # iterations for each combination of parameters
k_folds = 10
degree = 111325 # this is from rasterNet function in blockCV. it internally scales everything,
# so if we want ranges of 0.025, 0.05, 0.075, they need to be relative to this

# facets here will by cv (fixed range, kappa)
# variable range / kappa with cv = "random"
df = expand.grid(est_model = c("null","covariate"), oper_model = c("null","covariate"),
  "sim"=1:n_sims, "kappa" = c(2),
  "range" = c(0.025,0.05,0.075)*degree, cv = c("random","systematic","none"))
df = rbind(df,
  expand.grid(est_model = c("null","covariate"), oper_model = c("null","covariate"),
  "sim"=1:n_sims, "kappa" = c(2),
  "range" = c(0.025,0.075)*degree, cv = c("random")),
  expand.grid(est_model = c("null","covariate"), oper_model = c("null","covariate"),
    "sim"=1:n_sims, "kappa" = c(1,4),
    "range" = c(0.05)*degree, cv = c("random")))

saveRDS(df, "sims/example_1/grid.rds")

for(i in 1:nrow(df)) {

set.seed(df$sim[i]) # this ensures all the different facets are using same data

# create spatially structured environmental variable X without much observation error
x = stats::runif(300, 0, 10)
y = stats::runif(300, 0, 10)
# small values of kappa create a very smooth field. sigma_O adjusts magnitude of variance.
# fixed at 1 so estimated have roughly mean 0 and sd = 1
X = sdmTMB::sim(x = x, y = y, sigma_O =1 , kappa = df$kappa[i], phi=0.0001)
X = dplyr::rename(X, covariate = observed)
#ggplot(X, aes(x,y,col=real,fill=real))+ geom_point()

# start with logistic example. latent field (e.g. intrinsic population dynamics dominates spatial variation, but some
# variability due to external driver X)
latent_field = sdmTMB::sim(x = x, y = y, sigma_O = 0.5 , kappa = 4, phi=0.0001)
B0 = 1 #plogis(2+c(-1,0,1))
# don't include covariate effect in null models

X$pred = ifelse(df$oper_model[i]=="covariate",1,0)*0.1*X$covariate + (latent_field$observed+B0)
#ggplot(X, aes(x,y,col=pred,fill=pred))+ geom_powhich(int()
X$y01 = ifelse(runif(nrow(X)) < plogis(X$pred), 1, 0)

# convert to sf object
pa_data <- sf::st_as_sf(X, coords = c("x", "y"))

if(df$cv[i] %in% c("random","systematic")) {
sb <- spatialBlock(speciesData = pa_data,
  species = "y01",
  theRange = df$range[i],
  k = k_folds,
  selection = df$cv[i],
  showBlocks = FALSE)
} else {
  sb = data.frame(foldID = sample(1:k_folds, size=nrow(X),replace=T))
}

X$folds = sb$foldID

# fit the cross-validation model to the data. mesh is recalculated for each fold
#cv = sdmTMB::sdmTMB_cv(form, time = NULL,
#  spatial_only = TRUE,
#  family = binomial(), x = "x", y="y",
#  data=X, fold_ids = "folds", k_folds = max(X$folds), n_knots = 30)

#ggplot(X, aes(x,y,col=as.factor(fold_ids),fill=as.factor(fold_ids)))+ geom_point()

X$pred_cv = NA
X$converged = NA
coef_est = list()
coef_covar = list()
fit_model = list()
sd_covar = list()

for(k in 1:k_folds) {
  sub = dplyr::filter(X, folds!=k)
  #bnd <- INLA::inla.nonconvex.hull(as.matrix(cbind(sub[,c("x","y")])), convex = -0.05)
  #mesh <- INLA::inla.mesh.2d(
  #  boundary = bnd, max.edge = c(3)
  #)
  #mesh = INLA::inla.mesh.create(loc = sub[,c("x","y")])
  #spde <- make_spde(sub$x, sub$y, n_knots = n_knots, mesh = mesh)

  max.edge = diff(range(sub$x))/10
  bound.outer=2
  mesh = INLA::inla.mesh.2d(loc=cbind(sub$x,sub$y),
    max.edge = c(5,5)*max.edge,
    cutoff = max.edge,
    offset = c(max.edge, bound.outer))
  spde <- make_spde(sub$x, sub$y, n_knots = 0, mesh = mesh)

  form = y01 ~ covariate
  if(df$est_model[i]=="null") form = y01 ~ 1
  fit = sdmTMB::sdmTMB(form, time = NULL, spatial_only = TRUE,
    spde = spde, family = binomial(), data=sub)
  # predict to new data
  test_dat = dplyr::filter(X, folds == k)
  pred = predict(fit, newdata=test_dat, xy_cols = c("x","y"))
  X$pred_cv[which(X$folds==k)] = pred$est

  X$converged[which(X$folds==k)] = 1
  eigval <- try(1 / eigen(fit$sd_report$cov.fixed)$values, silent = TRUE) # from sdmtmb
  if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
    X$converged[which(X$folds==k)] = 0
  }
  if(!fit$sd_report$pdHess) {
    X$converged[which(X$folds==k)] = 0
  }

  coef_est[[k]] = fit$sd_report$par.fixed
  coef_covar[[k]] = fit$sd_report$cov.fixed
  sd_covar[[k]] = fit$sd_report$sd
  fit_model[[k]] = fit$model
}
print(paste0(i, ": ",sum(unlist(lapply(fit_model,getElement,3)))))
# save coefficients
save(fit_model, coef_est,coef_covar, sd_covar,file=paste0("sims/example_1/binomial_coefs_",i,".Rdata"))
saveRDS(X, paste0("sims/example_1/binomial_pred_",i,".rds"))

}
