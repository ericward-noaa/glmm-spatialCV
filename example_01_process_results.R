library(glmmfields)
library(blockCV)
library(tidyverse)
library(sdmTMB)
library(ggplot2)

# drop cases where CV is random wrt space, but range is != 5000
df = readRDS(df, file="sims/example_1/grid.rds")

df$LL = df$est_int = df$est_cov = df$est_tau_O = df$est_kappa = NA

for(i in 1:nrow(df)) {

# load log-likelihood, only retaining saved models
ll = readRDS(file = paste0("sims/example_1/binomial_pred_",i,".rds"))
load(file=paste0("sims/example_1/binomial_coefs_",i,".Rdata"))

if(sum(unlist(lapply(fit_model, getElement, 3)))==0) {
  #df$LL[i] = sum(dbinom(ll$y01, size=1, prob=plogis(ll$pred_cv),log=TRUE))
  load(file=paste0("sims/example_1/binomial_coefs_",i,".Rdata"))
  df$LL[i] = mean(unlist(lapply(fit_model,getElement,2)))
  load(file=paste0("sims/example_1/binomial_coefs_",i,".Rdata"))
  df$est_int[i] = mean(unlist(lapply(coef_est,getElement,1)))
  indx = 0
  if(df$oper_model[i]=="covariate") {
    df$est_cov[i] = mean(unlist(lapply(coef_est,getElement,2)))
    indx = 1
  }
  df$est_tau_O[i] = mean(exp(unlist(lapply(coef_est,getElement,indx+1))))
  df$est_kappa[i] = mean(exp(unlist(lapply(coef_est,getElement,indx+2))))
}

}

# facets here will by cv (fixed range, kappa)
# first facet here is looking at effect of cv sampling
df_1 = dplyr::filter(df, range==unique(df$range)[2],kappa==2) %>%
  dplyr::group_by(sim, kappa, range, cv) %>%
  dplyr::summarize(ll_cov_cov = (LL[which(oper_model=="covariate" & est_model=="covariate")]),
    ll_null_cov = (LL[which(oper_model=="null" & est_model=="covariate")]),
    ll_cov_null = (LL[which(oper_model=="covariate" & est_model=="null")]),
    ll_null_null = (LL[which(oper_model=="null" & est_model=="null")]),
    est_cov_cov = mean(est_cov[which(oper_model=="covariate" & est_model=="covariate")],na.rm=T),
    est_cov_null = mean(LL[which(oper_model=="covariate" & est_model=="null")],na.rm=T))


# compare likelihoods when true model includes covariate. negative values mean null model
# gets more support, and that support for the null increases as we do spatial block sampling
ggplot(df_1, aes(cv, ll_cov_cov - ll_cov_null)) + geom_boxplot()

# look at range effect
df_2 = dplyr::filter(df, cv=="random", model=="covariate",kappa==2) %>%
  dplyr::group_by(sim, kappa, range, cv) %>%
  dplyr::summarize(ll_cov_cov = mean(LL[which(covariate==TRUE & model=="covariate")],na.rm=T),
    ll_null_cov = mean(LL[which(covariate==FALSE & model=="covariate")],na.rm=T),
    ll_cov_null = mean(LL[which(covariate==TRUE & model=="null")],na.rm=T),
    ll_null_null = mean(LL[which(covariate==FALSE & model=="null")],na.rm=T),
    est_cov_cov = mean(est_cov[which(covariate==TRUE & model=="covariate")],na.rm=T),
    est_cov_null = mean(LL[which(covariate==TRUE & model=="null")],na.rm=T))

# compare likelihoods when true model includes covariate. negative values mean null model
# gets more support, and that support for the null increases as we do spatial block sampling
ggplot(df_2, aes(range, ll_cov_cov - ll_null_cov)) + geom_boxplot()



# look at kappa effect
df_3 = dplyr::filter(df, range==5000, model=="covariate", cv=="random")


df_summary = dplyr::filter(df, cv!="none") %>%
  dplyr::group_by(sim, kappa, range, cv) %>%
  dplyr::summarize(ll_true = LL[which(covariate==TRUE)],
    ll_false = LL[which(covariate==FALSE)],
    est_cov = mean(est_cov,na.rm=T))
df_summary$ll_diff = df_summary$ll_true - df_summary$ll_false

ggplot(df_summary, aes(as.factor(kappa), ll_diff)) + geom_boxplot()
ggplot(df_summary, aes(as.factor(cv), ll_diff)) + geom_boxplot()
ggplot(df_summary, aes(as.factor(range), ll_diff)) + geom_boxplot()

ggplot(df_summary, aes(as.factor(range), est_cov)) + geom_boxplot()
ggplot(df_summary, aes(as.factor(cv), est_cov)) + geom_boxplot()
ggplot(df_summary, aes(as.factor(kappa), est_cov)) + geom_boxplot()


df_summary = dplyr::filter(df, range=="5000") %>%
dplyr::group_by(sim, kappa, range, cv) %>%
  dplyr::summarize(ll_true = LL[which(covariate==TRUE)],
    ll_false = LL[which(covariate==FALSE)],
    est_cov = mean(est_cov,na.rm=T))
df_summary$ll_diff = df_summary$ll_true - df_summary$ll_false


ggplot(df_summary, aes(as.factor(cv), ll_diff)) + geom_boxplot()
ggplot(df_summary, aes(as.factor(cv), est_cov)) + geom_boxplot()
