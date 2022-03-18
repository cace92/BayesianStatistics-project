# Loading and preparing data------------------
library(rstan);
library(bayesplot);
library(ggplot2);
library(plotrix);
library(emulator);
library(sn);
source("functions.r");
load("data_noNA.rda");
load("Chmodel_trunc.rda");
# Convergence assessment---------------------
var_names <-colnames(data)[22:41];
mumi <-extract(fit, pars = "mumi");
colnames(mumi$mumi) <-var_names;
mubv <-extract(fit, pars = "mubv");
colnames(mubv$mubv) <-var_names;
muoth <-extract(fit, pars = "muoth");
colnames(muoth$muoth) <-var_names;
w <-extract(fit, pars = "w");
colnames(w$w) <-var_names;
N <-dim(mumi$mumi)[1];
# Traceplots
x11();
mcmc_trace(mumi);
x11();
mcmc_trace(mubv);
x11();
mcmc_trace(muoth);
x11();
mcmc_trace(w);
# Autocorrelation plots
x11();
mcmc_acf(mumi);
x11();
mcmc_acf(mubv);
x11();
mcmc_acf(muoth);
x11();
mcmc_acf(w);
# Inference----------------------------------------------------------
# Density plot mu
mq <-data.frame(mumi$mumi)$mq20; # Change here the variable
dat <-data.frame(question = mq);
x11();
dmu <-ggplot(data = dat, mapping = aes(x = question)) + 
  geom_density(color = "black", fill = "blue", size = 1, alpha = 0.4);
dmu
# Densities plot w
w_ <-vector();
w_q <-vector();
for (i in 1:20) {
  w_ <-c(w_, w$w[, i]);
  w_q <-c(w_q, rep(var_names[i], times = length(w$w[, i])));
}
dat <-data.frame(k = w_, quest = w_q);
dat$quest <-factor(dat$quest, levels = var_names);
x11();
ggplot(dat, aes(x = k, fill = quest)) + geom_density(color = "black", size = 1, alpha = 0.4);
# Differences of the means
alpha <-0.10;
dat <-data.frame(posterior = rep(0.00, times = 60), question = rep(var_names, times = 3),
                 campus = c(rep("MI", times = 20), rep("BV", times = 20), rep("OTH", times = 20)),
                 lower = rep(0, times = 60), upper = rep(0, times = 60));
for (i in 1:20) {
  dat$posterior[i] <-mean(mumi$mumi[, i]);
  res <-quantile(mumi$mumi[, i], probs = c(alpha/2, 1-alpha/2));
  dat$lower[i] <-res[1];
  dat$upper[i] <-res[2]; 
}
for (i in 21:40) {
  dat$posterior[i] <-mean(mubv$mubv[, i-20]);
  res <-quantile(mubv$mubv[, i-20], probs = c(alpha/2, 1-alpha/2));
  dat$lower[i] <-res[1];
  dat$upper[i] <-res[2]; 
}
for (i in 41:60) {
  dat$posterior[i] <-mean(muoth$muoth[, i-40]);
  res <-quantile(muoth$muoth[, i-40], probs = c(alpha/2, 1-alpha/2));
  dat$lower[i] <-res[1];
  dat$upper[i] <-res[2]; 
}
dat$question <-factor(dat$question, levels = var_names);
dat$campus <-factor(dat$campus, levels = c("MI", "BV", "OTH"));
x11();
pmu <-ggplot(data = dat, mapping = aes(x = question, y = posterior, color = campus)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 1, alpha = 0.6);
pmu

