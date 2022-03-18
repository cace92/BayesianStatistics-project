rm(list = ls());
# Loading and preparing data------------------
library(rstan);
library(bayesplot);
library(ggplot2);
source("functions.r");
load("data_noNA.rda");
X <-data[,22:41];
m_indx <-which(data$lev=="1 liv.");
f_indx <-which(data$lev=="Mag.");
Xm <-X[m_indx, ];
Xf <-X[f_indx, ];
p <-dim(Xm)[2];
load("BMhmodel_trunc.rda");
# !!! Since the model is the same of the previous Male-Female, to better understand variables' names Mlae<->Bachelor and Female<->Master
# Convergence assessment---------------------
var_names <-colnames(data)[22:41];
mum <-extract(fit, pars = "mum");
colnames(mum$mum) <-var_names;
muf <-extract(fit, pars = "muf");
colnames(muf$muf) <-var_names;
w <-extract(fit, pars = "w");
colnames(w$w) <-var_names;
N <-dim(mum$mum)[1];
# Traceplots
x11();
mcmc_trace(mum);
x11();
mcmc_trace(muf);
x11();
mcmc_trace(w);
# Autocorrelation plots
x11();
mcmc_acf(mum);
x11();
mcmc_acf(muf);
x11();
mcmc_acf(w);
# Inference----------------------------------------------------------
# Posterior density plot mu
mq <-data.frame(mum$mum)$mq20; # Change here the variable
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
# Hypothesis testing-differences MF: H0: mum-muf<=0 vs H1: mum-muf>0 read Bachelor vs Master
alpha <-0.10;
mu_diff <-data.frame(mum$mum - muf$muf); colnames(mu_diff) <-var_names;
dat <-data.frame(postrior = rep(0.00, times = 20), question = var_names,
                 lower = rep(0, times = 20), upper = rep(0, times = 20));
for (i in 1:20) {
  dat$posterior[i] <-mean(mu_diff[, i]);
  res <-quantile(mu_diff[, i], probs = c(alpha/2, 1-alpha/2));
  dat$lower[i] <-res[1];
  dat$upper[i] <-res[2]; 
}
dat$question <-factor(dat$question, levels = dat$question);
x11();
pmu <-ggplot(data = dat, mapping = aes(x = question, y = posterior, color = question)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 1) + geom_hline(aes(yintercept = 0));
pmu
# Posterior density plot of mu difference
mq <-mu_diff$mq20; # Change here the variable
dat <- with(density(mq), data.frame(mq = x, d = y));
x11();
dmudiff <-ggplot(data = dat, mapping = aes(x = mq, y = d)) +
  geom_line(color = "black", size = 1) +
  geom_area(mapping = aes(x = ifelse(mq <= 0, mq, 0)), fill = "magenta", alpha = 0.5) +
  geom_area(mapping = aes(x = ifelse(mq > 0, mq, 0)), fill = "blue", alpha = 0.7) + coord_cartesian(ylim = c(0, 25));
dmudiff