rm(list = ls());
# Loading and preparing data------------------
library(rstan);
library(bayesplot);
library(ggplot2);
library(plotrix);
library(emulator);
library(sn);
source("functions.r");
load("data_noNA.rda");
X <-data[,22:41];
m_indx <-which(data$sex=="M");
f_indx <-which(data$sex=="F");
Xm <-X[m_indx, ];
Xf <-X[f_indx, ];
p <-dim(Xm)[2];
# Initialization------------------------------
# dat <-list(Nm = dim(Xm)[1], Nf = dim(Xf)[1], p = p, Xm = Xm, Xf = Xf);
# init <-list(list(s2_0 = 0.1, g = 0.1, w = rep(1, times = p),
#                  mum = rep(2.5, times = p), muf = rep(2.5, times = p), mu_0 = 2.5,
#                  am = rep(0, times = p), af = rep(0, times = p), c = 0.5));
# Fitting in stan
# fit <- stan(file = "MFhmodel_mix.stan",
#             data = dat, init = init,
#             iter = 20, thin = 2, warmup = 10, chains = 1,
#             pars = c("L"), include = FALSE,
#             algorithm = "NUTS",
#             control = list(adapt_delta = 0.90));
# save(fit, file = "MFhmodel.rda");
load("MFhmodel_true.rda");
fit1 <-fit;
load("MFhmodel_trunc1.rda");
# Comparing models---------------------------
LPLM_1 <-LPML(fit1);
WAIC_1 <-WAIC(fit1);
LPLM_2 <-LPML(fit);
WAIC_2 <-WAIC(fit);
# Monte Carlo estimate of "goodness of truncation"
fit <-fit1;
mum_hat <-as.vector(get_posterior_mean(fit, par = "mum"));
muf_hat <-as.vector(get_posterior_mean(fit, par = "muf"));
am_hat <-as.vector(get_posterior_mean(fit, par = "am"));
af_hat <-as.vector(get_posterior_mean(fit, par = "af"));
w_hat <-as.vector(get_posterior_mean(fit, par = "w"));
R_hat <-matrix(get_posterior_mean(fit, par = "Sigma"), nrow = 20, ncol = 20);
S_hat <-diag(as.numeric(w_hat))%*%R_hat%*%diag(as.numeric(w_hat));
likm <-rmsn(n = 10000, xi = mum_hat, Omega = S_hat, alpha = am_hat);
likf <-rmsn(n = 10000, xi = muf_hat, Omega = S_hat, alpha = af_hat);
goftm <-mtc_cdf(likm); # TMale data
goftf <-mtc_cdf(likf); # Female data
# The first element of got is probability < 1 while the second probability > 4
# QUA ABBIAMO SCELTO DI ANDARE AVANTI NELLA ANALISI CON IL MODELLO PIU SEMPLICE NEL FILE "MFhmodel_trunc.R"
# Tuttavia la parte successiva di inferenza funziona


# Convergence assessment---------------------
fit <-fit_mix;
var_names <-colnames(data)[22:41];
mum <-extract(fit, pars = "mum");
colnames(mum$mum) <-var_names;
muf <-extract(fit, pars = "muf");
colnames(muf$muf) <-var_names;
am <-extract(fit, pars = "am");
colnames(am$am) <-var_names;
af <-extract(fit, pars = "af");
colnames(af$af) <-var_names;
w <-extract(fit, pars = "w");
colnames(w$w) <-var_names;
N <-dim(mum$mum)[1];
R <-extract(fit, pars = "Sigma");
# Traceplots
x11();
mcmc_trace(mum);
x11();
mcmc_trace(muf);
x11();
mcmc_trace(am);
x11();
mcmc_trace(af);
x11();
mcmc_trace(w);
# Autocorrelation plots
x11();
mcmc_acf(mum);
x11();
mcmc_acf(muf);
x11();
mcmc_acf(am);
x11();
mcmc_acf(af);
x11();
mcmc_acf(w);

# Inference ----------------------------------------------------
# Posterior of Means
Mm <-Mf <-matrix(rep(0, times = 20*N), nrow = N, ncol = 20);
for (i in 1:N) {
  wi <-diag(w$w[i, ]);
  Ri <-matrix(R$Sigma[i, ,], nrow = 20, ncol = 20);
  Mm[i, ] <-t(skMean(mum$mum[i, ], wi, am$am[i, ], Ri));
  Mf[i, ] <-t(skMean(muf$muf[i, ], wi, af$af[i, ], Ri));
}
M_ <-vector();
M_q <-vector();
M_fill <-vector();
for (i in 1:20) {
  M_ <-c(M_, Mm[ , i]);
  M_q <-c(M_q, rep(var_names[i], times = length(Mm[, i])));
  M_fill <-c(M_fill, rep("M", times = length(Mf[, i])));
}
for (i in 1:20) {
  M_ <-c(M_, Mf[ , i]);
  M_q <-c(M_q, rep(var_names[i], times = length(Mf[, i])));
  M_fill <-c(M_fill, rep("F", times = length(Mf[, i])));
}
dat <-data.frame(Mean = M_, quest = M_q, sex = M_fill);
x11();
ggplot(dat, aes(x = Mean, color = quest, fill = sex)) + geom_density(size = 1, alpha = 0.4);
# Correlation matrix
R_post <-matrix(get_posterior_mean(fit, par = "Sigma"), nrow = 20, ncol = 20);
R_post[lower.tri(R_post, diag = FALSE)] <-0;
x11();
color2D.matplot(R_post, cs1 = c(0, 1), cs2 = c(0, 1), cs3 = c(0, 1), axes = TRUE,
                show.values = TRUE, show.legend = FALSE, xlab = "mq", ylab = "mq");
# Posterior densities of standard deviation of each question
w_ <-vector();
w_q <-vector();
for (i in 1:20) {
  w_ <-c(w_, w$w[, i]);
  w_q <-c(w_q, rep(var_names[i], times = length(w$w[, i])));
}
dat <-data.frame(w = w_, quest = w_q);
x11();
ggplot(dat, aes(x = w, fill = quest)) + geom_density(color = "black", size = 1, alpha = 0.4);
# Hypothesis testing-differences MF: H0: mum-muf<=0 vs H1: mum-muf>0 THIS ARE THE LOCATION PARAMETERS, NOT THE MEANS
mu_diff <-data.frame(mum$mum - muf$muf); colnames(mu_diff) <-var_names;
mq <-mu_diff$mq12; # Change here the variable
dat <- with(density(mq), data.frame(mq = x, d = y));
x11();
pmu <-ggplot(data = dat, mapping = aes(x = mq, y = d)) +
  geom_line(color = "black", size = 1) +
  geom_area(mapping = aes(x = ifelse(mq <= 0, mq, 0)), fill = "magenta", alpha = 0.5) +
  geom_area(mapping = aes(x = ifelse(mq > 0, mq, 0)), fill = "blue", alpha = 0.7) + coord_cartesian(ylim = c(0, 25)) +
  annotate("text", x = mean(mq), y = 25, label = paste("2logBF01 = ", toString(BF01sp(mq, 0))));
pmu
BFmu <-data.frame(matrix(rep(0, times = 20), nrow = 1, ncol = 20));
colnames(BFmu) <-var_names;
for(i in 1:20) {
  BFmu[i] <-BF01sp(mu_diff[, i], 0);
}
# Hypothesis testing H0: am <= 0 vs H1: am > 0
a_m <-data.frame(am$am); colnames(am) <-var_names;
a <-a_m$mq1; # Change here the variable
dat <- with(density(a), data.frame(a = x, d = y));
x11();
pam <-ggplot(data = dat, mapping = aes(x = a, y = d)) +
  geom_line(color = "black", size = 1) +
  geom_area(mapping = aes(x = ifelse(a <= 0, a, 0)), fill = "green", alpha = 0.6) +
  geom_area(mapping = aes(x = ifelse(a > 0, a, 0)), fill = "black", alpha = 0.4) + coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = mean(a), y = 1, label = paste("2logBF01 = ", toString(BF01sp(a, 0))));
pam
BFam <-data.frame(matrix(rep(0, times = 20), nrow = 1, ncol = 20));
colnames(BFam) <-var_names;
for(i in 1:20) {
  BFam[i] <-BF01sp(am[, i], 0);
}
# Hypothesis testing H0: af <= 0 vs H1: af > 0
a_f <-data.frame(af$af); colnames(af) <-var_names;
a <-a_f$mq20; # Change here the variable
dat <- with(density(a), data.frame(a = x, d = y));
x11();
paf <-ggplot(data = dat, mapping = aes(x = a, y = d)) +
  geom_line(color = "black", size = 1) +
  geom_area(mapping = aes(x = ifelse(a <= 0, a, 0)), fill = "green", alpha = 0.6) +
  geom_area(mapping = aes(x = ifelse(a > 0, a, 0)), fill = "black", alpha = 0.4) + coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = mean(a), y = 1, label = paste("2logBF01 = ", toString(BF01sp(a, 0))));
paf
BFaf <-data.frame(matrix(rep(0, times = 20), nrow = 1, ncol = 20));
colnames(BFaf) <-var_names;
for(i in 1:20) {
  BFaf[i] <-BF01sp(af[, i], 0);
}
# Hypothesis testing differences MF: H0: |am|-|af|<=0 vs H1: |am|-|af|>0
a_diff <-data.frame(abs(am$am) - abs(af$af)); colnames(a_diff) <-var_names;
a <-a_diff$mq1; # Change here the variable
dat <- with(density(a), data.frame(a = x, d = y));
x11();
padiff <-ggplot(data = dat, mapping = aes(x = a, y = d)) +
  geom_line(color = "black", size = 1) +
  geom_area(mapping = aes(x = ifelse(a <= 0, a, 0)), fill = "magenta", alpha = 0.5) +
  geom_area(mapping = aes(x = ifelse(a > 0, a, 0)), fill = "blue", alpha = 0.7) + coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = mean(a), y = 1, label = paste("2logBF01 = ", toString(BF01sp(a, 0))));
padiff
BFadiff <-data.frame(matrix(rep(0, times = 20), nrow = 1, ncol = 20));
colnames(BFadiff) <-var_names;
for(i in 1:20) {
  BFadiff[i] <-BF01sp(a_diff[, i], 0);
}
# Hypothesis testing differences MF: H0: Mm - Mf<=0 vs H1: Mm - Mf>0
M_diff <-data.frame(Mm - Mf); colnames(M_diff) <-var_names;
M <-M_diff$mq20; # Change here the variable
dat <- with(density(M), data.frame(M = x, d = y));
x11();
pMdiff <-ggplot(data = dat, mapping = aes(x = M, y = d)) +
  geom_line(color = "black", size = 1) +
  geom_area(mapping = aes(x = ifelse(M <= 0, M, 0)), fill = "magenta", alpha = 0.5) +
  geom_area(mapping = aes(x = ifelse(M > 0, M, 0)), fill = "blue", alpha = 0.7) + coord_cartesian(ylim = c(0, 20)) +
  annotate("text", x = mean(M), y = 20, label = paste("2logBF01 = ", toString(BF01sp(M, 0))));
pMdiff
BFMdiff <-data.frame(matrix(rep(0, times = 20), nrow = 1, ncol = 20));
colnames(BFMdiff) <-var_names;
for(i in 1:20) {
  BFMdiff[i] <-BF01sp(M_diff[, i], 0);
}
