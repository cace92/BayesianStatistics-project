rm(list = ls());
# Loading and preparing data------------------
library(rstan);
source("functions.r");
load("data_noNA.rda");
X <-data[,22:41];
m_indx <-which(data$sex=="M");
f_indx <-which(data$sex=="F");
Xm <-X[m_indx, ];
Xf <-X[f_indx, ];
p <-dim(Xm)[2];
# Initialization------------------------------
dat <-list(Nm = dim(Xm)[1], Nf = dim(Xf)[1], p = p, Xm = Xm, Xf = Xf);
init <-list(list(s2_0 = 0.1, g = 0.1, w = rep(1, times = p),
                 mum = rep(2.5, times = p), muf = rep(2.5, times = p), mu_0 = 2.5,
                 am = rep(0, times = p), af = rep(0, times = p), c = 0.5));
# Fitting in stan
fit <- stan(file = "MFhmodel_true.stan",
            data = dat, init = init,
            iter = 20, thin = 2, warmup = 10, chains = 1,
            pars = c("L"), include = FALSE,
            algorithm = "NUTS",
            control = list(adapt_delta = 0.90));
save(fit, file = "MFhmodel_true.rda")
system("shutdown -s")