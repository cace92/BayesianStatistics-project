library(ggplot2);
library(sn);
library(truncnorm);

n <-1000;
p <-c(0.1, 0.25, 0.5, 0.75, 0.9);
np <-length(p);
mu <-2;
alpha <--2;
sigma <-2;
omega <-1;

delta <-1/sqrt(1+alpha^2*omega)*omega*alpha;
chi <-mu-sqrt(2/pi)*sigma*delta;

x <-matrix(rep(0, times = n*np), nrow = n, ncol = np);
for (i in 1:np) {
  for (j in 1:n) {
    if (runif(1, 0, 1) < p[i]) {
      x[j, i] <-rmsn(n = 1, xi = chi, Omega = omega, alpha = alpha);
    } else {
      x[j, i] <-rtruncnorm(1, a=1, b=4, mean = mu, sd = sigma);
    }
  }
}
out_prob <-rep(0, times = np);
for (i in 1:np) {
  out_prob[i] <-(length(which(x[,i] < 1))+length(which(x[,i] > 4))) / n;
}
x_vec <-vector();
x_lab <-vector();
lab <-strwrap(p);
for (i in 1:np) {
  x_vec <-c(x_vec, x[, i]);
  x_lab <-c(x_lab, lab[i]);
}
xd <-data.frame(data = x_vec, label = x_lab);
x <-data.frame(x); colnames(x) <-lab;
x11();
ggplot(xd, aes(x = data, color = label)) + geom_density(size = 1) + xlim(-2, 6);

x11();
ggplot(x, aes(x = x[,1], color = "black")) + geom_density(size = 1)