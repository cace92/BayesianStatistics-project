data {
	int<lower = 0> Nm;
	int<lower = 0> Nf;
	int<lower = 0> p;
	vector<lower = 1, upper = 4>[p] Xm[Nm];
	vector<lower = 1, upper = 4>[p] Xf[Nf];
}

parameters {
	// Parameters
	vector<lower = 1, upper = 4>[p] mum;
	vector<lower = 1, upper = 4>[p] muf;
	vector<lower = 0>[p] w;
	// Hyperparameters
	real<lower = 1, upper = 4> mu_0;
	real<lower = 0> s_0;
	real<lower = 0> g;
}

model {
	vector[p] temp;
	// Priors
	for (j in 1:p) {
		mum[j] ~ normal(mu_0, s_0)T[1, 4];
		muf[j] ~ normal(mu_0, s_0)T[1, 4];
	}
	w ~ cauchy(0, g);
	// Hyperpriors
	mu_0 ~ normal(2.5, 2.25)T[1, 4];
	s_0 ~ cauchy(0, 5);
	g ~ cauchy(0, 5);
	// Likelihood: M part
	for (i in 1:Nm) {
		temp = Xm[i];
		for (j in 1:p) {
			temp[j] ~ normal(mum[j], w[j])T[1, 4];
		}
	}
	// Likelihood: F part
	for (i in 1:Nf) {
		temp = Xf[i];
		for (j in 1:p) {
			temp[j] ~ normal(muf[j], w[j])T[1, 4];
		}
	}
}

generated quantities {
	vector[p] temp;
	vector[Nm+Nf] lik;
	vector[Nm] likM = rep_vector(1.0, Nm);
	vector[Nf] likF = rep_vector(1.0, Nf);
	for (i in 1:Nm) {
		temp = Xm[i];
		for (j in 1:p) {
			likM[i] *= exp(normal_lpdf(temp[j] | mum[j], w[j]))/(w[j]*(normal_cdf(4, mum[j], w[j]) - normal_cdf(1, mum[j], w[j])));
		}
	}
	for (i in 1:Nf) {
		temp = Xf[i];
		for (j in 1:p) {
			likF[i] *= exp(normal_lpdf(temp[j] | muf[j], w[j]))/(w[j]*(normal_cdf(4, muf[j], w[j]) - normal_cdf(1, muf[j], w[j])));
		}
	}
	lik = append_row(likM, likF);
}
