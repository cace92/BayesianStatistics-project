data {
	int<lower = 0> Nmi;
	int<lower = 0> Nbv;
	int<lower = 0> Noth;
	int<lower = 0> p;
	vector<lower = 1, upper = 4>[p] Xmi[Nmi];
	vector<lower = 1, upper = 4>[p] Xbv[Nbv];
	vector<lower = 1, upper = 4>[p] Xoth[Noth];
}

parameters {
	// Parameters
	vector<lower = 1, upper = 4>[p] mumi;
	vector<lower = 1, upper = 4>[p] mubv;
	vector<lower = 1, upper = 4>[p] muoth;
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
		mumi[j] ~ normal(mu_0, s_0)T[1, 4];
		mubv[j] ~ normal(mu_0, s_0)T[1, 4];
		muoth[j] ~ normal(mu_0, s_0)T[1, 4];
	}
	w ~ cauchy(0, g);
	// Hyperpriors
	mu_0 ~ normal(2.5, 2.25)T[1, 4];
	s_0 ~ cauchy(0, 5);
	g ~ cauchy(0, 5);
	// Likelihood: MI part
	for (i in 1:Nmi) {
		temp = Xmi[i];
		for (j in 1:p) {
			temp[j] ~ normal(mumi[j], w[j])T[1, 4];
		}
	}
	// Likelihood: BV part
	for (i in 1:Nbv) {
		temp = Xbv[i];
		for (j in 1:p) {
			temp[j] ~ normal(mubv[j], w[j])T[1, 4];
		}
	}
	// Likelihood: OTHER part
	for (i in 1:Noth) {
		temp = Xoth[i];
		for (j in 1:p) {
			temp[j] ~ normal(muoth[j], w[j])T[1, 4];
		}
	}
}
