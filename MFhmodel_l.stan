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
	cholesky_factor_corr[p] L;
	vector[p] am;
	vector[p] af;
	// Hyperparameters
	real<lower = 1, upper = 4> mu_0;
	real<lower = 0> s2_0;
	real<lower = 0> g;
	real<lower = 0, upper = 10> l;
	real a_0;
	real<lower = 0> t2_0;
}

model {
	matrix[p, p] W = diag_matrix(w);
	// Priors
	mum ~ normal(mu_0, s2_0);
	muf ~ normal(mu_0, s2_0);
	w ~ cauchy(0, g);
	l ~ cauchy(1, 10);
	L ~ lkj_corr_cholesky(l);
	am ~ normal(a_0, t2_0);
	af ~ normal(a_0, t2_0);
	// Hyperpriors
	mu_0 ~ normal(2.5, 2.25)T[1, 4];
	s2_0 ~ cauchy(0, 0.5);
	g ~ cauchy(0, 0.5);
	a_0 ~ normal(0, 10);
	t2_0 ~ cauchy(0, 5);
	// Likelihood: M part
	target += multi_normal_cholesky_lpdf(Xm | mum, W*L);
	for (i in 1:Nm) {
		target += normal_lcdf(dot_product(am, W\(Xm[i] - mum)) | 0, 1);
	}
	// Likelihood: F part
	target += multi_normal_cholesky_lpdf(Xf | muf, W*L);
	for (i in 1:Nf) {
		target += normal_lcdf(dot_product(af, W\(Xf[i] - muf)) | 0, 1);
	}
}

generated quantities {
	vector[Nm+Nf] lik;
	vector[Nm] likM;
	vector[Nf] likF;
	corr_matrix[p] Sigma;
	Sigma = multiply_lower_tri_self_transpose(L);
	for (i in 1:Nm) {
		likM[i] = exp(log(2) + multi_normal_cholesky_lpdf(Xm[i] | mum, diag_matrix(w)*L) + normal_lcdf(dot_product(am, diag_matrix(w)\(Xm[i] - mum)) | 0, 1));
	}
	for (i in 1:Nf) {
		likF[i] = exp(log(2) + multi_normal_cholesky_lpdf(Xf[i] | muf, diag_matrix(w)*L) + normal_lcdf(dot_product(af, diag_matrix(w)\(Xf[i] - muf)) | 0, 1));
	}
	lik = append_row(likM, likF);
}
