functions {
	real mix_lik_lpdf(vector[] X, int N, int c, vector mu, vector a, matrix W, matrix L, real p) {
		real log_lik = 0;
		real skew_log_lik = 0;
		real trunc_log_lik = 0;
		vector[c] chi;
		vector[c] delta;
		matrix[c, c] S;
		matrix[c, c] Omega;
		vector[c] temp;
		Omega = multiply_lower_tri_self_transpose(L);
		delta = (1/sqrt(1+quad_form_sym(Omega, a)))*Omega*a;
		chi = mu + sqrt(2/pi())*W*delta;
		S = W*(Omega-(2/pi())*(delta*delta'))*W;
		for (i in 1:N) {
			trunc_log_lik = 0;
			skew_log_lik = multi_normal_cholesky_lpdf(X[i] | mu, W*L) + normal_lcdf(dot_product(a, W\(X[i] - mu)) | 0, 1);
			temp = X[i];
			for (j in 1:c) {
				trunc_log_lik += normal_lpdf(temp[j] | chi[j], sqrt(S[j, j])) -
								 log(sqrt(S[j, j])) - log(normal_cdf(4, chi[j], sqrt(S[j, j])) - normal_cdf(1, chi[j], sqrt(S[j, j])));
			}
			log_lik += log_mix(p, skew_log_lik, trunc_log_lik);
		}
		return log_lik;
	}
	
	real mix_lik_un_lpdf(vector X, int c, vector mu, vector a, matrix W, matrix L, real p) {
		real log_lik;
		real skew_log_lik;
		real trunc_log_lik;
		vector[c] chi;
		vector[c] delta;
		matrix[c, c] S;
		matrix[c, c] Omega;
		Omega = multiply_lower_tri_self_transpose(L);
		delta = (1/sqrt(1+quad_form_sym(Omega, a)))*Omega*a;
		chi = mu + sqrt(2/pi())*W*delta;
		S = W*(Omega-(2/pi())*(delta*delta'))*W;
		trunc_log_lik = 0;
		skew_log_lik = multi_normal_cholesky_lpdf(X | mu, W*L) + normal_lcdf(dot_product(a, W\(X - mu)) | 0, 1);
		for (j in 1:c) {
			trunc_log_lik += normal_lpdf(X | chi[j], sqrt(S[j, j])) -
							 log(sqrt(S[j, j])) - log(normal_cdf(4, chi[j], sqrt(S[j, j])) - normal_cdf(1, chi[j], sqrt(S[j, j])));
		}
		log_lik = log_mix(p, skew_log_lik, trunc_log_lik);
		return log_lik;
	}
}

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
	real<lower = 0, upper = 1> c;
	// Hyperparameters
	real<lower = 1, upper = 4> mu_0;
	real<lower = 0> s2_0;
	real<lower = 0> g;
	real a_0;
	real<lower = 0> t2_0;
}

model {
	matrix[p, p] W = diag_matrix(w);
	// Priors
	mum ~ normal(mu_0, s2_0);
	muf ~ normal(mu_0, s2_0);
	w ~ cauchy(0, g);
	L ~ lkj_corr_cholesky(1);
	am ~ normal(a_0, t2_0);
	af ~ normal(a_0, t2_0);
	c ~ beta(1, 3);
	// Hyperpriors
	mu_0 ~ normal(2.5, 2.25)T[1, 4];
	s2_0 ~ cauchy(0, 0.5);
	g ~ cauchy(0, 0.5);
	a_0 ~ normal(0, 10);
	t2_0 ~ cauchy(0, 5);
	// Likelihood: M part
	target += mix_lik_lpdf(Xm | Nm, p, mum, am, W, L, c);
	// Likelihood: F part
	target += mix_lik_lpdf(Xf | Nf, p, muf, af, W, L, c);
}

generated quantities {
	vector[Nm+Nf] lik;
	vector[Nm] likM;
	vector[Nf] likF;
	matrix[p, p] W = diag_matrix(w);
	for (i in 1:Nm) {
		likM[i] = exp(mix_lik_un_lpdf(Xm[i] | p, mum, am, W, L, c));
	}
	for (i in 1:Nf) {
		likF[i] = exp(mix_lik_un_lpdf(Xf[i] | p, muf, af, W, L, c));
	}
	lik = append_row(likM, likF);
}