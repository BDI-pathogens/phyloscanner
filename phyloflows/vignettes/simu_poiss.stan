functions {
  matrix L_cov_exp_quad_ARD(vector[] x,
                            real alpha,
                            vector rho) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + 1e-10;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha
                      * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + 1e-10;
    return cholesky_decompose(K);
  }
}
data {
  int<lower=1> N;
  int<lower=1> D;
  vector[D] x[N];  
 
  vector<lower=0>[D] rho;
  real alpha;
  real mu;

  vector[N] xi;


}

transformed data {
  matrix[N, N] L_cov = L_cov_exp_quad_ARD(x, alpha, rho);
}

parameters {}
model {}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
  vector[N] y;
  for (n in 1:N)
    y[n] = poisson_log_rng((f[n]+mu)+log(xi[n]));
}
