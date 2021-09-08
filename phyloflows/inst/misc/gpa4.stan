functions {
  vector spd_cov_exp_quad(vector[] x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real constant = square(sdgp) * (sqrt(2 * pi()) * lscale[1])^D;
      real neg_half_lscale2 = -0.5 * square(lscale[1]);
      for (m in 1:NB) {
        out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
      }
    } else {
      // multi-dimensional non-isotropic GP
      real constant = square(sdgp) * sqrt(2 * pi())^D * prod(lscale);
      vector[Dls] neg_half_lscale2 = -0.5 * square(lscale);
      for (m in 1:NB) {
        out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m])));
      }
    }
    return out;
  }

  vector gpa_coefficient(int M_nD, real sdgp, vector lscale, vector zgp, vector[] slambda) { 
    vector[M_nD] diag_spd = sqrt(spd_cov_exp_quad(slambda, sdgp, lscale));
    return diag_spd .* zgp;
  }
}
data {
  int<lower=1> N;
  int<lower=1> N_group;
  int<lower=1> N_per_group;
  int<lower=1> D;
  matrix[N_per_group,D] x;
  int y[N_per_group, N_group];
	
  int<lower=1> M_nD;
  matrix[N_per_group, M_nD] Xgp;
  vector[D] slambda[M_nD];			

  int<lower=1> N_xi;	
  matrix[N_xi,2] shape;
  int xi_id_src[N_per_group,N_group];
  int xi_id_rec[N_per_group,N_group];
  int id_mf_h[6];
  int id_mf_l[6];
  int id_fm_h[6];
  int id_fm_l[6];
  int id_mf[12];
  int id_fm[12];
  int id_hh[4];
  int id_hl[4];	
  int id_lh[4];	
  int id_eh[4];
  int id_el[4];	
	
}

transformed data {
}

parameters {
  vector<lower=0>[D] rho[4];
  vector<lower=0>[4] alpha;
  vector[6] mu;
  real baseline;
  vector[M_nD] beta[4];
  vector<lower=0,upper=1>[N_xi] xi;
}

transformed parameters {
  vector[N_per_group] lxi_pair[N_group];
  vector[M_nD] coefficient[4];  
  vector[N_per_group] llbd[N_group];  
  vector[N_group] muv;
  vector[N_per_group] f_mf_h;
  vector[N_per_group] f_mf_l;
  vector[N_per_group] f_fm_h;
  vector[N_per_group] f_fm_l;
  vector[N_per_group] f[N_group];


  muv = rep_vector(0,N_group);
  muv[id_mf] = muv[id_mf] + mu[1];
  muv[id_hh] = muv[id_hh] + mu[2];
  muv[id_hl] = muv[id_hl] + mu[3];
  muv[id_lh] = muv[id_lh] + mu[4];
  muv[id_eh] = muv[id_eh] + mu[5];
  muv[id_el] = muv[id_el] + mu[6];


  coefficient[1] = gpa_coefficient(M_nD, alpha[1], rho[1], beta[1], slambda);
  f_mf_h = Xgp * coefficient[1];
  coefficient[2] = gpa_coefficient(M_nD, alpha[2], rho[2], beta[2], slambda);
  f_mf_l = Xgp * coefficient[2];
  coefficient[3] = gpa_coefficient(M_nD, alpha[3], rho[3], beta[3], slambda);
  f_fm_h = Xgp * coefficient[3];
  coefficient[4] = gpa_coefficient(M_nD, alpha[4], rho[4], beta[4], slambda);
  f_fm_l = Xgp * coefficient[4];

  for (i in 1:6){
    f[id_mf_h[i]]=f_mf_h;
    f[id_mf_l[i]]=f_mf_l;
    f[id_fm_h[i]]=f_fm_h;
    f[id_fm_l[i]]=f_fm_l;
  }

  for (i in 1:N_group){
    lxi_pair[i] = log(xi[xi_id_src[,i]]) + log(xi[xi_id_rec[,i]]);
    llbd[i] = baseline + muv[i] + f[i] + lxi_pair[i];
  }

}


model {
  for (i in 1:N_xi){
    xi[i] ~ beta(shape[i,1],shape[i,2]);
  }
  for (i in 1:4){
    beta[i] ~ normal(0,1);
  } 
  for (i in 1:2){
    rho[i][1] ~ inv_gamma(2.93854,8.30173);
    rho[i][2] ~ inv_gamma(2.60877,7.73394); 
  }
  for (i in 3:4){
    rho[i][1] ~ inv_gamma(3.48681,9.21604);
    rho[i][2] ~ inv_gamma(2.60877,7.73394);
  }

  
  alpha ~ normal(0, 10);
  mu ~ normal(0, 10);
  baseline ~ normal(0, 10);

  for (i in 1:N_group){
    y[,i] ~ poisson_log(llbd[i]);
  }
}

generated quantities{
  int y_predict[N_per_group,N_group];
  real log_lik[N_per_group,N_group];
  for(i in 1:N_group){
    for(j in 1:N_per_group){
      y_predict[j,i] = poisson_log_rng(llbd[i][j]);
      log_lik[j,i] = poisson_log_lpmf(y[j,i]| llbd[i][j]);
    }
  }
}


