functions {
    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
    
  matrix gp(int N_rows, int N_columns, real[] rows_idx, real[] columns_index,
            real delta0,
            real alpha_gp,
            real rho_gp1, real rho_gp2,
            matrix z1)
  {
    
    matrix[N_rows,N_columns] GP;
    
    matrix[N_rows, N_rows] K1;
    matrix[N_rows, N_rows] L_K1;
    
    matrix[N_columns, N_columns] K2;
    matrix[N_columns, N_columns] L_K2;
    
    K1 = cov_exp_quad(rows_idx, sqrt(alpha_gp), rho_gp1) + diag_matrix(rep_vector(delta0, N_rows));
    K2 = cov_exp_quad(columns_index, sqrt(alpha_gp), rho_gp2) + diag_matrix(rep_vector(delta0, N_columns));

    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    
    GP = kron_mvprod(L_K2, L_K1, z1);

    return(GP);
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
	
	//splines
	int A;
  int num_basis_rows;
  int num_basis_columns;
  matrix[num_basis_rows, A] BASIS_ROWS; 
  matrix[num_basis_columns, A] BASIS_COLUMNS; 
  
  // GP
  real IDX_BASIS_ROWS[num_basis_rows];
  real IDX_BASIS_COLUMNS[num_basis_columns];

}

transformed data
{   
    real delta0 = 1e-9;  
}

parameters {
  vector<lower=0>[D] rho[4];
  vector<lower=0>[4] alpha;
  vector[6] mu;
  real baseline;
  matrix[num_basis_rows,num_basis_columns] z1[4];
  vector<lower=0,upper=1>[N_xi] xi;
}

transformed parameters {
  vector[N_per_group] lxi_pair[N_group];
  vector[N_per_group] llbd[N_group];  
  vector[N_group] muv;
  vector[N_per_group] f_mf_h;
  vector[N_per_group] f_mf_l;
  vector[N_per_group] f_fm_h;
  vector[N_per_group] f_fm_l;
  vector[N_per_group] f[N_group];
  matrix[num_basis_rows,num_basis_columns] beta[4]; 

  muv = rep_vector(0,N_group);
  muv[id_mf] = muv[id_mf] + mu[1];
  muv[id_hh] = muv[id_hh] + mu[2];
  muv[id_hl] = muv[id_hl] + mu[3];
  muv[id_lh] = muv[id_lh] + mu[4];
  muv[id_eh] = muv[id_eh] + mu[5];
  muv[id_el] = muv[id_el] + mu[6];

  for(i in 1:4){
    beta[i] = gp(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0,
              alpha[i], rho[i][1], rho[i][2], z1[i]);
  }
  f_mf_h = to_vector(((BASIS_ROWS') * beta[1] * BASIS_COLUMNS)');
  f_mf_l = to_vector(((BASIS_ROWS') * beta[2] * BASIS_COLUMNS)');
  f_fm_h = to_vector(((BASIS_ROWS') * beta[3] * BASIS_COLUMNS)');
  f_fm_l = to_vector(((BASIS_ROWS') * beta[4] * BASIS_COLUMNS)'); 

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
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
      z1[:,i,j] ~ normal(0,1);
    }
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