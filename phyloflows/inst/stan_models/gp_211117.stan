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
  int<lower=1> N_group; // number of directions
  int<lower=1> N_per_group; // age-age entries
  int y[N_per_group, N_group]; // count of transmissions for each age-age entry
  int<lower=0,upper=1> is_mf[N_group]; // is the direction mf
	
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
  real<lower=0> rho_gp1[N_group];
  real<lower=0> rho_gp2[N_group];
  real<lower=0> alpha_gp[N_group];
  real nu;
  real mu;
  matrix[num_basis_rows,num_basis_columns] z1[N_group];
}

transformed parameters {
  vector[N_per_group] log_lambda[N_group];  
  vector[N_per_group] f[N_group];
  matrix[num_basis_rows,num_basis_columns] beta[N_group]; 

  for(i in 1:N_group){
    beta[i] = gp(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0,
              alpha_gp[i], rho_gp1[i], rho_gp2[i], z1[i]);
    f[i] = to_vector(((BASIS_ROWS') * beta[i] * BASIS_COLUMNS)');

    log_lambda[i] = mu + f[i];
    
    if(is_mf[i])
      log_lambda[i] += nu;
  }
  

}

model {
  alpha_gp ~ cauchy(0,1);
  rho_gp1 ~ inv_gamma(5, 5);
  rho_gp2 ~ inv_gamma(5, 5);
  mu ~ normal(0, 10);
  nu ~ normal(0, 10);
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
      z1[:,i,j] ~ normal(0,1);
    }
  }

  for (i in 1:N_group){
    y[,i] ~ poisson_log(log_lambda[i]);
  }
}

generated quantities{
  int y_predict[N_per_group,N_group];
  real log_lik[N_per_group,N_group];
  for(i in 1:N_group){
    for(j in 1:N_per_group){
      y_predict[j,i] = poisson_log_rng(log_lambda[i][j]);
      log_lik[j,i] = poisson_log_lpmf(y[j,i]| log_lambda[i][j]);
    }
  }
}