data{
    int<lower=1> N;
    int SUC[N];
    int TRIAL[N];
    int COMM_NUM_B[N];
    int AGE1[N];
    int AGE2[N];
    int MALE[N];
}
parameters{
    real a;
    vector[2] comm;
    real<lower=0> dispersion;
    real<lower=0> sig_comm;
    real male;
    real age1;
    real age2;
}
model{
    vector[N] p_suc;
    age2 ~ normal( 0 , 10 );
    age1 ~ normal( 0 , 10 );
    male ~ normal( 0 , 10 );
    dispersion ~ exponential( 1 );
    sig_comm ~ cauchy( 0 , 1 );
    comm ~ normal( 0 , sig_comm );
    a ~ normal( 0 , 100 );
    for ( i in 1:N ) {
        p_suc[i] = a + comm[COMM_NUM_B[i]] + male * MALE[i] + age1 * AGE1[i] + age2 * AGE2[i];
        p_suc[i] = inv_logit(p_suc[i]);
    }
    SUC ~ beta_binomial( TRIAL , p_suc*dispersion , (1-p_suc)*dispersion );
}
generated quantities{
    vector[N] p_suc;
    for ( i in 1:N ) {
        p_suc[i] = a + comm[COMM_NUM_B[i]] + male * MALE[i] + age1 * AGE1[i] + age2 * AGE2[i];
        p_suc[i] = inv_logit(p_suc[i]);
    }
}
