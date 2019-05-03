data{
    int<lower=1> N;
    int SUC[N];
    int TRIAL[N];
    int GROUP[N];
    int MALE[N];
}
parameters{
    real a;
    real grouptwo;
    real male;
}
model{
    vector[N] p_suc;
    male ~ normal( 0 , 10 );
    grouptwo ~ normal( 0 , 10 );
    a ~ normal( 0 , 100 );
    for ( i in 1:N ) {
        p_suc[i] = a + grouptwo * GROUP[i] + male * MALE[i];
    }
    SUC ~ binomial_logit( TRIAL , p_suc );
}

