data{    int<lower=1> N;    int SUC[N];    int TRIAL[N];    int AGE1[N];    int AGE2[N];    real INMIGRANT[N];    int MALE[N];    int COMM_TYPE_F[N];}parameters{    real a;    real<lower=0> dispersion;    real fishing;
    real male_age1_inmigrant; 
    real male_age2_inmigrant; 
    real male_age3_inmigrant; 
    real female_age1_inmigrant; 
    real female_age2_inmigrant; 
    real female_age3_inmigrant;
    real male_age1_resident; 
    real male_age2_resident; 
    real male_age3_resident; 
    real female_age1_resident; 
    real female_age2_resident;}model{    vector[N] p_part;
    dispersion ~ exponential( 1 );    a ~ normal( 0 , 100 );
    fishing ~ normal( 0 , 10 );
    male_age1_inmigrant ~ normal( 0 , 10 );
    male_age2_inmigrant ~ normal( 0 , 10 ); 
    male_age3_inmigrant ~ normal( 0 , 10 );
    female_age1_inmigrant ~ normal( 0 , 10 );
    female_age2_inmigrant ~ normal( 0 , 10 ); 
    female_age3_inmigrant ~ normal( 0 , 10 );
    male_age1_resident ~ normal( 0 , 10 ); 
    male_age2_resident ~ normal( 0 , 10 );
    male_age3_resident ~ normal( 0 , 10 );
    female_age1_resident ~ normal( 0 , 10 );
    female_age2_resident ~ normal( 0 , 10 );    for ( i in 1:N ) {        p_part[i] = a + fishing*COMM_TYPE_F[i] + 
			male_age1_inmigrant*MALE[i]*AGE1[i]*INMIGRANT[i] + 
			male_age2_inmigrant*MALE[i]*AGE2[i]*INMIGRANT[i] + 
			male_age3_inmigrant*MALE[i]*INMIGRANT[i] + 
			female_age1_inmigrant*(1-MALE[i])*AGE1[i]*INMIGRANT[i] + 
			female_age2_inmigrant*(1-MALE[i])*AGE2[i]*INMIGRANT[i] + 
			female_age3_inmigrant*(1-MALE[i])*INMIGRANT[i] + 
			male_age1_resident*MALE[i]*AGE1[i]*(1-INMIGRANT[i]) + 
			male_age2_resident*MALE[i]*AGE2[i]*(1-INMIGRANT[i]) + 
			male_age3_resident*MALE[i]*(1-INMIGRANT[i]) + 
			female_age1_resident*(1-MALE[i])*AGE1[i]*(1-INMIGRANT[i]) + 
			female_age2_resident*(1-MALE[i])*AGE2[i]*(1-INMIGRANT[i]);        p_part[i] = inv_logit(p_part[i]);    }    SUC ~ beta_binomial( TRIAL , p_part*dispersion , (1-p_part)*dispersion );}
