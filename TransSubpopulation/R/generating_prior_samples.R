# input dp: individual covariates, PART, ELIG
#       dg: individual covariates, HIV, SEQ
#       dc: transmission pairs, covariates for transmitters' group and recipients' group
#       method: 'empirical' and 'betaapprox'
#       biastype: 'sequence' and 'sampling'
samples.from.GLM.prior<-function(dp=NULL,dg,dc,iteration,method,biastype){
  # number of prior samples
  nprior<-ceiling(iteration/(nrow(dc)+1))+1

  require(data.table)
  require(rethinking)
  require(fitdistrplus)

  S.DTL.prior<-list()
  S.DTL.prior$TR_SEQ_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
  S.DTL.prior$TR_SEQ_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
  S.DTL.prior$REC_SEQ_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
  S.DTL.prior$REC_SEQ_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))

  # predict sampling rates for stratum a using logistic regression
  ms 	<- map2stan(
    alist(
      SEQ ~ dbinom(HIV, p_seq),
      logit(p_seq) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T +
        fishing*COMM_TYPE_F + male*MALE + young*AGE_YOUNG + midage*AGE_MID,
      a ~ dnorm(0, 100),
      comm ~ dnorm(0, sig_comm),
      sig_comm ~ dcauchy(0,1),
      c(trading, fishing, male, young, midage) ~ dnorm(0,10)
    ),
    data=as.data.frame(dg),
    start=list(a=0, comm=rep(0,36), sig_comm=1, trading=0, fishing=0, male=0, young=0, midage=0),
    warmup=5e2, iter=1e4, chains=1, cores=4
  )		 # iteration should be larger than the number of samples needed

  # # posterior predictive check
  # sims 	<- sim(ms, n=1000)
  # tmp		<- apply(sims, 2, median)
  # dpp 	<- apply(sims, 2, PI, prob=0.95)
  # dpp		<- rbind(tmp, dpp)
  # rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
  # dpp		<- as.data.table(t(dpp))
  # dpp		<- cbind(dg, dpp)
  # dpp[, c(mean(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95), sum(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95))]

  # extract posterior samples
  tmp		<- extract.samples(ms)

  if(method=='empirical'){
    for (i in 1:nrow(dc)){
      tmp2<-tmp$a + tmp$comm[,dc$TR_COMM_NUM_B[i]] + tmp$trading*dc$TR_COMM_TYPE_T[i] + tmp$fishing*dc$TR_COMM_TYPE_F[i] + tmp$male*dc$TR_MALE[i] +
        tmp$young*dc$TR_AGE_YOUNG[i] + tmp$midage*dc$TR_AGE_MID[i]
      tmpp<-exp(tmp2)/(1+exp(tmp2))
      S.DTL.prior$TR_SEQ_P[,i]<-sample(tmpp,nprior,replace=FALSE)
      d.tr.seq<-approxfun(density(tmpp))
      S.DTL.prior$TR_SEQ_P_LOGD[,i]<-d.tr.seq(S.DTL.prior$TR_SEQ_P[,i])

      tmp2<-tmp$a + tmp$comm[,dc$REC_COMM_NUM_B[i]] + tmp$trading*dc$REC_COMM_TYPE_T[i] + tmp$fishing*dc$REC_COMM_TYPE_F[i] + tmp$male*dc$REC_MALE[i] +
        tmp$young*dc$REC_AGE_YOUNG[i] + tmp$midage*dc$REC_AGE_MID[i]
      tmpp<-exp(tmp2)/(1+exp(tmp2))
      S.DTL.prior$REC_SEQ_P[,i]<-sample(tmpp,nprior,replace=FALSE)
      d.rec.seq<-approxfun(density(tmpp))
      S.DTL.prior$REC_SEQ_P_LOGD[,i]<-d.rec.seq(S.DTL.prior$REC_SEQ_P[,i])
    }
  }else if(method=='betaapprox'){
    for (i in 1:nrow(dc)){
      tmp2<-tmp$a + tmp$comm[,dc$TR_COMM_NUM_B[i]] + tmp$trading*dc$TR_COMM_TYPE_T[i] + tmp$fishing*dc$TR_COMM_TYPE_F[i] + tmp$male*dc$TR_MALE[i] +
        tmp$young*dc$TR_AGE_YOUNG[i] + tmp$midage*dc$TR_AGE_MID[i]
      tmpp<-exp(tmp2)/(1+exp(tmp2))
      betapar<-fitdist(as.vector(tmpp),"beta")
      S.DTL.prior$TR_SEQ_P[,i]<-rbeta(nprior,betapar$estimate[1],betapar$estimate[2])
      S.DTL.prior$TR_SEQ_P_LOGD[,i]<-dbeta(S.DTL.prior$TR_SEQ_P[,i],shape1 = betapar$estimate[1], shape2 = betapar$estimate[2],log=TRUE)

      tmp2<-tmp$a + tmp$comm[,dc$REC_COMM_NUM_B[i]] + tmp$trading*dc$REC_COMM_TYPE_T[i] + tmp$fishing*dc$REC_COMM_TYPE_F[i] + tmp$male*dc$REC_MALE[i] +
        tmp$young*dc$REC_AGE_YOUNG[i] + tmp$midage*dc$REC_AGE_MID[i]
      tmpp<-exp(tmp2)/(1+exp(tmp2))
      betapar<-fitdist(as.vector(tmpp),"beta")
      S.DTL.prior$REC_SEQ_P[,i]<-rbeta(nprior,betapar$estimate[1],betapar$estimate[2])
      S.DTL.prior$REC_SEQ_P_LOGD[,i]<-dbeta(S.DTL.prior$REC_SEQ_P[,i],shape1 = betapar$estimate[1], shape2 = betapar$estimate[2],log=TRUE)
    }
  }

  if (biastype=='sampling'){
    S.DTL.prior$TR_PART_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
    S.DTL.prior$TR_PART_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
    S.DTL.prior$REC_PART_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
    S.DTL.prior$REC_PART_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))

    ms <- map2stan(
      alist(
        PART ~ dbetabinom(ELIG, p_part, dispersion),
        logit(p_part) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F +
          male*MALE + young*AGE_YOUNG + midage*AGE_MID,
        a ~ dnorm(0, 100),
        comm ~ dnorm(0, sig_comm),
        sig_comm ~ dcauchy(0,1),
        dispersion ~ dexp(1),
        c(trading, fishing, male, young, midage) ~ dnorm(0,10)
      ),
      data=as.data.frame(dp),
      start=list(	a=0, comm=rep(0,36), dispersion=10, sig_comm=1, trading=0, fishing=0, male=0, young=0, midage=0),
      warmup=5e2, iter=1e4, chains=1, cores=4
    )		# note ms2 and ms4 don't work well

    # # posterior predictive check
    # sims 	<- sim(ms, n=1000)
    # tmp		<- apply(sims, 2, median)
    # dpp 	<- apply(sims, 2, PI, prob=0.95)
    # dpp		<- rbind(tmp, dpp)
    # rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
    # dpp		<- as.data.table(t(dpp))
    # dpp		<- cbind(dp, dpp)
    # dpp[, c(mean(PART<predicted_obs_l95 | PART>predicted_obs_u95), sum(PART<predicted_obs_l95 | PART>predicted_obs_u95))]
    # #	  0.01388889 6.00000000

    tmp		<- extract.samples(ms)

    if(method=='empirical'){
      for (i in 1:nrow(dc)){
        tmp2<-tmp$a + tmp$comm[,dc$TR_COMM_NUM_B[i]] + tmp$trading*dc$TR_COMM_TYPE_T[i] + tmp$fishing*dc$TR_COMM_TYPE_F[i] + tmp$male*dc$TR_MALE[i] +
          tmp$young*dc$TR_AGE_YOUNG[i] + tmp$midage*dc$TR_AGE_MID[i]
        tmpp<-exp(tmp2)/(1+exp(tmp2))
        S.DTL.prior$TR_PART_P[,i]<-sample(tmpp,nprior,replace=FALSE)
        d.tr.part<-approxfun(density(tmpp))
        S.DTL.prior$TR_PART_P_LOGD[,i]<-d.tr.part(S.DTL.prior$TR_PART_P[,i])

        tmp2<-tmp$a + tmp$comm[,dc$REC_COMM_NUM_B[i]] + tmp$trading*dc$REC_COMM_TYPE_T[i] + tmp$fishing*dc$REC_COMM_TYPE_F[i] + tmp$male*dc$REC_MALE[i] +
          tmp$young*dc$REC_AGE_YOUNG[i] + tmp$midage*dc$REC_AGE_MID[i]
        tmpp<-exp(tmp2)/(1+exp(tmp2))
        S.DTL.prior$REC_PART_P[,i]<-sample(tmpp,nprior,replace=FALSE)
        d.rec.part<-approxfun(density(tmpp))
        S.DTL.prior$REC_PART_P_LOGD[,i]<-d.rec.part(S.DTL.prior$REC_PART_P[,i])
      }
    }else if(method=='betaapprox'){
      for (i in 1:nrow(dc)){
        tmp2<-tmp$a + tmp$comm[,dc$TR_COMM_NUM_B[i]] + tmp$trading*dc$TR_COMM_TYPE_T[i] + tmp$fishing*dc$TR_COMM_TYPE_F[i] + tmp$male*dc$TR_MALE[i] +
          tmp$young*dc$TR_AGE_YOUNG[i] + tmp$midage*dc$TR_AGE_MID[i]
        tmpp<-exp(tmp2)/(1+exp(tmp2))
        betapar<-fitdist(as.vector(tmpp),"beta")
        S.DTL.prior$TR_PART_P[,i]<-rbeta(nprior,betapar$estimate[1],betapar$estimate[2])
        S.DTL.prior$TR_PART_P_LOGD[,i]<-dbeta(S.DTL.prior$TR_PART_P[,i],shape1 = betapar$estimate[1], shape2 = betapar$estimate[2],log=TRUE)

        tmp2<-tmp$a + tmp$comm[,dc$REC_COMM_NUM_B[i]] + tmp$trading*dc$REC_COMM_TYPE_T[i] + tmp$fishing*dc$REC_COMM_TYPE_F[i] + tmp$male*dc$REC_MALE[i] +
          tmp$young*dc$REC_AGE_YOUNG[i] + tmp$midage*dc$REC_AGE_MID[i]
        tmpp<-exp(tmp2)/(1+exp(tmp2))
        betapar<-fitdist(as.vector(tmpp),"beta")
        S.DTL.prior$REC_PART_P[,i]<-rbeta(nprior,betapar$estimate[1],betapar$estimate[2])
        S.DTL.prior$REC_PART_P_LOGD[,i]<-dbeta(S.DTL.prior$REC_PART_P[,i],shape1 = betapar$estimate[1], shape2 = betapar$estimate[2],log=TRUE)
      }
    }
  }
  S.DTL.prior$EMP_S<-dc$S
  return(S.DTL.prior)
}

samples.from.SARWS.prior<-function(dc,iteration,biastype){
  nprior<-ceiling(iteration/(nrow(dc)+1))+1

  S.DTL.prior<-list()
  S.DTL.prior$TR_SEQ_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
  S.DTL.prior$TR_SEQ_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
  S.DTL.prior$REC_SEQ_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
  S.DTL.prior$REC_SEQ_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))

  for (i in 1:nrow(dc)){
    S.DTL.prior$TR_SEQ_P[,i]<-rbeta(nprior,dc$TR_P_SEQ_ALPHA[i],dc$TR_P_SEQ_BETA[i])
    S.DTL.prior$TR_SEQ_P_LOGD[,i]<-dbeta(S.DTL.prior$TR_SEQ_P[,i],dc$TR_P_SEQ_ALPHA[i],dc$TR_P_SEQ_BETA[i],log=TRUE)
    S.DTL.prior$REC_SEQ_P[,i]<-rbeta(nprior,dc$REC_P_SEQ_ALPHA[i],dc$REC_P_SEQ_BETA[i])
    S.DTL.prior$REC_SEQ_P_LOGD[,i]<-dbeta(S.DTL.prior$REC_SEQ_P[,i],dc$REC_P_SEQ_ALPHA[i],dc$REC_P_SEQ_BETA[i],log=TRUE)
  }

  if (biastype=='sampling'){
    S.DTL.prior$TR_PART_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
    S.DTL.prior$TR_PART_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
    S.DTL.prior$REC_PART_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))
    S.DTL.prior$REC_PART_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(dc))

    for (i in 1:nrow(dc)){
      S.DTL.prior$TR_PART_P[,i]<-rbeta(nprior,dc$TR_P_PART_ALPHA[i],dc$TR_P_PART_BETA[i])
      S.DTL.prior$TR_PART_P_LOGD[,i]<-dbeta(S.DTL.prior$TR_PART_P[,i],dc$TR_P_PART_ALPHA[i],dc$TR_P_PART_BETA[i],log=TRUE)
      S.DTL.prior$REC_PART_P[,i]<-rbeta(nprior,dc$REC_P_PART_ALPHA[i],dc$REC_P_PART_BETA[i])
      S.DTL.prior$REC_PART_P_LOGD[,i]<-dbeta(S.DTL.prior$REC_PART_P[,i],dc$REC_P_PART_ALPHA[i],dc$REC_P_PART_BETA[i],log=TRUE)
    }
  }
  S.DTL.prior$EMP_S<-dc$S

  return(S.DTL.prior)
}
