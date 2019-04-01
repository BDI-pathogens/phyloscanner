# input df.sampling: individual covariates, the number of trials and success
#       dc: transmission pairs, covariates for transmitters' group and recipients' group
#       iteration: the number of iterations
#       sample.method: 'empirical' and 'betaapprox'
#       glm.model: choice of binomial or beta-binomial model for inferring sampling rate
#                  it has the same length as the number of data table in df.sampling
#       seed: random seeds
samples.from.GLM.prior<-function(df.sampling,dc,iteration,sample.method,glm.model,seed=NULL){
  library(data.table)
  #library(rstan)
  #library(fitdistrplus)
  set.seed(seed)

  # take unique subpopulation that appears in transmitters and recipients in dc
  tr<-dc[, grep("TR_", names(dc)), with = FALSE]
  setnames(tr, colnames(tr), gsub('TR_','',colnames(tr)))
  tr<-unique(tr)
  rec<-dc[, grep("REC_", names(dc)), with = FALSE]
  setnames(rec, colnames(rec), gsub('REC_','',colnames(rec)))
  rec<-unique(rec)
  unique.group<-unique(rbind(tr,rec))
  # setkey(unique.group,COMM_NUM_A,SEX,AGE_AT_MID_C,INMIGRANT)
  setkey(unique.group,CATEGORY)
  unique.group<-merge(subset(unique.group,select = 'CATEGORY'),dg,by='CATEGORY')
  unique.group[,CATEGORY_ID:=1:nrow(unique.group)]

  # take unique pairs of subpopulations for transmission flows
  unique.pairs<-subset(dc,select=c('COUNT_ID','TR_CATEGORY','REC_CATEGORY'))
  tmp<-subset(unique.group,select = c('CATEGORY','CATEGORY_ID'))
  setnames(tmp,colnames(tmp),paste0('TR_',colnames(tmp)))
  unique.pairs<-merge(unique.pairs,tmp,by='TR_CATEGORY')
  setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
  unique.pairs<-merge(unique.pairs,tmp,by='REC_CATEGORY')
  setkey(unique.pairs,COUNT_ID)

  # number of prior samples
  nprior<-ceiling(iteration/(nrow(unique.group)+1))+1

  S.DTL.prior<-list()
  fit.list<-list() #model
  S.DTL.prior$SAM_P<-matrix(1,nrow=nprior,ncol=nrow(unique.group)) # sampling rate
  S.DTL.prior$SAM_P_LOGD<-matrix(0,nrow=nprior,ncol=nrow(unique.group)) # log density
  S.DTL.prior$EMP_S<-1

  for (i in 1L:length(df.sampling)){
    # record priors of sampling rate for the ith source of bias
    SAM_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(unique.group))
    SAM_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(unique.group))
    EMP_S<-(sum(df.sampling[[i]]$SUC)/sum(df.sampling[[i]]$TRIAL))^2

    # predict sampling rates for stratum a using logistic regression
    ######################## will depend on covariates??? not general???
    if (glm.model[i]=='binomial'){
      data.glm<-as.list(subset(df.sampling[[i]],select=c('COMM_NUM_B','TRIAL','AGE1','AGE2','MALE','SUC')))
      data.glm$N<-nrow(df.sampling[[i]])
      fit <- stan(file = 'glm_age3sex2comm2bin.stan', data = data.glm, iter=max(1e4,5e2+nprior),warmup = 5e2,
                  cores = 4,chains = 1,init = list(list(a=0, comm=rep(0,2), sig_comm=1, male=0, age1=0, age2=0)))
    }else if(glm.model[i]=='beta_binomial'){
      data.glm<-as.list(subset(df.sampling[[i]],select=c('COMM_NUM_B','TRIAL','AGE1','AGE2','MALE','SUC')))
      data.glm$N<-nrow(df.sampling[[i]])
      fit <- stan(file = 'glm_age3sex2comm2betabin.stan', data = data.glm, iter=max(1e4,5e2+nprior),warmup = 5e2,
                  cores = 4,chains = 1,init = list(a=0, comm=rep(0,2), sig_comm=1, male=0, age1=0, age2=0,dispersion=1))
    }
    # extract samples for the parameters
    tmp		<- extract(fit,permute=TRUE)
    if(sample.method=='empirical'){
      for (j in 1:nrow(unique.group)){
        # posterior samples of p_suc
        # tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]] + tmp$trading*unique.group$COMM_TYPE_T[j] + tmp$fishing*unique.group$COMM_TYPE_F[j] + tmp$male*unique.group$MALE[j] +
        #   tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j] + tmp$age3*unique.group$AGE3[j] +  tmp$age4*unique.group$AGE4[j] +
        #   tmp$age5*unique.group$AGE5[j]  + tmp$age6*unique.group$AGE6[j] + tmp$inmigrant*unique.group$INMIGRANT[j]
        tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]]  + tmp$male*unique.group$MALE[j] +
          tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j]
        tmpp<-exp(tmp2)/(1+exp(tmp2))
        # draw samples from posterior samples of p_suc
        SAM_P[,j]<-sample(tmpp,nprior,replace=TRUE)
        # calculate the log empirical density estimate
        d.sam<-approxfun(density(tmpp))
        SAM_P_LOGD[,j]<-log(d.sam(SAM_P[,j]))
      }
    }else if(sample.method=='betaapprox'){
      for (j in 1:nrow(unique.group)){
        # posterior samples of p_suc
        # tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]] + tmp$trading*unique.group$COMM_TYPE_T[j] + tmp$fishing*unique.group$COMM_TYPE_F[j] + tmp$male*unique.group$MALE[j] +
        #   tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j] + tmp$age3*unique.group$AGE3[j] +  tmp$age4*unique.group$AGE4[j] +
        #   tmp$age5*unique.group$AGE5[j]  + tmp$age6*unique.group$AGE6[j] + tmp$inmigrant*unique.group$INMIGRANT[j]
        tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]]  + tmp$male*unique.group$MALE[j] +
          tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j]
        tmpp<-exp(tmp2)/(1+exp(tmp2))
        # fit a beta distribution for p_suc
        betapar<-fitdist(as.vector(tmpp),"beta")
        # draw samples from the fitted beta distribution
        SAM_P[,j]<-rbeta(nprior,betapar$estimate[1],betapar$estimate[2])
        # calculate the log density from the fitted beta distribution
        SAM_P_LOGD[,j]<-dbeta(S.DTL.prior$SEQ_P[,j],shape1 = betapar$estimate[1], shape2 = betapar$estimate[2],log=TRUE)
      }
    }
    S.DTL.prior$SAM_P<-S.DTL.prior$SAM_P*SAM_P
    S.DTL.prior$SAM_P_LOGD<-S.DTL.prior$SAM_P+SAM_P_LOGD
    S.DTL.prior$EMP_S<-S.DTL.prior$EMP_S*EMP_S
    fit.list[[i]]<-fit
  }
  S.DTL.prior$method<-'GLM'
  TR.OBS<-data.table(OBS=dc$OBS,TR_CATEGORY_ID=unique.pairs$TR_CATEGORY_ID,
                     REC_CATEGORY_ID=unique.pairs$REC_CATEGORY_ID)
  save(S.DTL.prior,fit.list,TR.OBS,file=paste0('InputsMCMC_GLM',iteration,'.rda'))
}



samples.from.SARWS.prior<-function(df.sampling,dc,iteration,seed=NULL){
  library(data.table)

  # take unique subpopulation that appears in transmitters and recipients in dc
  tr<-dc[, grep("TR_", names(dc)), with = FALSE]
  setnames(tr, colnames(tr), gsub('TR_','',colnames(tr)))
  tr<-unique(tr)
  rec<-dc[, grep("REC_", names(dc)), with = FALSE]
  setnames(rec, colnames(rec), gsub('REC_','',colnames(rec)))
  rec<-unique(rec)
  unique.group<-unique(rbind(tr,rec))
  # setkey(unique.group,COMM_NUM_A,SEX,AGE_AT_MID_C)
  setkey(unique.group,CATEGORY)
  unique.group[,CATEGORY_ID:=1:nrow(unique.group)]

  # take unique pairs of subpopulations for transmission flows
  unique.pairs<-subset(dc,select=c('COUNT_ID','TR_CATEGORY','REC_CATEGORY'))
  tmp<-subset(unique.group,select = c('CATEGORY','CATEGORY_ID'))
  setnames(tmp,colnames(tmp),paste0('TR_',colnames(tmp)))
  unique.pairs<-merge(unique.pairs,tmp,by='TR_CATEGORY')
  setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
  unique.pairs<-merge(unique.pairs,tmp,by='REC_CATEGORY')
  setkey(unique.pairs,COUNT_ID)

  # number of prior samples
  nprior<-ceiling(iteration/(nrow(unique.group)+1))+1

  S.DTL.prior<-list()
  fit.list<-list() #model
  S.DTL.prior$SAM_P<-matrix(1,nrow=nprior,ncol=nrow(unique.group)) # sampling rate
  S.DTL.prior$SAM_P_LOGD<-matrix(0,nrow=nprior,ncol=nrow(unique.group)) # log density
  S.DTL.prior$EMP_S<-1

  for (i in 1L:length(df.sampling)){
    # record priors of sampling rate for the ith source of bias
    SAM_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(unique.group))
    SAM_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(unique.group))
    EMP_S<-NA_real_
    #unique.group.tmp<-merge(unique.group,df.sampling[[i]],by=c('COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT'))
    unique.group.tmp<-merge(unique.group,df.sampling[[i]],by=c('CATEGORY'))

    unique.group.tmp[,ALPHA:=SUC+1]
    unique.group.tmp[,BETA:=TRIAL-SUC+1]
    EMP_S<-(sum(df.sampling[[i]]$SUC)/sum(df.sampling[[i]]$TRIAL))^2

    for (j in 1:nrow(unique.group)){
      SAM_P[,j]<-rbeta(nprior,unique.group.tmp$ALPHA[j],unique.group.tmp$BETA[j])
      SAM_P_LOGD[,j]<-dbeta(SAM_P[,j],unique.group.tmp$ALPHA[j],unique.group.tmp$BETA[j],log=TRUE)
    }

    S.DTL.prior$SAM_P<-S.DTL.prior$SAM_P*SAM_P
    S.DTL.prior$SAM_P_LOGD<-S.DTL.prior$SAM_P+SAM_P_LOGD
    S.DTL.prior$EMP_S<-S.DTL.prior$EMP_S*EMP_S
  }
  S.DTL.prior$method<-'SARWS'
  # extract corresponding groups for pairs
  TR.OBS<-data.table(OBS=dc$OBS,TR_CATEGORY_ID=unique.pairs$TR_CATEGORY_ID,
                     REC_CATEGORY_ID=unique.pairs$REC_CATEGORY_ID)
  save(S.DTL.prior,TR.OBS,file=paste0('InputsMCMC_SARWS',iteration,'.rda'))
}
