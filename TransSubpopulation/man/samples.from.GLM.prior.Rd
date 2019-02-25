\name{samples.from.GLM.prior}
\alias{samples.from.GLM.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drawing samples from the prior distribution given by generalised linear model.
}
\description{
This function provides the inputs for mcmc sampler. Samples are drawn from the prior distribution of the sampling rate and the corresponding densities are evaluated. The prior distribution is informative and given by a posterior distribution inferred from Bayesian hierarchical model using empirical sampling rate and group covariates.
}
\usage{
samples.from.GLM.prior(df.sampling,dc,iteration,sample.method,glm.model,seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df.sampling}{
a list of data tables each of which includes data used for infering the distribution of sampling rates through a hierarchical model, i.e. the number of eligible individuals and individual sampled in each subpopulation, the covariates of individuals in a subpopulation.
}
\item{dc}{
 a data table includes transmission counts between subpopulations and covariates of transmitters' and recipients' subpopulations.
}
\item{iteration}{
the number of iterations for mcmc sampler, which decides the number of prior samples needed.
}
\item{sample.method}{
the methods for evaluating density of each samples. method='empirical' indicates density of samples come from an empirical distribution, method='betaapprox' indicates a beta distribution is used to approximate the distribution of samples.
}
\item{glm.model}{
the type of glm model for infering the sampling rates. glm.model='binomial' indicates using binomial model, glm.model='beta-binomial' indicates using beta-binomial model. Choices of glm.model needed to be determined by posterior predictive checks.
}
\item{seed}{
the random seed.
}
}

\examples{
library(data.table)
library(Boom)
set.seed(123)
# assume 120 infected individuals in 1 and 80 infected individuals in 2
# 70\% are female in 1 and 60\% are females in 2
# in 1 F: 25\% 15-24 60\% 25-34 15\% 35+
# in 1 M: 20\% 15-24 50\% 25-34 30\% 35+
# in 2 F: 25\% 15-24 50\% 25-34 25\% 35+
# in 2 M: 20\% 15-24 40\% 25-34 40\% 35+
dg<-data.table(COMM_NUM_B=c(rep(1,120),rep(2,80)),INF=rep(1,200))
num.f.1<-round(120*0.7)
num.m.1<-120-round(120*0.7)
num.f.2<-round(80*0.6)
num.m.2<-80-round(80*0.6)

dg$SEX<-c(rep('F',num.f.1),
          rep('M',num.m.1),
          rep('F',num.f.2),
          rep('M',num.m.2))
dg$AGE_AT_MID_C<-c(rep('15-24',round(num.f.1*0.25)),
                   rep('25-34',round(num.f.1*0.6)),
                   rep('35+',num.f.1-round(num.f.1*0.25)-round(num.f.1*0.6)),
                   rep('15-24',round(num.m.1*0.2)),
                   rep('25-34',round(num.m.1*0.5)),
                   rep('35+',num.m.1-round(num.m.1*0.2)-round(num.m.1*0.5)),
                   rep('15-24',round(num.f.2*0.25)),
                   rep('25-34',round(num.f.2*0.5)),
                   rep('35+',num.f.2-round(num.f.2*0.25)-round(num.f.2*0.5)),
                   rep('15-24',round(num.m.2*0.2)),
                   rep('25-34',round(num.m.2*0.4)),
                   rep('35+',num.m.2-round(num.m.2*0.2)-round(num.m.2*0.4))
                   )

# regression coefficients
dg		<- dg[, list(TRIAL=length(INF)),by=c('AGE_AT_MID_C','SEX','COMM_NUM_B')]
comm<-c(0.5,-0.5)
a<-0.4
male<-0.4
young<-0.3
midage<-0.2

#	binarize age, sex
dg[, AGE1:= as.integer(AGE_AT_MID_C=='15-24')]
dg[, AGE2:= as.integer(AGE_AT_MID_C=='25-34')]
dg[, MALE:= as.integer(SEX=='M')]

dg[,LOGIT_P_SEQ:=a + comm[COMM_NUM_B] + male*MALE + young*AGE1 + midage*AGE2]
dg[,P_SEQ:=exp(LOGIT_P_SEQ)/(1+exp(LOGIT_P_SEQ))]
dg[,SUC:=rbinom(nrow(dg),TRIAL,P_SEQ)]

dg[,CATEGORY:=paste0(COMM_NUM_B,':',SEX,':',AGE_AT_MID_C)]
dg[,COMM_NUM_B:=as.integer(COMM_NUM_B)]

#dc
dc<-as.data.table(expand.grid(TR_CATEGORY=dg$CATEGORY,REC_CATEGORY=dg$CATEGORY))
dc[,COUNT_ID:=1:nrow(dc)]
dc[, c("TR_COMM_NUM_B","TR_SEX", "TR_AGE_AT_MID_C") := tstrsplit(TR_CATEGORY, ":", fixed=TRUE)]
dc[, c("REC_COMM_NUM_B","REC_SEX", "REC_AGE_AT_MID_C") := tstrsplit(REC_CATEGORY, ":", fixed=TRUE)]

dc<-dc[TR_SEX!=REC_SEX,]

# transmission proportions
lbd<-matrix(NA_real_,nrow = 1, ncol = nrow(dc))
for (i in 1:nrow(dc)){
  if (dc$TR_COMM_NUM_B[i]==dc$REC_COMM_NUM_B[i]){
    if ( (dc$TR_AGE_AT_MID_C[i]=='15-24' & dc$REC_AGE_AT_MID_C[i]=='25-34') |
         (dc$TR_AGE_AT_MID_C[i]=='15-24' & dc$REC_AGE_AT_MID_C[i]=='35+') |
         (dc$TR_AGE_AT_MID_C[i]=='25-34' & dc$REC_AGE_AT_MID_C[i]=='35+')){
      if (dc$TR_SEX[i]=='M'){
        lbd[i]<-rnorm(1,1,0.2)
      }else{
        lbd[i]<-rnorm(1,2,0.2)
      }
    }else if ( (dc$REC_AGE_AT_MID_C[i]=='15-24' & dc$TR_AGE_AT_MID_C[i]=='25-34') |
               (dc$REC_AGE_AT_MID_C[i]=='15-24' & dc$TR_AGE_AT_MID_C[i]=='35+') |
               (dc$REC_AGE_AT_MID_C[i]=='25-34' & dc$TR_AGE_AT_MID_C[i]=='35+')){
      if (dc$TR_SEX[i]=='M'){
        lbd[i]<-rnorm(1,2,0.2)
      }else{
        lbd[i]<-rnorm(1,1,0.2)
      }
    }else if(dc$TR_AGE_AT_MID_C[i]==dc$REC_AGE_AT_MID_C[i]){
      lbd[i]<-rnorm(1,2,0.2)
    }
  }else{
    lbd[i]<-rnorm(1,1,0.2)
  }
}
TRUE_PI<-rdirichlet(1,lbd)

tmp<-subset(dg,select=c('CATEGORY','P_SEQ'))
setnames(tmp,colnames(tmp),paste0('TR_',colnames(tmp)))
dc<-merge(dc,tmp,by='TR_CATEGORY')
setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
dc<-merge(dc,tmp,by='REC_CATEGORY')
dc[,S:=TR_P_SEQ*REC_P_SEQ]

N<-rpois(1,20/mean(dc$S))
z<-rmultinom(1,size=N,prob=TRUE_PI)
n<-matrix(NA_integer_,ncol=1,nrow=length(TRUE_PI))
for (i in 1:length(TRUE_PI)){
  n[i]<-rbinom(1,size=z[i],dc$S[i])
}

dc$OBS<-n
dc.obs.zero<-which(dc$OBS==0)
dc<-dc[which(dc$OBS!=0),]
dc$COUNT_ID<-1:nrow(dc)

# remove irrelevant terms
dg$P_SEQ<-NULL
dg$LOGIT_P_SEQ<-NULL
dc$TR_P_SEQ<-NULL
dc$REC_P_SEQ<-NULL
dc$S<-NULL

df.sampling<-list()
df.sampling[[1]]<-dg

samples.from.GLM.prior(df.sampling=df.sampling,dc=dc,iteration=1e6,sample.method='empirical',glm.model='binomial',seed=123)
}
