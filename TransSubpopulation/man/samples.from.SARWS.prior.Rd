\name{samples.from.SARWS.prior}
\alias{samples.from.SARWS.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drawing samples from the prior distribution under sampling at random within strata (SARWS) assumption.
}
\description{
This function provides the inputs for mcmc sampler. Samples are drawn from the prior distribution of the sampling rate and evaluate corresponding density. The prior distribution is informative and given by a posterior distribution inferred from Bayesian hierarchical model using the number of individuals who are eligible and sampled under the assumption of missing at random within strata (SARWS).
}
\usage{
samples.from.SARWS.prior(df.sampling,dc,iteration,seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df.sampling}{
a list of data tables each of which includes data used for infering the distribution of sampling rates through a hierarchical model, i.e. the number of eligible individuals and individual sampled in each subpopulation.
}
\item{iteration}{
the number of iterations for mcmc sampler, which decides the number of prior samples needed.
}
\item{dc}{
 a data table includes transmission counts between subpopulations and covariates of transmitters' and recipients' subpopulations.
}
\item{seed}{
the random seed.
}
}

\examples{
library(data.table)
set.seed(123)
# dg
dg<-data.table(CATEGORY=c(1,2),TRIAL=c(2000,2500),P_SEQ_EMP=c(0.6,0.45))
dg[,SUC:=TRIAL*P_SEQ_EMP]
# dc
dc<-data.table(TR_CATEGORY=c(1,2,1,2),REC_CATEGORY=c(1,1,2,2))
tmp<-subset(dg,select=c('CATEGORY','P_SEQ_EMP'))
setnames(tmp,colnames(tmp),paste0('TR_',colnames(tmp)))
dc<-merge(dc,tmp,by='TR_CATEGORY')
setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
dc<-merge(dc,tmp,by='REC_CATEGORY')
dc[,S:=TR_P_SEQ_EMP*REC_P_SEQ_EMP]

# transmission proportions
TRUE_PI<-c(0.36,0.04,0.06,0.54)

# assume observed transmission counts are 300
N<-rpois(1,300/mean(dc$S))
z<-rmultinom(1,size=N,prob=TRUE_PI)
n<-matrix(NA_integer_,ncol=1,nrow=length(TRUE_PI))
for (i in 1:length(TRUE_PI)){
  n[i]<-rbinom(1,size=z[i],dc$S[i])
}

dc$OBS<-n
dc$COUNT_ID<-1:nrow(dc)

# remove irrelevant terms
dg$P_SEQ_EMP<-NULL
dc$TR_P_SEQ_EMP<-NULL
dc$REC_P_SEQ_EMP<-NULL
dc$S<-NULL

df.sampling<-list()
df.sampling[[1]]<-dg

samples.from.SARWS.prior(df.sampling=df.sampling,dc=dc,iteration=1e6,seed=123)
}
