library(data.table)
library(Boom)
# load simulated counts from ode
load(paste0('data.rda'))

# set id for transmission
dobs[,TRM_CAT_PAIR_ID:=seq_len(nrow(dobs))]

# adjusted flow
mcmc.file	<- 'mcmc_sarws'

nprior<-1000
dprior<-dss[,list(P=rbeta(nprior,SUC+0.5,TRIAL-SUC+0.5),
                  SAMPLE=1:nprior),
            by='SAMPLING_CATEGORY']
dprior<-merge(dprior,
              subset(dss,select=c('SAMPLING_CATEGORY','TRIAL','SUC')),
              by='SAMPLING_CATEGORY')
dprior[,LP:=dbeta(P,SUC+0.5,TRIAL-SUC+0.5,log=TRUE)]
set(dprior,NULL,c('TRIAL','SUC'),NULL)

control		<- list(seed=42, 
                 mcmc.n=4 * 1e5,							
                 verbose=0, 
                 outfile=mcmc.file,
                 sweep_group=1e5L)
source.attribution.mcmc(dobs, dprior, control=control)

# unadjusted flow
mcmc.file	<- 'mcmc_sar'

dss <- data.table(SAMPLING_CATEGORY=c(1,2),TRIAL=c(50,50),P_SEQ_EMP=c(0.5,0.5))
dss[,SUC:=TRIAL*P_SEQ_EMP]
dss[,P_SEQ_EMP:=NULL]

nprior<-1000
dprior<-dss[,list(P=rbeta(nprior,SUC+0.5,TRIAL-SUC+0.5),
                  SAMPLE=1:nprior),
            by='SAMPLING_CATEGORY']
dprior<-merge(dprior,
              subset(dss,select=c('SAMPLING_CATEGORY','TRIAL','SUC')),
              by='SAMPLING_CATEGORY')
dprior[,LP:=dbeta(P,SUC+0.5,TRIAL-SUC+0.5,log=TRUE)]
set(dprior,NULL,c('TRIAL','SUC'),NULL)


control=list(seed=42,
             mcmc.n=4*1e5, 
             verbose=0,
             outfile=mcmc.file,
             sweep_group=1e5L)

source.attribution.mcmc(dobs, dprior, control)

