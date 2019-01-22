require(data.table)

# define sequencing rate
df<-data.table(CATEGORY=c(1,2),P_SEQ_EMP=c(0.6,0.5),HIV_YES=c(2000,2500))
df$DEEP_SEQ<-df[,HIV_YES*P_SEQ_EMP]
df$P_SEQ_ALPHA<-df$DEEP_SEQ+1
df$P_SEQ_BETA<-df$HIV_YES-df$DEEP_SEQ+1

# define transmission pairs
dc<-data.table(REC_CATEGORY=c(1,1,2,2),TR_CATEGORY=c(1,2,1,2))

setnames(df, colnames(df), paste0('TR_',colnames(df)))
dc	<- merge(dc, df, by=c('TR_CATEGORY'))
setnames(df, colnames(df), gsub('TR_','REC_',colnames(df)))
dc	<- merge(dc, df, by=c('REC_CATEGORY'))

dc[,S:=TR_P_SEQ_EMP*REC_P_SEQ_EMP]

# define the proportion of transmissions for each pair
TRUE_PI<-c(0.36,0.04,0.06,0.54)

# simulate the number of transmissions for each pair
N<-rpois(1,300/mean(dc$S))
z<-rmultinom(1,size=N,prob=TRUE_PI)
n<-matrix(NA_integer_,ncol=1,nrow=length(TRUE_PI))
for (i in 1:length(TRUE_PI)){
  n[i]<-rbinom(1,size=z[i],dc$S[i])
}

# make the data frame dc
dc$TR_OBS<-n
dc$TRUE_OBS<-z
dc$TRUE_PI<-TRUE_PI
dc$COUNT_ID<-1:nrow(dc)


# example
mcmc.sarws.biass(dc,1e3)

load("core_inference_SNPIZ_mcmcEachCount.rda")

mcmc.sarws.diagnostics(mc, 1e3, 2e2)
