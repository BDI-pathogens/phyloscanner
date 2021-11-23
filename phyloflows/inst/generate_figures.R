library(rstan)
library(data.table)
library(ggplot2)
library(viridis)
library(coda)
# library(phyloflows)
# library(devtools)
# install_github("BDI-pathogens/phyloscanner/phyloflows",dependencies=T)

# load
dir.create('~/gpxx_result')
setwd('~/gpxx_result')
load('~/gpxx.rda')
load('~/ageanalysis/input.rda')
load('~/input_dobs.rda')

# change formats
mc <- list()
params <- extract(fit,pars=c(paste0('llbd[',rep(1:data.fit$N_group,each=data.fit$N_per_group),',',rep(1:data.fit$N_per_group,times=data.fit$N_group),']'),
                             paste0('lxi_pair[',rep(1:data.fit$N_group,each=data.fit$N_per_group),',',rep(1:data.fit$N_per_group,times=data.fit$N_group),']')))
params <- do.call(rbind,params)
log_lambda <- params[1:(nrow(params)/2),] - params[(1+nrow(params)/2):nrow(params),]
mc$pars$LOG_LAMBDA <- t(log_lambda)
rm(fit)
gc()

# aggregate by gender and source age
aggregate.file	<- 'aggregated_bygendersourceage.csv'
daggregateTo	<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
daggregateTo[, TR_TARGETCAT:= paste0(gsub('^(.)\\:(.)\\:(.+)$','\\2',TR_TRM_CATEGORY),':',
                                     gsub('^(.)\\:(.)\\:(.+)$','\\3',TR_TRM_CATEGORY))]
daggregateTo[, REC_TARGETCAT:= paste0(gsub('^(.)\\:(.)\\:(.+)$','\\2',REC_TRM_CATEGORY))]
set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
control			<- list(	burnin.p=0,
                   thin=1,
                   regex_pars='*',
                   outfile=aggregate.file)

source.attribution.mcmc.aggregateToTarget(mcmc.file=NULL,
                                          mc=mc,
                                          daggregateTo=daggregateTo,
                                          control=control)

control			<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                   flowratios= list( ))
pars <- data.table(read.csv(aggregate.file))
pars <- pars[VARIABLE=='PI']
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
tmp <- pars[,list(SUM=sum(VALUE)),by=c('TR_SEX','SAMPLE')]
pars <- merge(pars, tmp, by=c('TR_SEX','SAMPLE'))
pars[,VALUE:=VALUE/SUM]
z        <- pars[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
tmp <- pars[,list(P=c('mean','sd'),Q=c(mean(VALUE),sd=sd(VALUE))),by=c('TR_TARGETCAT','REC_TARGETCAT')]
z <- rbind(z,tmp)
z        <- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
setkey(z, TR_TARGETCAT, REC_TARGETCAT )
z[, STAT:='flows']
ans        <- copy(z)
gc()
write.csv(ans, row.names=FALSE,file=gsub('.csv','_normbygender.csv',aggregate.file))



# aggregated by age genders
aggregate.file	<- 'aggregated_byagegender.csv'
daggregateTo	<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
daggregateTo[, TR_TARGETCAT:= paste0(gsub('^(.)\\:(.)\\:(.+)$','\\2',TR_TRM_CATEGORY),':',
                                     gsub('^(.)\\:(.)\\:(.+)$','\\3',TR_TRM_CATEGORY))]
daggregateTo[, REC_TARGETCAT:= paste0(gsub('^(.)\\:(.)\\:(.+)$','\\2',REC_TRM_CATEGORY),':',
                                      gsub('^(.)\\:(.)\\:(.+)$','\\3',REC_TRM_CATEGORY))]
set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)

control	<- list(burnin.p=0,
                thin=1,
                regex_pars='*',
                outfile=aggregate.file)

source.attribution.mcmc.aggregateToTarget(mcmc.file=NULL,
                                          mc=mc,
                                          daggregateTo=daggregateTo,
                                          control=control)

control			<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                   flowratios= list( ))
pars <- data.table(read.csv(aggregate.file))
pars <- pars[VARIABLE=='PI']
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]  
tmp <- pars[,list(SUM=sum(VALUE)),by=c('TR_SEX','SAMPLE')]
pars <- merge(pars, tmp, by=c('TR_SEX','SAMPLE'))
pars[,VALUE:=VALUE/SUM]
z  <- pars[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
tmp <- pars[,list(P=c('mean','sd'),Q=c(mean(VALUE),sd=sd(VALUE))),by=c('TR_TARGETCAT','REC_TARGETCAT')]
z <- rbind(z,tmp)
z <- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
setkey(z, TR_TARGETCAT, REC_TARGETCAT )
z[, STAT:='flows']
ans        <- copy(z)
gc()
write.csv(ans, row.names=FALSE,file=gsub('.csv','_normbygender.csv',aggregate.file))


control			<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                   flowratios= list( ))

pars <- data.table(read.csv(aggregate.file))
pars <- pars[VARIABLE=='PI']

pars[,REC_SEX:=gsub('^(.)\\:(.+)$','\\1',REC_TARGETCAT)]
pars[,REC_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',REC_TARGETCAT)]
tmp <- pars[,list(SUM=sum(VALUE)),by=c('REC_SEX','REC_AGE_AT_MID_C','SAMPLE')]
pars <- merge(pars, tmp, by=c('REC_SEX','REC_AGE_AT_MID_C','SAMPLE'))
pars[,VALUE:=VALUE/SUM]
z        <- pars[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles,na.rm=TRUE))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
tmp <- pars[,list(P=c('mean','sd'),Q=c(mean(VALUE),sd=sd(VALUE))),by=c('TR_TARGETCAT','REC_TARGETCAT')]
z <- rbind(z,tmp)
z        <- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
setkey(z, TR_TARGETCAT, REC_TARGETCAT )
z[, STAT:='flows']
ans        <- copy(z)
gc()
write.csv(ans, row.names=FALSE,file=gsub('.csv','_normbygender_recage.csv',aggregate.file))

control			<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                   flowratios= list( ))
pars <- data.table(read.csv(aggregate.file))
pars <- pars[VARIABLE=='PI']
pars[,TR_COMM_TYPE:=gsub('^(.)\\:(.)\\:(.+)$','\\1',TR_TARGETCAT)]  
pars[,TR_SEX:=gsub('^(.)\\:(.)\\:(.+)$','\\2',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.)\\:(.+)$','\\3',TR_TARGETCAT)]
tmp <- pars[,list(SUM=sum(VALUE)),by=c('TR_COMM_TYPE','TR_SEX','TR_AGE_AT_MID_C','SAMPLE')]
pars <- merge(pars, tmp, by=c('TR_COMM_TYPE','TR_SEX','TR_AGE_AT_MID_C','SAMPLE'))
pars[,VALUE:=VALUE/SUM]
z        <- pars[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
tmp <- pars[,list(P=c('mean','sd'),Q=c(mean(VALUE),sd=sd(VALUE))),by=c('TR_TARGETCAT','REC_TARGETCAT')]
z <- rbind(z,tmp)
z        <- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
setkey(z, TR_TARGETCAT, REC_TARGETCAT )
z[, STAT:='flows']
ans        <- copy(z)
gc()
write.csv(ans, row.names=FALSE,file=gsub('.csv','_normbygender_sourceage.csv',aggregate.file))

# plots
tmp <- data.table(read.csv('aggregated_byagegender_normbygender.csv'))
tmp[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
tmp[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
tmp[,REC_SEX:=gsub('^(.)\\:(.+)$','\\1',REC_TARGETCAT)]
tmp[,REC_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',REC_TARGETCAT)]
tmp[,TR_AGE:=as.numeric(substr(TR_AGE_AT_MID_C,1,2))+0.5]
tmp[,REC_AGE:=as.numeric(substr(REC_AGE_AT_MID_C,1,2))+0.5]

load("~/RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
rtr <- copy(rtr3)
rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)

ggplot(tmp[TR_SEX=='F',], aes(x=TR_AGE, y=REC_AGE))+
  geom_tile(aes(fill = M)) +
  scale_fill_viridis(breaks=seq(0,0.005,by=0.001),
                     limits=c(0,0.005),labels = function(x) format(x, scientific = TRUE))+
  scale_x_continuous(expand = c(0,0), limits = c(15.5,50.5))+
  scale_y_continuous(expand = c(0,0), limits = c(15.5,50.5))+
  labs(x='\n ages of female sources', y='ages of male recipients \n', 
       fill='flows')+
  theme_classic()+
  geom_point(data=rtr[TR_SEX=='F'], aes(x=TR_AGE_AT_MID,y=REC_AGE_AT_MID,z=NULL), color='red' )+
  geom_density_2d(data=rtr[TR_SEX=='F'], aes(x=TR_AGE_AT_MID,y=REC_AGE_AT_MID,z=NULL), color='red',bins=10)+
  theme(legend.position = 'bottom',
        legend.box = 'horizontal', legend.direction = 'horizontal') +
  geom_abline(intercept = 0, slope = 1, color='white', linetype=2)+
  theme(text = element_text(size = 20),legend.text = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))
ggsave(file='median_flow_plot_fm.pdf',width = 6.5, height = 7.5)

ggplot(tmp[TR_SEX=='M',], aes(x=TR_AGE, y=REC_AGE))+
  geom_tile(aes(fill = M)) +
  scale_fill_viridis(breaks=seq(0,0.005,by=0.001),
                     limits=c(0,0.005),labels = function(x) format(x, scientific = TRUE))+
  scale_x_continuous(expand = c(0,0), limits = c(15.5,50.5))+
  scale_y_continuous(expand = c(0,0), limits = c(15.5,50.5))+
  labs(x='\n ages of male sources', y='ages of female recipients \n', 
       fill='flows')+
  theme_classic()+
  geom_point(data=rtr[TR_SEX=='M'], aes(x=TR_AGE_AT_MID,y=REC_AGE_AT_MID,z=NULL), color='red' )+
  geom_density_2d(data=rtr[TR_SEX=='M'], aes(x=TR_AGE_AT_MID,y=REC_AGE_AT_MID,z=NULL), color='red',bins=10)+
  theme(legend.position = 'bottom',
        legend.box = 'horizontal', legend.direction = 'horizontal') +
  geom_abline(intercept = 0, slope = 1, color='white', linetype=2)+
  theme(text = element_text(size = 20),legend.text = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5))
ggsave(file='median_flow_plot_mf.pdf',width = 6.5, height = 7.5)


# 
df <- data.table(read.csv('aggregated_bygendersourceage_normbygender.csv'))
df[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
df[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
df[,TR_AGE:=as.numeric(substr(TR_AGE_AT_MID_C,1,2))+0.5]

ggplot(df[TR_SEX=='M',],aes(TR_AGE,M))+
  geom_point()+
  geom_errorbar(aes(ymin=CL,ymax=CU))+
  scale_x_continuous(expand = c(0,0),limits = c(15,50), breaks = seq(15,50,5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1), breaks = seq(0,0.1,0.02), labels = scales::percent) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x='\n ages of male sources', y='proportions \n')
ggsave('sources_mf.pdf',width=6, height=5)

ggplot(df[TR_SEX=='F',],aes(TR_AGE,M))+
  geom_point()+
  geom_errorbar(aes(ymin=CL,ymax=CU)) +
  scale_x_continuous(expand = c(0,0),limits = c(15,50), breaks = seq(15,50,5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1), breaks = seq(0,0.1,0.02), labels = scales::percent) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x='\n ages of female sources', y='proportions \n')
ggsave('source_fm.pdf',width=6, height=5)


# 
pars <- data.table(read.csv('aggregated_byagegender_normbygender_sourceage.csv'))
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
pars[,REC_SEX:=gsub('^(.)\\:(.+)$','\\1',REC_TARGETCAT)]
pars[,REC_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',REC_TARGETCAT)]
pars[,TR_AGE:=as.numeric(substr(TR_AGE_AT_MID_C,1,2))+0.5]
pars[,REC_AGE:=as.numeric(substr(REC_AGE_AT_MID_C,1,2))+0.5]

df <- data.table()
for(sex in unique(pars$TR_SEX)){
  for (age in unique(pars$TR_AGE)) {
    tmp <- pars[TR_SEX==sex & TR_AGE==age,]
    setkey(tmp, REC_AGE)
    samples <- sample(tmp$REC_AGE,prob = tmp$M, size=1e5, replace = T)
    hpd9 <- HPDinterval(as.mcmc(samples),prob = 0.95)
    hpd5 <- HPDinterval(as.mcmc(samples),prob = 0.5)
    hp <- tmp$REC_AGE[which.max(tmp$M)]
    tmp <- data.table(M=hp, CL=hpd9[1], IL=hpd5[1], CU=hpd9[2], IU=hpd5[2],SEX=sex, AGE=age)
    df <- rbind(df, tmp)
  } 
}

ggplot(df)+
  geom_line(aes(x=AGE,y=M)) + 
  guides(col=F)+
  geom_ribbon(aes(x=AGE,ymin=IL, ymax=IU),alpha=0.6)+
  geom_ribbon(aes(x=AGE,ymin=CL, ymax=CU),alpha=0.2)+
  scale_x_continuous(expand = c(0,0), limits = c(15,50), breaks = seq(15,50,5)) + 
  scale_y_continuous(expand = c(0,0), limits = c(15,50), breaks = seq(15,50,5)) + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'bottom',
        panel.spacing = unit(2, "lines"))+    
  theme(text = element_text(size = 20)) +
  labs(x='ages of sources', y='ages of recipients',fill='genders of sources') +
  facet_grid(.~factor(SEX,c('F','M'),c('female-to-male','male-to-female'))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")
ggsave('recipient_age_by_sourceagegender.pdf',width=10, height=5.5)


# 
pars <- data.table(read.csv('aggregated_byagegender.csv'))
pars <- pars[VARIABLE=='PI']
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
tmp <- pars[,list(SUM=sum(VALUE)),by=c('TR_SEX','TR_AGE_AT_MID_C','SAMPLE')]
pars <- merge(pars, tmp, by=c('TR_SEX','TR_AGE_AT_MID_C','SAMPLE'))
pars[,VALUE:=VALUE/SUM]
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
pars[,REC_SEX:=gsub('^(.)\\:(.+)$','\\1',REC_TARGETCAT)]
pars[,REC_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',REC_TARGETCAT)]
pars[,TR_AGE:=as.numeric(substr(TR_AGE_AT_MID_C,1,2))+0.5]
pars[,REC_AGE:=as.numeric(substr(REC_AGE_AT_MID_C,1,2))+0.5]

df <- data.table()
for (age in unique(pars$TR_AGE)) {
  pars_s <- pars[TR_SEX=='M' & TR_AGE==age,]
  tmp5 <- pars_s[REC_AGE + 5 < age]
  tmp5 <- tmp5[,list(VALUE=sum(VALUE)), by='SAMPLE']
  tmp10 <- pars_s[REC_AGE +10 < age]
  tmp10 <- tmp10[,list(VALUE=sum(VALUE)), by='SAMPLE']
  tmp <- data.table(M=median(tmp5$VALUE),
                    CL=quantile(tmp5$VALUE,probs = 0.025),
                    CU=quantile(tmp5$VALUE,probs = 0.975),
                    SEX='M', AGE=age, GAP=5)
  df <- rbind(df, tmp)
  tmp <- data.table(M=median(tmp10$VALUE),
                    CL=quantile(tmp10$VALUE,probs = 0.025),
                    CU=quantile(tmp10$VALUE,probs = 0.975),
                    SEX='M', AGE=age, GAP=10)
  df <- rbind(df, tmp)
}


ggplot(df[SEX=='M',])+
  geom_line(aes(x=AGE,y=M,col=factor(GAP))) + 
  guides(col=F)+
  geom_ribbon(aes(x=AGE,ymin=CL, ymax=CU,fill=factor(GAP)),alpha=0.6)+
  scale_x_continuous(expand = c(0,0),limits = c(15,50), breaks = seq(15,50,5)) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,1), breaks = seq(0,1,0.2), labels = scales::percent) + 
  scale_fill_manual(values = c('10'='#00008B',
                               '5'='#87CEEB')) +
  scale_color_manual(values = c('10'='#00008B',
                                '5'='#87CEEB'))+
  theme_bw() +    theme(text = element_text(size = 20)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'bottom')+
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x='age of male sources', y='transmissions to \n younger women',
       fill='age differences between \n  male sources and female \n recipients greater than ')
ggsave('agediff_mf.pdf',width=6, height=5)


pars <- data.table(read.csv('aggregated_byagegender.csv'))
pars <- pars[VARIABLE=='PI']
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
tmp <- pars[,list(SUM=sum(VALUE)),by=c('TR_SEX','TR_AGE_AT_MID_C','SAMPLE')]
pars <- merge(pars, tmp, by=c('TR_SEX','TR_AGE_AT_MID_C','SAMPLE'))
pars[,VALUE:=VALUE/SUM]
pars[,TR_SEX:=gsub('^(.)\\:(.+)$','\\1',TR_TARGETCAT)]
pars[,TR_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',TR_TARGETCAT)]
pars[,REC_SEX:=gsub('^(.)\\:(.+)$','\\1',REC_TARGETCAT)]
pars[,REC_AGE_AT_MID_C:=gsub('^(.)\\:(.+)$','\\2',REC_TARGETCAT)]
pars[,TR_AGE:=as.numeric(substr(TR_AGE_AT_MID_C,1,2))+0.5]
pars[,REC_AGE:=as.numeric(substr(REC_AGE_AT_MID_C,1,2))+0.5]

df <- data.table()
for (age in unique(pars$TR_AGE)) {
  pars_s <- pars[TR_SEX=='F' & TR_AGE==age,]
  tmp5 <- pars_s[REC_AGE + 5 < age]
  tmp5 <- tmp5[,list(VALUE=sum(VALUE)), by='SAMPLE']
  tmp10 <- pars_s[REC_AGE +10 < age]
  tmp10 <- tmp10[,list(VALUE=sum(VALUE)), by='SAMPLE']
  tmp <- data.table(M=median(tmp5$VALUE),
                    CL=quantile(tmp5$VALUE,probs = 0.025),
                    CU=quantile(tmp5$VALUE,probs = 0.975),
                    SEX='F', AGE=age, GAP=5)
  df <- rbind(df, tmp)
  tmp <- data.table(M=median(tmp10$VALUE),
                    CL=quantile(tmp10$VALUE,probs = 0.025),
                    CU=quantile(tmp10$VALUE,probs = 0.975),
                    SEX='F', AGE=age, GAP=10)
  df <- rbind(df, tmp)
}
ggplot(df[SEX=='F',])+
  geom_line(aes(x=AGE,y=M,col=factor(GAP))) + 
  guides(col=F)+
  geom_ribbon(aes(x=AGE,ymin=CL, ymax=CU,fill=factor(GAP)),alpha=0.6)+
  scale_x_continuous(expand = c(0,0),limits = c(15,50), breaks = seq(15,50,5)) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,1), breaks = seq(0,1,0.2), labels = scales::percent) + 
  scale_fill_manual(values = c('10'='#B22222',
                               '5'='#FFA07A')) +
  scale_color_manual(values = c('10'='#B22222',
                                '5'='#FFA07A'))+
  theme_bw() +    theme(text = element_text(size = 20)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'bottom')+
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x='age of female sources', y='transmissions to \n younger men',
       fill='age differences between \n  female sources and male \n recipients greater than ')
ggsave('agediff_fm.pdf',width=6, height=5)

