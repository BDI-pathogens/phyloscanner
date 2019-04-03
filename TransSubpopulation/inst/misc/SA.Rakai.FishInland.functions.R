Rakai190327.participitation.differences.betabinomialmodel3<- function()
{
	require(data.table)
	require(rstan)
	require(extraDistr)
	
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.participation.stan.model <- file.path(indir,"180322_glm_participation.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
			
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
	
	#	binarize age, sex
	des[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	des[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
			
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate participation
	dp		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), by=c('AGE_AT_MID_C','SEX','AGE_YOUNG','AGE_MID','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	dp[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(dp,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE_YOUNG','AGE_MID','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(dp)
	tmp$N_COMM	<- length(unique(dp$COMM_NUM_A))
	fit.par 	<- stan(	file = infile.participation.stan.model, 
						data = tmp, 
						iter = 10e3,
						warmup = 5e2,
						cores = 1,
						chains = 1,
						init = list(list(a=0, comm=rep(0,36), sig_comm=1, dispersion=1, male=0, midage=0, female_young=0, male_young=0, inmigrant_young=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','dispersion','sig_comm','female_young','male_young','inmigrant_young','male','midage','inmigrant','trading','fishing')
	any(rhat(fit.par, pars=fit.pars)>1.02)
	any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.par) 
	po              <- po[, colnames(po)[!grepl('p_suc|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_participation_model_marginalposteriors.pdf'), w=7, h=20)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.par, pars= fit.pars, probs = c(0.025, 0.975))$summary,
				file=paste0(outfile.base,'190327_participation_model_summary.csv')
				)
	
	# extract samples for unique strata levels
	nprior		<- 1e3
	dps			<- data.table(CATEGORY=c("aam:F:15-24:0", "aam:F:25-34:1", "aam:M:25-34:1", "abj:F:15-24:1", "abj:F:25-34:0", "abj:M:25-34:0", "abm:F:15-24:1", "abm:F:25-34:1", "abm:F:35+:0", 
					"abm:M:25-34:0", "abm:M:25-34:1", "abm:M:35+:0", "abo:F:25-34:0", "abo:F:25-34:1", "abo:M:25-34:0", "abo:M:35+:0", "adi:F:15-24:0", "adi:F:15-24:1", "adi:F:25-34:0", 
					"adi:F:25-34:1", "adi:F:35+:0", "adi:M:15-24:1", "adi:M:25-34:0", "adi:M:25-34:1", "adi:M:35+:0", "adi:M:35+:1", "aev:F:15-24:0", "aev:F:25-34:0", "aev:F:25-34:1", 
					"aev:F:35+:0", "aev:M:15-24:0", "aev:M:15-24:1", "aev:M:25-34:0", "afr:F:15-24:0", "afr:F:25-34:0", "afr:F:35+:0", "afr:M:25-34:0", "agf:F:15-24:0", "agf:F:15-24:1", 
					"agf:F:25-34:0", "agf:F:35+:0", "agf:M:15-24:0", "agf:M:25-34:0", "agf:M:35+:0", "agv:F:15-24:0", "agv:F:25-34:0", "agv:F:25-34:1", "agv:M:25-34:0", "agv:M:25-34:1", 
					"agw:F:15-24:0", "agw:F:15-24:1", "agw:F:25-34:0", "agw:F:25-34:1", "agw:F:35+:1", "agw:M:25-34:0", "agw:M:35+:0", "ait:F:25-34:1", "ait:F:35+:0", "ait:M:25-34:1", 
					"akh:F:15-24:0", "akh:F:25-34:0", "akh:F:25-34:1", "akh:F:35+:1", "akh:M:15-24:0", "akh:M:25-34:0", "akh:M:35+:0", "akh:M:35+:1", "ald:F:25-34:0", "ald:F:25-34:1", 
					"ald:F:35+:0", "ald:M:25-34:0", "ald:M:25-34:1", "ald:M:35+:0", "ald:M:35+:1", "alr:F:15-24:0", "alr:F:15-24:1", "alr:F:25-34:0", "alr:F:35+:0", "alr:M:25-34:0", 
					"alr:M:25-34:1", "alr:M:35+:0", "amf:F:25-34:0", "amf:F:25-34:1", "amf:F:35+:0", "amf:F:35+:1", "amf:M:15-24:0", 
					"amf:M:25-34:0", "amf:M:35+:0", "aop:F:15-24:0", "aop:F:25-34:1", "aop:M:15-24:0", "aqi:F:35+:0", "aqj:F:15-24:0", "aqj:F:15-24:1", "aqj:F:25-34:0", 
					"aqj:F:25-34:1", "aqj:M:15-24:0", "aqj:M:25-34:0", "aqj:M:35+:0", "aqj:M:35+:1", "asm:F:15-24:0", "asm:F:15-24:1", "asm:F:25-34:0", "asm:F:35+:0", "asm:F:35+:1", 
					"asm:M:15-24:1", "asm:M:25-34:0", "asm:M:25-34:1", "asm:M:35+:0", "aur:F:15-24:1", "aur:F:25-34:0", "aur:F:25-34:1", "aur:M:25-34:1", "ave:F:15-24:1", "ave:F:25-34:0", 
					"ave:F:25-34:1", "ave:F:35+:0", "ave:M:15-24:0", "ave:M:25-34:0", "ave:M:35+:0", "avl:F:15-24:0", "avl:F:25-34:0", "avl:F:35+:0", "avl:M:25-34:0", "avl:M:35+:0", "awa:F:15-24:0", 
					"awa:F:15-24:1", "awa:F:25-34:0", "awa:M:25-34:0", "awa:M:25-34:1", "awa:M:35+:0", "awr:F:25-34:1", "awr:M:15-24:0", "awr:M:25-34:0", "awr:M:35+:0", "awr:M:35+:1", "axn:F:15-24:1", 
					"fno:F:15-24:0", "fno:F:15-24:1", "fno:F:25-34:0", "fno:F:25-34:1", "fno:F:35+:0", "fno:M:15-24:0", "fno:M:25-34:0", "fno:M:25-34:1", "fno:M:35+:0", "fno:M:35+:1", "fpt:F:15-24:0", 
					"fpt:F:15-24:1", "fpt:F:25-34:0", "fpt:F:25-34:1", "fpt:F:35+:0", "fpt:F:35+:1", "fpt:M:15-24:0", "fpt:M:15-24:1", "fpt:M:25-34:0", "fpt:M:25-34:1", "fpt:M:35+:0", "fpt:M:35+:1", 
					"fus:F:15-24:0", "fus:F:15-24:1", "fus:F:25-34:0", "fus:F:25-34:1", "fus:F:35+:0", "fus:F:35+:1", "fus:M:15-24:0", "fus:M:15-24:1", "fus:M:25-34:0", "fus:M:25-34:1", "fus:M:35+:0", 
					"fus:M:35+:1", "fwd:F:15-24:0", "fwd:F:15-24:1", "fwd:F:25-34:0", "fwd:F:25-34:1", "fwd:F:35+:0", "fwd:F:35+:1", "fwd:M:15-24:0", "fwd:M:15-24:1", "fwd:M:25-34:0", "fwd:M:25-34:1", 
					"fwd:M:35+:0", "fwd:M:35+:1", "thd:F:15-24:0", "thd:F:25-34:0", "thd:F:25-34:1", "thd:F:35+:0", "thd:M:25-34:0", "thd:M:35+:0", "thd:M:35+:1", "tpq:F:15-24:1", "tpq:F:25-34:0", 
					"tpq:M:15-24:0", "tpq:M:25-34:0", "tpq:M:25-34:1", "tpq:M:35+:0", "ttp:F:25-34:0", "ttp:F:35+:1", "ttp:M:35+:0", "twr:F:15-24:0", "twr:F:15-24:1", "twr:F:25-34:0", "twr:F:25-34:1", 
					"twr:F:35+:0", "twr:M:15-24:1", "twr:M:25-34:0", "twr:M:25-34:1", "twr:M:35+:0"))
	dps[, COMM_NUM_A:= dps[, gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\1',CATEGORY)]]
	dps[, SEX:= dps[,gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\2',CATEGORY)]]
	dps[, AGE_AT_MID_C:= dps[, gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\3',CATEGORY)]]
	dps[, INMIGRANT:= dps[, as.integer(gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\4',CATEGORY))]]	
	dps[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	dps[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	dps[, MALE:= as.integer(SEX=='M')]
	dps[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	dps[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	dps			<- merge(dps, unique(subset(dp, select=c(COMM_NUM_A,COMM_NUM_B))), by='COMM_NUM_A')
	
	fit.e		<- extract(fit.par)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dps			<- dps[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + inmigrant_young*INMIGRANT*AGE_YOUNG +
								male_young*AGE_YOUNG*MALE + female_young*AGE_YOUNG*(1-MALE) + midage*AGE_MID)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
			},	
			by=c('CATEGORY')]
	dps[, P:= exp(ETA)/(1+exp(ETA))]
	#	fit kernel density to prior samples with Gourierous and Monfort, 2006 macrobetakernel bounded density estimator, and then estimate log density
	#tmp	<- subset(dps, CATEGORY=='aam:F:15-24:0')	
	#tmp2<- bde(tmp$P, dataPointsCache=tmp$P, b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
	#tmp3<- data.table(DX=seq(0,1,0.001), DY=density(tmp2, seq(0,1,0.001)))
	#tmp4<- data.table(DX=tmp$P, DY=density(tmp2, tmp$P))
	#ggplot(tmp) + geom_histogram(aes(x=P, y=stat(density)), bins=50) +
	#		geom_line(data=tmp3, aes(x=DX, y=DY)) +
	#		geom_point(data=tmp4, aes(x=DX, y=DY), colour='red') +
	#		coord_cartesian(xlim=c(0,1))
	require(bde)
	tmp	<- dps[, {
				bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
				list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
			}, by='CATEGORY']
	dps	<- merge(dps, tmp, by=c('CATEGORY','SAMPLE'))
	set(dps, NULL, c('ETA'), NULL)		
	save(dp, dps, fit.par, file=paste0(outfile.base,'190327_participation_model_samples.rda'))
}
Rakai190327.participitation.differences.betabinomialmodel4<- function()
{
	require(data.table)
	require(rstan)
	require(extraDistr)
	
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.participation.stan.model <- file.path(indir,"180322_glm_participation_age6.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	
	#	binarize age, sex
	des[, AGE1:= as.integer(AGE_AT_MID_C=='20-')]
	des[, AGE2:= as.integer(AGE_AT_MID_C=='20-24')]
	des[, AGE3:= as.integer(AGE_AT_MID_C=='25-29')]
	des[, AGE4:= as.integer(AGE_AT_MID_C=='30-34')]
	des[, AGE5:= as.integer(AGE_AT_MID_C=='35-39')]
	des[, AGE6:= as.integer(AGE_AT_MID_C=='40-44')]
	des[, MALE:= as.integer(SEX=='M')]	
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate participation
	dp		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	dp[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(dp,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(dp)
	tmp$N_COMM	<- length(unique(dp$COMM_NUM_A))
	fit.par 	<- stan(	file = infile.participation.stan.model, 
			data = tmp, 
			iter = 10e3,
			warmup = 5e2,
			cores = 1,
			chains = 1,
			init = list(list(a=0, comm=rep(0,36), sig_comm=1, dispersion=1, male=0, female_age1=0, male_age1=0, female_age2=0, male_age2=0, age3=0, age4=0, age5=0, age6=0, inmigrant_age1=0, inmigrant_age2=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','dispersion','sig_comm','female_age1','male_age1','female_age2','male_age2','inmigrant_age1','inmigrant_age2','age3','age4','age5','age6','male','inmigrant','trading','fishing')
	any(rhat(fit.par, pars=fit.pars)>1.02)
	any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.par) 
	po              <- po[, colnames(po)[!grepl('p_suc|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_participation_modelage6_marginalposteriors.pdf'), w=7, h=20)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_modelage6_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.par, pars= fit.pars, probs = c(0.025, 0.975))$summary,
			file=paste0(outfile.base,'190327_participation_modelage6_summary.csv')
	)
	
	# extract samples for unique strata levels
	nprior		<- 1e3
	fit.e		<- extract(fit.par)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dps			<- dp[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + inmigrant_age1*INMIGRANT*AGE1 + inmigrant_age2*INMIGRANT*AGE2+
								male_age1*AGE1*MALE + female_age1*AGE1*(1-MALE) +
								male_age2*AGE2*MALE + female_age2*AGE2*(1-MALE) +
								age3*AGE3 + age4*AGE4 + age5*AGE5 + age6*AGE6)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]), DISPERSION=as.numeric(fit.e$dispersion[tmp]))
			},	
			by=c('CATEGORY','TRIAL','SUC')]
	dps[, P:= exp(ETA)/(1+exp(ETA))]
	dps[, LP:= dbbinom(SUC, TRIAL, alpha= P*DISPERSION, beta= (1-P)*DISPERSION, log=TRUE)]	
	set(dps, NULL, c('TRIAL','SUC','ETA','DISPERSION'), NULL)
	
	save(dp, dps, fit.par, file=paste0(outfile.base,'190327_participation_modelage6_samples.rda'))
}


Rakai190327.sequencing.differences.binomialmodel1<- function()
{
	require(data.table)
	require(rstan)
	require(bayesplot)
	
	opt.exclude.onART.from.denominator	<- 1
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.sequencing.stan.model <- file.path(indir,"180322_glm_sequencing.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
	
	#	binarize age, sex
	des[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	des[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate sequencing
	des		<- subset(des, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
		stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
		des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
	}
	stopifnot( all(unique(sort(des$COMM_NUM_B))==1:36) )
	ds		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))), by=c('AGE_AT_MID_C','SEX','AGE_YOUNG','AGE_MID','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	ds[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(ds,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE_YOUNG','AGE_MID','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(ds)
	tmp$N_COMM	<- length(unique(ds$COMM_NUM_A))
	fit.seq 	<- stan(	file = infile.sequencing.stan.model, 
			data = tmp, 
			iter = 10e3,
			warmup = 5e2,
			cores = 1,
			chains = 1,
			init = list(list(a=0, comm=rep(0,36), sig_comm=1, male=0, age_mid=0, age_young=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','sig_comm','male','age_mid','age_young','inmigrant','trading','fishing')
	any(rhat(fit.seq, pars=fit.pars)>1.02)
	any(neff_ratio(fit.seq, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.seq) 
	po              <- po[, colnames(po)[!grepl('p_suc_logit|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_sequencing_model_marginalposteriors.pdf'), w=7, h=12)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_sequencing_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.seq, pars= fit.pars, probs = c(0.025, 0.975))$summary,
				file=paste0(outfile.base,'190327_sequencing_model_summary.csv')
				)	
	
	# extract samples for unique strata levels
	nprior		<- 1e3
	dss			<- data.table(CATEGORY=c("aam:F:15-24:0", "aam:F:25-34:1", "aam:M:25-34:1", "abj:F:15-24:1", "abj:F:25-34:0", "abj:M:25-34:0", "abm:F:15-24:1", "abm:F:25-34:1", "abm:F:35+:0", 
					"abm:M:25-34:0", "abm:M:25-34:1", "abm:M:35+:0", "abo:F:25-34:0", "abo:F:25-34:1", "abo:M:25-34:0", "abo:M:35+:0", "adi:F:15-24:0", "adi:F:15-24:1", "adi:F:25-34:0", 
					"adi:F:25-34:1", "adi:F:35+:0", "adi:M:15-24:1", "adi:M:25-34:0", "adi:M:25-34:1", "adi:M:35+:0", "adi:M:35+:1", "aev:F:15-24:0", "aev:F:25-34:0", "aev:F:25-34:1", 
					"aev:F:35+:0", "aev:M:15-24:0", "aev:M:15-24:1", "aev:M:25-34:0", "afr:F:15-24:0", "afr:F:25-34:0", "afr:F:35+:0", "afr:M:25-34:0", "agf:F:15-24:0", "agf:F:15-24:1", 
					"agf:F:25-34:0", "agf:F:35+:0", "agf:M:15-24:0", "agf:M:25-34:0", "agf:M:35+:0", "agv:F:15-24:0", "agv:F:25-34:0", "agv:F:25-34:1", "agv:M:25-34:0", "agv:M:25-34:1", 
					"agw:F:15-24:0", "agw:F:15-24:1", "agw:F:25-34:0", "agw:F:25-34:1", "agw:F:35+:1", "agw:M:25-34:0", "agw:M:35+:0", "ait:F:25-34:1", "ait:F:35+:0", "ait:M:25-34:1", 
					"akh:F:15-24:0", "akh:F:25-34:0", "akh:F:25-34:1", "akh:F:35+:1", "akh:M:15-24:0", "akh:M:25-34:0", "akh:M:35+:0", "akh:M:35+:1", "ald:F:25-34:0", "ald:F:25-34:1", 
					"ald:F:35+:0", "ald:M:25-34:0", "ald:M:25-34:1", "ald:M:35+:0", "ald:M:35+:1", "alr:F:15-24:0", "alr:F:15-24:1", "alr:F:25-34:0", "alr:F:35+:0", "alr:M:25-34:0", 
					"alr:M:25-34:1", "alr:M:35+:0", "amf:F:25-34:0", "amf:F:25-34:1", "amf:F:35+:0", "amf:F:35+:1", "amf:M:15-24:0", 
					"amf:M:25-34:0", "amf:M:35+:0", "aop:F:15-24:0", "aop:F:25-34:1", "aop:M:15-24:0", "aqi:F:35+:0", "aqj:F:15-24:0", "aqj:F:15-24:1", "aqj:F:25-34:0", 
					"aqj:F:25-34:1", "aqj:M:15-24:0", "aqj:M:25-34:0", "aqj:M:35+:0", "aqj:M:35+:1", "asm:F:15-24:0", "asm:F:15-24:1", "asm:F:25-34:0", "asm:F:35+:0", "asm:F:35+:1", 
					"asm:M:15-24:1", "asm:M:25-34:0", "asm:M:25-34:1", "asm:M:35+:0", "aur:F:15-24:1", "aur:F:25-34:0", "aur:F:25-34:1", "aur:M:25-34:1", "ave:F:15-24:1", "ave:F:25-34:0", 
					"ave:F:25-34:1", "ave:F:35+:0", "ave:M:15-24:0", "ave:M:25-34:0", "ave:M:35+:0", "avl:F:15-24:0", "avl:F:25-34:0", "avl:F:35+:0", "avl:M:25-34:0", "avl:M:35+:0", "awa:F:15-24:0", 
					"awa:F:15-24:1", "awa:F:25-34:0", "awa:M:25-34:0", "awa:M:25-34:1", "awa:M:35+:0", "awr:F:25-34:1", "awr:M:15-24:0", "awr:M:25-34:0", "awr:M:35+:0", "awr:M:35+:1", "axn:F:15-24:1", 
					"fno:F:15-24:0", "fno:F:15-24:1", "fno:F:25-34:0", "fno:F:25-34:1", "fno:F:35+:0", "fno:M:15-24:0", "fno:M:25-34:0", "fno:M:25-34:1", "fno:M:35+:0", "fno:M:35+:1", "fpt:F:15-24:0", 
					"fpt:F:15-24:1", "fpt:F:25-34:0", "fpt:F:25-34:1", "fpt:F:35+:0", "fpt:F:35+:1", "fpt:M:15-24:0", "fpt:M:15-24:1", "fpt:M:25-34:0", "fpt:M:25-34:1", "fpt:M:35+:0", "fpt:M:35+:1", 
					"fus:F:15-24:0", "fus:F:15-24:1", "fus:F:25-34:0", "fus:F:25-34:1", "fus:F:35+:0", "fus:F:35+:1", "fus:M:15-24:0", "fus:M:15-24:1", "fus:M:25-34:0", "fus:M:25-34:1", "fus:M:35+:0", 
					"fus:M:35+:1", "fwd:F:15-24:0", "fwd:F:15-24:1", "fwd:F:25-34:0", "fwd:F:25-34:1", "fwd:F:35+:0", "fwd:F:35+:1", "fwd:M:15-24:0", "fwd:M:15-24:1", "fwd:M:25-34:0", "fwd:M:25-34:1", 
					"fwd:M:35+:0", "fwd:M:35+:1", "thd:F:15-24:0", "thd:F:25-34:0", "thd:F:25-34:1", "thd:F:35+:0", "thd:M:25-34:0", "thd:M:35+:0", "thd:M:35+:1", "tpq:F:15-24:1", "tpq:F:25-34:0", 
					"tpq:M:15-24:0", "tpq:M:25-34:0", "tpq:M:25-34:1", "tpq:M:35+:0", "ttp:F:25-34:0", "ttp:F:35+:1", "ttp:M:35+:0", "twr:F:15-24:0", "twr:F:15-24:1", "twr:F:25-34:0", "twr:F:25-34:1", 
					"twr:F:35+:0", "twr:M:15-24:1", "twr:M:25-34:0", "twr:M:25-34:1", "twr:M:35+:0"))
	dss[, COMM_NUM_A:= dss[, gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\1',CATEGORY)]]
	dss[, SEX:= dss[,gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\2',CATEGORY)]]
	dss[, AGE_AT_MID_C:= dss[, gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\3',CATEGORY)]]
	dss[, INMIGRANT:= dss[, as.integer(gsub('^(.+)\\:(.)\\:(.+)\\:(.)$','\\4',CATEGORY))]]	
	dss[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	dss[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	dss[, MALE:= as.integer(SEX=='M')]
	dss[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	dss[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	dss			<- merge(dss, unique(subset(ds, select=c(COMM_NUM_A,COMM_NUM_B))), by='COMM_NUM_A')
	
	fit.e		<- extract(fit.seq)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dss			<- dss[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + age_young*AGE_YOUNG + age_mid*AGE_MID)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
			},	
			by=c('CATEGORY')]
	dss[, P:= exp(ETA)/(1+exp(ETA))]
	#	fit kernel density to prior samples with Gourierous and Monfort, 2006 macrobetakernel bounded density estimator, and then estimate log density
	require(bde)
	tmp	<- dss[, {
				bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
				list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
			}, by='CATEGORY']
	dss	<- merge(dss, tmp, by=c('CATEGORY','SAMPLE'))
	set(dss, NULL, c('ETA'), NULL)		
	save(ds, dss, fit.seq, file=paste0(outfile.base,'190327_sequencing_model_samples.rda'))
}

Rakai190327.sequencing.differences.binomialmodel2<- function()
{
	require(data.table)
	require(rstan)
	require(bayesplot)
	
	opt.exclude.onART.from.denominator	<- 1
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.sequencing.stan.model <- file.path(indir,"180322_glm_sequencing_age6.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
		
	#	binarize age, sex
	des[, AGE1:= as.integer(AGE_AT_MID_C=='20-')]
	des[, AGE2:= as.integer(AGE_AT_MID_C=='20-24')]
	des[, AGE3:= as.integer(AGE_AT_MID_C=='25-29')]
	des[, AGE4:= as.integer(AGE_AT_MID_C=='30-34')]
	des[, AGE5:= as.integer(AGE_AT_MID_C=='35-39')]
	des[, AGE6:= as.integer(AGE_AT_MID_C=='40-44')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate sequencing
	des		<- subset(des, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
		stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
		des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
	}
	stopifnot( all(unique(sort(des$COMM_NUM_B))==1:36) )
	ds		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))), by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	ds[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(ds,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(ds)
	tmp$N_COMM	<- length(unique(ds$COMM_NUM_A))
	fit.seq 	<- stan(	file = infile.sequencing.stan.model, 
			data = tmp, 
			iter = 10e3,
			warmup = 5e2,
			cores = 1,
			chains = 1,
			init = list(list(a=0, comm=rep(0,36), sig_comm=1, male=0, age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','sig_comm','male','age1','age2','age3','age4','age5','age6','inmigrant','trading','fishing')
	any(rhat(fit.seq, pars=fit.pars)>1.02)
	any(neff_ratio(fit.seq, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.seq) 
	po              <- po[, colnames(po)[!grepl('p_suc_logit|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_sequencing_modelage6_marginalposteriors.pdf'), w=7, h=12)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_sequencing_modelage6_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.seq, pars= fit.pars, probs = c(0.025, 0.975))$summary,
				file=paste0(outfile.base,'190327_sequencing_modelage6_summary.csv')
				)	
	
	# extract samples for unique strata levels
	nprior		<- 1e3
	fit.e		<- extract(fit.seq)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dss			<- ds[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + 
								age1*AGE1 + age2*AGE2 + age3*AGE3 + age4*AGE4 + age5*AGE5 + age6*AGE6)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
			},	
			by=c('CATEGORY','TRIAL','SUC')]
	dss[, P:= exp(ETA)/(1+exp(ETA))]
	dss[, LP:= dbinom(SUC, TRIAL, prob=P, log=TRUE)]	
	set(dss, NULL, c('TRIAL','SUC','ETA'), NULL)
	
	save(ds, dss, fit.seq, file=paste0(outfile.base,'190327_sequencing_modelage6_samples.rda'))
}



Rakai190327.analysispipeline.age3model<- function(infile.inference=NULL, infile.participation.prior.samples=NULL, infile.sequencing.prior.samples=NULL, opt=NULL)
{
	require(data.table)	
	require(TransSubpopulation)
	
	#
	#	input args
	#
	if(is.null(opt))
	{
		opt									<- list()
		opt$adjust.sequencing.bias			<- 0
		opt$adjust.participation.bias		<- 1
		opt$migration.def.code				<- '24'
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
	}
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
	}
	if(is.null(infile.participation.prior.samples))
	{
		infile.participation.prior.samples	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_participation_model_samples.rda"	
	}
	if(is.null(infile.sequencing.prior.samples))
	{
		infile.sequencing.prior.samples		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_sequencing_model_samples.rda"
	}
	cat('\ninfile.inference=',infile.inference)
	cat('\ninfile.participation.prior.samples=',infile.participation.prior.samples)
	cat('\ninfile.sequencing.prior.samples=',infile.sequencing.prior.samples)
	cat('\nopt=',unlist(opt))			
	indir					<- dirname(infile.inference)	
	outfile.base			<- gsub('data_with_inmigrants.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	prepare data on observed transmission flows
	#
	#	subset to variables needed, using RTR3	
	rtr	<- copy(rtr3)
	if(opt$migration.def.code=='06')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_05YR, REC_INMIGRATE_05YR, TR_COMM_NUM_A_MIG_05YR')
		setnames(rtr, 'TR_INMIGRATE_05YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_05YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_05YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='12')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_1YR, REC_INMIGRATE_1YR, TR_COMM_NUM_A_MIG_1YR')
		setnames(rtr, 'TR_INMIGRATE_1YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_1YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_1YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='24')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_2YR, REC_INMIGRATE_2YR, TR_COMM_NUM_A_MIG_2YR')
		setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_2YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='36')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_3YR, REC_INMIGRATE_3YR, TR_COMM_NUM_A_MIG_3YR')
		setnames(rtr, 'TR_INMIGRATE_3YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_3YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_3YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='48')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_4YR, REC_INMIGRATE_4YR, TR_COMM_NUM_A_MIG_4YR')
		setnames(rtr, 'TR_INMIGRATE_4YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_4YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_4YR', 'TR_COMM_NUM_A_MIG')
	}
	
	rtr	<- subset(rtr, select=c(	'PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
							'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE',
							'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
							'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE'))
	
	# inmigrant status
	rtr[, TR_INMIGRANT:= as.integer(TR_INMIGRATE!='resident')]
	rtr[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
	set(rtr, NULL, 'TR_COMM_NUM_A_MIG', rtr[, gsub('[0-9]+','',TR_COMM_NUM_A_MIG)])
	#	set unknown origin to either fishing or inland
	tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
	if(opt$set.missing.migloc.to.inland)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
	}		
	if(opt$set.missing.migloc.to.fishing)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fisherfolk')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
	}
	
			
	# add age 
	rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
	rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
	# impute age
	tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
	set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
	tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
	set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
	# fixup from latest surveillance data
	set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
	set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
	set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
	set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
	#	stratify age
	rtr[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	rtr[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
	stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
			
	#	build category to match with sampling data tables 
	rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_NUM_A,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	#	build transmission flow category 
	rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_NUM_A_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	
	#
	#	calculate observed number of transmissions
	#
	dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
	setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY )	
	dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
	
	#	load samples from prior
	load(infile.participation.prior.samples)
	load(infile.sequencing.prior.samples)
	dprior	<- merge(dps, dss, by=c('CATEGORY','SAMPLE'))
	if(opt$adjust.sequencing.bias==0)
	{
		set(dprior, NULL, 'P.y', 1)
		set(dprior, NULL, 'LP.y', 0)
	}
	if(opt$adjust.participation.bias==0)
	{
		set(dprior, NULL, 'P.x', 1)
		set(dprior, NULL, 'LP.x', 0)
	}
	dprior[, P:= P.x*P.y]		# multiply participation and sequencing probabilities 
	dprior[, LP:= LP.x+LP.y]	# add log posterior densities
	set(dprior, NULL, c('P.x','LP.x','P.y','LP.y'), NULL)
	setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')

	#	run MCMC
	mcmc.file	<- paste0(outfile.base,"samcmc190327_nsweep1e5_opt",paste0(opt, collapse=''),".rda")
	control		<- list(	seed=42, 
							mcmc.n=238000*1, 
							verbose=0, 
							outfile=mcmc.file)
	source.attribution.mcmc(dobs, dprior, control=control)
	
	#	run diagnostics
	#mcmc.file	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_SAMCMCv190327_mcmc.rda"
	#mcmc.file	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_samcmc190327_nsweep1e5_opt11101.rda'
	control		<- list(	burnin.p=0.05, 
							regex_pars='*', 
							credibility.interval=0.95, 
							pdf.plot.all.parameters=FALSE, 
							pdf.plot.n.worst.case.parameters=10, 
							pdf.height.per.par=1.2, 
							outfile.base=gsub('\\.rda','',mcmc.file))
	source.attribution.mcmc.diagnostics(mcmc.file, control=control)
	
	#	aggregate MCMC output to fish<->inland
	aggregate.file	<- gsub('\\.rda','_aggregatedFishInland.csv',mcmc.file)
	daggregateTo	<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
	daggregateTo[, TR_TARGETCAT:= gsub('^e$','external',gsub('^f$','fishing',gsub('^a|i|t$','inland',substring(daggregateTo$TR_TRM_CATEGORY,1,1))))]
	daggregateTo[, REC_TARGETCAT:= gsub('^e$','external',gsub('^f$','fishing',gsub('^a|i|t$','inland',substring(daggregateTo$REC_TRM_CATEGORY,1,1))))]
	set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)	
	control			<- list(	burnin.p=0.05, 
								thin=NA_integer_, 
								regex_pars='*', 
								outfile=aggregate.file)
	source.attribution.mcmc.aggregateToTarget(mcmc.file, daggregateTo, control=control)
	
	#	calculate flows sources WAIFM flow_ratio overall		
	control			<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
								flowratios= list( c('inland/fishing', 'inland fishing', 'fishing inland')),
								outfile=gsub('\\.csv','_flowsetc.csv',aggregate.file))
	Rakai190327.getKeyQuantitiesFromAggregatedMCMC(aggregate.file, control)
		
	#	aggregate MCMC output to fish<->inland by gender
	aggregate.file	<- gsub('\\.rda','_aggregatedFishInlandByGender.csv',mcmc.file)
	daggregateTo	<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
	daggregateTo[, TR_TARGETCAT:= paste0(	gsub('^e$','external',gsub('^f$','fishing',gsub('^a|i|t$','inland',substring(daggregateTo$TR_TRM_CATEGORY,1,1)))),
											':',
											gsub('^[a-z]+:([FM]):.*','\\1',daggregateTo$TR_TRM_CATEGORY))]
	daggregateTo[, REC_TARGETCAT:= paste0(	gsub('^e$','external',gsub('^f$','fishing',gsub('^a|i|t$','inland',substring(daggregateTo$REC_TRM_CATEGORY,1,1)))),
											':',
											gsub('^[a-z]+:([FM]):.*','\\1',daggregateTo$REC_TRM_CATEGORY))]
	set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)	
	control		<- list(	burnin.p=0.05, 
							thin=NA_integer_, 
							regex_pars='*', 
							outfile=aggregate.file)
	source.attribution.mcmc.aggregateToTarget(mcmc.file, daggregateTo, control=control)
	
	#	calculate flows sources WAIFM flow_ratio by gender		
	control		<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
							flowratios= list( c('inland:M/fishing:M', 'inland:M fishing:F', 'fishing:M inland:F'), c('inland:F/fishing:F', 'inland:F fishing:M', 'fishing:F inland:M')),
							outfile=gsub('\\.csv','_flowsetc.csv',aggregate.file))
	Rakai190327.getKeyQuantitiesFromAggregatedMCMC(aggregate.file, control)				
}

Rakai190327.getKeyQuantitiesFromAggregatedMCMC<- function(infile, control)
{
	cat('\nReading aggregated MCMC output...')
	pars		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	pars		<- subset(pars, VARIABLE=='PI')
	
	cat('\nComputing flows...')
	#	calculate flows
	z		<- pars[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
	z		<- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_TARGETCAT, REC_TARGETCAT )
	z[, STAT:='flows']
	ans		<- copy(z)		
	gc()
	
	cat('\nComputing WAIFM...')
	#	calculate WAIFM	
	z		<- pars[, list(REC_TARGETCAT=REC_TARGETCAT, VALUE=VALUE/sum(VALUE)), by=c('TR_TARGETCAT','SAMPLE')]
	z		<- z[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
	z		<- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_TARGETCAT, REC_TARGETCAT )
	z[, STAT:='waifm']
	ans		<- rbind(ans,z)	
	gc()
	
	cat('\nComputing sources...')
	#	calculate sources
	z		<- pars[, list(TR_TARGETCAT=TR_TARGETCAT, VALUE=VALUE/sum(VALUE)), by=c('REC_TARGETCAT','SAMPLE')]
	z		<- z[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('TR_TARGETCAT','REC_TARGETCAT')]
	z		<- dcast.data.table(z, TR_TARGETCAT+REC_TARGETCAT~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_TARGETCAT, TR_TARGETCAT )
	z[, STAT:='sources']
	ans		<- rbind(ans,z)	
	gc()
	
	if(length(control$flowratios)>0)
	{
		cat('\nComputing flow ratios...')
		#	calculate transmission flow ratios	
		z		<- copy(pars)
		z[, FLOW:=paste0(TR_TARGETCAT,' ',REC_TARGETCAT)]
		z		<- dcast.data.table(z, SAMPLE~FLOW, value.var='VALUE')	
		set(z, NULL, colnames(z)[!colnames(z)%in%c('SAMPLE',unlist(control$flowratios))], NULL)
		for(ii in seq_along(control$flowratios))
		{
			if(!all(c(control$flowratios[[ii]][2],control$flowratios[[ii]][3])%in%colnames(z)))
				stop('column names in control$flowratios are not in MCMC output')
			set(z, NULL, control$flowratios[[ii]][1], z[[ control$flowratios[[ii]][2] ]] /  z[[ control$flowratios[[ii]][3] ]])
			set(z, NULL, c(control$flowratios[[ii]][2],control$flowratios[[ii]][3]), NULL)
		}
		z		<- melt(z, id.vars='SAMPLE', value.name='VALUE', variable.name='FLOWRATIO_CAT')
		z		<- z[, list(P=names(control$quantiles), Q=unname(quantile(VALUE, p=control$quantiles))), by=c('FLOWRATIO_CAT')]
		z		<- dcast.data.table(z, FLOWRATIO_CAT~P, value.var='Q')
		z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
		z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
		z[, STAT:='flow_ratio']
		ans		<- rbind(ans,z, fill=TRUE)	
		gc()		
	}
	
	cat('\nWriting output to',control$outfile)
	write.csv(ans, row.names=FALSE,file=control$outfile)
}

Rakai190327.analysispipeline.age6model<- function(infile.inference=NULL, opt=NULL)
{
	require(data.table)	
	
	if(is.null(opt))
	{
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$migration.def.code				<- '24'
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
	}
	if(is.null(infile.inference))
	{
		indir				<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
		indir				<- "/rds/general/user/or105/home/WORK/Gates_2014/Rakai"
		infile.inference	<- file.path(indir,"todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
	}	
	indir					<- dirname(infile.inference)	
	outfile.base			<- gsub('data_with_inmigrants.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	prepare data on observed transmission flows
	#
	
	#	subset to variables needed, using RTR3	
	rtr	<- copy(rtr3)
	if(opt$migration.def.code=='06')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_05YR, REC_INMIGRATE_05YR, TR_COMM_NUM_A_MIG_05YR')
		setnames(rtr, 'TR_INMIGRATE_05YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_05YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_05YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='12')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_1YR, REC_INMIGRATE_1YR, TR_COMM_NUM_A_MIG_1YR')
		setnames(rtr, 'TR_INMIGRATE_1YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_1YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_1YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='24')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_2YR, REC_INMIGRATE_2YR, TR_COMM_NUM_A_MIG_2YR')
		setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_2YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='36')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_3YR, REC_INMIGRATE_3YR, TR_COMM_NUM_A_MIG_3YR')
		setnames(rtr, 'TR_INMIGRATE_3YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_3YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_3YR', 'TR_COMM_NUM_A_MIG')
	}
	if(opt$migration.def.code=='48')
	{
		cat('\nSource attribution analysis based on TR_INMIGRATE_4YR, REC_INMIGRATE_4YR, TR_COMM_NUM_A_MIG_4YR')
		setnames(rtr, 'TR_INMIGRATE_4YR', 'TR_INMIGRATE')
		setnames(rtr, 'REC_INMIGRATE_4YR', 'REC_INMIGRATE')
		setnames(rtr, 'TR_COMM_NUM_A_MIG_4YR', 'TR_COMM_NUM_A_MIG')
	}
	
	rtr	<- subset(rtr, select=c(	'PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
					'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE',
					'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
					'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE'))
	
	# inmigrant status
	rtr[, TR_INMIGRANT:= as.integer(TR_INMIGRATE!='resident')]
	rtr[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
	set(rtr, NULL, 'TR_COMM_NUM_A_MIG', rtr[, gsub('[0-9]+','',TR_COMM_NUM_A_MIG)])
	#	set unknown origin to either fishing or inland
	tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
	if(opt$set.missing.migloc.to.inland)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
	}		
	if(opt$set.missing.migloc.to.fishing)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fisherfolk')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
	}
	
	# add age 
	rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
	rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
	# impute age
	tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
	set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
	tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
	set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
	# fixup from latest surveillance data
	#de[RID=="C036808",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[TR_RID=="C036808",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
	#de[RID=="G036802",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[REC_RID=="G036802",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
	#de[RID=="H103745",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[REC_RID=="H103745",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
	set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
	#	stratify age
	rtr[, TR_AGE_AT_MID_C:= cut(TR_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	rtr[, REC_AGE_AT_MID_C:= cut(REC_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
	stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
	
	#	build category to match with sampling data tables 
	rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_NUM_A,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	#	build transmission flow category 
	rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_NUM_A_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	
	#
	#	calculate observed number of transmissions
	#
	dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
	setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY )	
	dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
	
	#	load samples from prior
	infile.participation.prior.samples	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_participation_modelage6_samples.rda"
	infile.sequencing.prior.samples		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_sequencing_modelage6_samples.rda"
	load(infile.participation.prior.samples)
	load(infile.sequencing.prior.samples)
	dprior	<- merge(dps, dss, by=c('CATEGORY','SAMPLE'))
	dprior[, P:= P.x*P.y]		# multiply participation and sequencing probabilities 
	dprior[, LP:= LP.x+LP.y]	# add log posterior densities
	set(dprior, NULL, c('P.x','LP.x','P.y','LP.y'), NULL)
	setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')
	
	#	run MCMC
	control	<- list(seed=42, mcmc.n=238000, verbose=0, outfile="~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_SAMCMCv190327_age6_mcmc.rda")
	mcmc.core.inference(dobs, dprior, control=control)
	
	#	diagnostics
	load(control$outfile)
	
	#	acceptance rate per site
	da	<- subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, list(ACC_RATE=mean(ACCEPT)), by='PAR_ID']
	setnames(da, 'PAR_ID', 'UPDATE_ID')
	tmp	<- mc$dl[, list(N_TRM_CAT_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
	da	<- merge(da, tmp, by='UPDATE_ID')
	ggplot(da, aes(x=N_TRM_CAT_PAIRS, y=ACC_RATE)) + geom_point()
	
	#	traces
	tmp	<- mc$pars$PI
	colnames(tmp)<- paste0('PI-', 1:ncol(tmp))
	p	<- mcmc_trace(tmp, pars=colnames(tmp), facet_args = list(ncol = 1))
	pdf(file=paste0(outfile.base,'samodel_marginaltraces.pdf'), w=7, h=300)
	p
	dev.off()
	
	#	effective sample size on target parameter
	require(coda)
	tmp	<- mcmc(mc$pars$PI)
	de	<- data.table(PI= 1:ncol(tmp), NEFF=as.numeric(effectiveSize(tmp)))
	ggplot(de, aes(x=PI, y=NEFF)) + geom_point()	
	
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
}

Rakai190327.sensitivitypipeline.age3model<- function()
{
	
	#
	#	make data set
	#
	if(0)
	{
		infiles		<- c("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_withmetadata.rda"
		)
		for(ii in seq_along(infiles))
		{
			infile.inference	<- infiles[[ii]]
			RakaiFull.phylogeography.181006.gender.mobility.data(infile.inference)		
		}
	}
	
	#
	#	core inference
	#
	if(0)
	{
		indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'
		indir		<- '/rds/general/user/or105/home/WORK/Gates_2014/Rakai'
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$migration.def.code				<- '24'
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infiles		<- c( "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_data_with_inmigrants.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_data_with_inmigrants.rda"												
		)
		for(ii in seq_along(infiles))
		{
			infile.inference	<- file.path(indir, infiles[[ii]])
			RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		}
		#
		opt									<- list()
		opt$adjust.sequencing.bias			<- 0
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 0
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 0
		opt$adjust.participation.bias		<- 0
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 0
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 1
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
	}
	
	if(0)
	{
		#
		#	extract target variables of interest
		#
		indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'	
		infiles		<- c(	"todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc_11110.rda"
				, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc_11101.rda"
				, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc_11001.rda"
				, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc_10101.rda"
				, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc_01101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_inference_mcmc_11101.rda"
				, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_inference_mcmc_11101.rda"
		)
		for(ii in seq_along(infiles))
		{
			print(ii)
			infile.inference	<- file.path(indir,infiles[[ii]])
			RakaiFull.phylogeography.181006.flows.fishinlandgender(infile.inference)
			#RakaiFull.phylogeography.181006.mcmc.assess(infile.inference)
			#RakaiFull.phylogeography.181006.flows.fishinlandmigrant(infile.inference)
			#RakaiFull.phylogeography.181006.flows.fishinlandmigrantgender(infile.inference)
			#RakaiFull.phylogeography.181006.flows.netflows(infile.inference)
		}
	}
	
	#
	#	predict flows between subdistrict areas
	#
	if(1)
	{
		indir							<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
		indir							<- '/rds/general/user/or105/home/WORK/Gates_2014/Rakai'		
		infile.inference.mcmc			<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_core_inference_mcmc_11101.rda")
		mc.thin							<- 5		
		RakaiFull.phylogeography.181006.gender.mobility.thin.mcmc(infile.inference.mcmc, mc.thin=mc.thin)
		
		infile.inference.mcmc.thinned	<- gsub('\\.rda',paste0('_thinned',mc.thin,'sweep.rda'),infile.inference.mcmc)
		infile.inference.data			<- file.path(indir, "RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")						
		infile.subdistricts				<- file.path(indir, "Subdistrict_Data_fromMapsAndSurvey.rda")		
		RakaiFull.phylogeography.181006.predict.areaflows(indir=indir, infile.inference.data=infile.inference.data, infile.inference.mcmc.thinned=infile.inference.mcmc.thinned, infile.subdistricts=infile.subdistricts)						
	}
	
}