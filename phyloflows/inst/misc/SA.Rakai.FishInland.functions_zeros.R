Rakai190910.participitation.age3.mp5<- function()
{
  require(data.table)
  require(rstan)
  require(extraDistr)
  require(bayesplot)
  
  indir <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
  infile.data <- file.path(indir,"190327_sampling_by_gender_age.rda")
  indir.stan <- '~/git/phyloscanner/phyloflows/inst/misc'
  infile.participation.stan.model <- file.path(indir.stan,"190910_glm_participation_age3_mp5.stan")
  outfile.base <- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
  
  #	set up variables for STAN
  load(infile.data)
  des <- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
  
  # remove individuals in the communities where no sequences were obtained successfully
  tmp <- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
  des <- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
  
  # set no INMIGRANT to 0
  set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
  
  # age group
  des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
  
  #	binarize age, sex
  des[, AGE1:= as.integer(AGE_AT_MID_C=='15-24')]
  des[, AGE2:= as.integer(AGE_AT_MID_C=='25-34')]
  des[, MALE:= as.integer(SEX=='M')]
  
  #	binarize community type
  des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
  
  #	vanilla community IDs
  tmp <- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
  des <- merge(des, tmp, by='COMM_NUM_A')
  
  #	aggregate participation
  dp <- des[, 
		  		list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), 
             	by=c('COMM_NUM_A','COMM_NUM_B','AGE_AT_MID_C','SEX','AGE1','AGE2','INMIGRANT','MALE','COMM_TYPE_F')]
  dp[, CATEGORY:= paste0(COMM_TYPE_F,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
  
  #	run STAN 
  tmp <- as.list(subset(dp,select=c('TRIAL','SUC','MALE','AGE1','AGE2','INMIGRANT','COMM_TYPE_F')))
  tmp$N <- nrow(dp)  
  fit.par <- stan(	file = infile.participation.stan.model, 
                    data = tmp, 
                    iter = 15e3,
                    warmup = 5e2,
                    cores = 1,
                    chains = 1,
                    control=list(max_treedepth=15),
                    init = list(list(dispersion=1, a=0, fishing=0, male_age1_inmigrant=0, male_age2_inmigrant=0, male_age3_inmigrant=0, female_age1_inmigrant=0, 
									female_age2_inmigrant=0, female_age3_inmigrant=0,
									male_age1_resident=0, male_age2_resident=0, male_age3_resident=0, female_age1_resident=0, 
									female_age2_resident=0)))  
  # assess convergence
  fit.pars <- c("a","fishing","male_age1_inmigrant","male_age2_inmigrant","male_age3_inmigrant","female_age1_inmigrant","female_age2_inmigrant","female_age3_inmigrant","male_age1_resident","male_age2_resident","male_age3_resident","female_age1_resident","female_age2_resident")
  any(rhat(fit.par, pars=fit.pars)>1.02)
  any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
  po <- as.matrix(fit.par) 
  po <- po[, colnames(po)[!grepl('p_suc|lp__',colnames(po))]]  
  # plot traces
  p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
  pdf(file=paste0(outfile.base,'190910_participation_age3_mp5_marginaltraces.pdf'), w=7, h=20)
  p
  dev.off()
  # plot marginal posteriors
  p <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + geom_vline(xintercept = 0)
  pdf(file=paste0(outfile.base,'190910_participation_age3_mp5_marginalposteriors.pdf'), w=7, h=20)
  p
  dev.off()
  
  # summary csv file
  write.csv(	summary(fit.par, pars= fit.pars, probs = c(0.025, 0.975))$summary,
              	file=paste0(outfile.base,'190910_participation_age3_mp5_summary.csv')
 				)
  
  
  # extract samples for unique strata levels
  nprior <- 1e3
  dps <- dp  
  set(dps, NULL, c('COMM_NUM_A','COMM_NUM_B','TRIAL','SUC'), NULL)
  dps <- unique(dps)  
  fit.e <- extract(fit.par)
  set.seed(42)
  tmp <- sample(length(fit.e$a), nprior)
  dps <- dps[,	
               {
                 z<- with(fit.e, a + fishing * COMM_TYPE_F +  
								 male_age1_inmigrant*MALE*AGE1*INMIGRANT + male_age2_inmigrant*MALE*AGE2*INMIGRANT + male_age3_inmigrant*MALE*INMIGRANT +      
								 female_age1_inmigrant*(1-MALE)*AGE1*INMIGRANT + female_age2_inmigrant*(1-MALE)*AGE2*INMIGRANT + female_age3_inmigrant*(1-MALE)*INMIGRANT +
								 male_age1_resident*MALE*AGE1*(1-INMIGRANT) + male_age2_resident*MALE*AGE2*(1-INMIGRANT) + male_age3_resident*MALE*(1-INMIGRANT) +      
								 female_age1_resident*(1-MALE)*AGE1*(1-INMIGRANT) + female_age2_resident*(1-MALE)*AGE2*(1-INMIGRANT)
                 			)
                 list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
               },	
               by=c('CATEGORY')]
  dps[, P:= exp(ETA)/(1+exp(ETA))]
  # estimate log density
  require(bde)
  tmp <- dps[,	{
    				bdest <- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
    				list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
  				}, 
				by='CATEGORY']
  
  dps <- merge(dps, tmp, by=c('CATEGORY','SAMPLE'))
  set(dps, NULL, c('ETA'), NULL)		
  
  save(dp, dps, fit.par, file=paste0(outfile.base,'190910_participation_age3_mp5_samples.rda'))
}

Rakai190910.participitation.modeldev<- function()
{
	require(data.table)
	require(rethinking)
	require(extraDistr)
	require(bayesplot)
	
	indir <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data <- file.path(indir,"190327_sampling_by_gender_age.rda")
	outfile.base <- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des <- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp <- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des <- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
	
	#	binarize age, sex
	des[, AGE1:= as.integer(AGE_AT_MID_C=='15-24')]
	des[, AGE2:= as.integer(AGE_AT_MID_C=='25-34')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	
	#	vanilla community IDs
	tmp <- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des <- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate participation but keep communities as units
	#	think of these as replicate observations
	dp <- des[, 
			list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), 
			by=c('COMM_NUM_A','COMM_NUM_B','AGE_AT_MID_C','SEX','AGE1','AGE2','INMIGRANT','MALE','COMM_TYPE_F')]
	dp[, CATEGORY:= paste0(COMM_TYPE_F,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	# XX model, but Binomial
	mp1 	<- map2stan(
			alist(
					SUC ~ dbinom(TRIAL, p_part),
					logit(p_part) <- a + fishing * COMM_TYPE_F + inmigrant * INMIGRANT + male * MALE +      
							age1 * AGE1 + age2 * AGE2 + 
							age1_female * AGE1 * (1 - MALE) +      
							age1_male * AGE1 * MALE + 
							age3_female * (1 - AGE1 - AGE2) * (1 - MALE) + 
							age3_male * (1 - AGE1 - AGE2) * MALE,
					
					a ~ dnorm(0, 100),
					c(fishing, inmigrant, male, age1, age2, age1_female, age1_male, age3_female, age3_male) ~ dnorm(0,10)
			),
			data=as.data.frame(dp),
			control=list(max_treedepth=15),
			start=list(	a=0, fishing=0, inmigrant=0, male=0, age1=0, age2=0, age1_female=0, age1_male=0, age3_female=0, age3_male=0),
			warmup=5e2, iter=15e3, chains=1, cores=1
		)
	plot( precis(mp1, depth=2, prob=0.95) )
	save(mp1, dp, des, file=paste0(outfile.base,'190910_participation_differences_mp1.rda'))
	#	posterior check
	sims 	<- sim(mp1, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dp, dpp)		
	dpp[, mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
	#	0.26 
	dpp[, COMM_TYPE_F2:= as.character(factor(COMM_TYPE_F, levels=c(0,1), labels=c('inland','fisherfolk')))]
	dpp[, INMIGRANT2:= as.character(factor(INMIGRANT, levels=c(0,1), labels=c('resident','inmigrant')))]
	dpp[, CLASS:= ifelse(SUC<predicted_obs_l95, 'below 95%', ifelse(SUC>predicted_obs_u95, 'above 95%', 'within 95%'))]
	ggplot(dpp, aes(x=interaction(SEX, AGE_AT_MID_C)) ) +
			geom_bar(aes(fill=CLASS)) +
			facet_grid( INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_participation_differences_mp1_pp.pdf'), w=8, h=8)
	
	# basic model, no interactions
	mp2 	<- map2stan(
			alist(
					SUC ~ dbinom(TRIAL, p_part),
					logit(p_part) <- a + fishing * COMM_TYPE_F + inmigrant * INMIGRANT + male * MALE +      
							age1 * AGE1 + age2 * AGE2,					
					a ~ dnorm(0, 100),
					c(fishing, inmigrant, male, age1, age2) ~ dnorm(0,10)
			),
			data=as.data.frame(dp),
			control=list(max_treedepth=15),
			start=list(	a=0, fishing=0, inmigrant=0, male=0, age1=0, age2=0),
			warmup=5e2, iter=15e3, chains=1, cores=1
		)
	plot( precis(mp2, depth=2, prob=0.95) )
	save(mp2, dp, des, file=paste0(outfile.base,'190910_participation_differences_mp2.rda'))
	#	posterior check
	sims 	<- sim(mp2, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dp, dpp)		
	dpp[, mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
	#	0.28 
	dpp[, COMM_TYPE_F2:= as.character(factor(COMM_TYPE_F, levels=c(0,1), labels=c('inland','fisherfolk')))]
	dpp[, INMIGRANT2:= as.character(factor(INMIGRANT, levels=c(0,1), labels=c('resident','inmigrant')))]
	dpp[, CLASS:= ifelse(SUC<predicted_obs_l95, 'below 95%', ifelse(SUC>predicted_obs_u95, 'above 95%', 'within 95%'))]
	ggplot(dpp, aes(x=interaction(SEX, AGE_AT_MID_C)) ) +
			geom_bar(aes(fill=CLASS)) +
			facet_grid( INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_participation_differences_mp2_pp.pdf'), w=8, h=8)
	
	
	
	# model with gender:age interactions
	mp3 	<- map2stan(
			alist(
					SUC ~ dbinom(TRIAL, p_part),
					logit(p_part) <- a + fishing * COMM_TYPE_F + inmigrant * INMIGRANT + 
							male_age1*MALE*AGE1 + male_age2*MALE*AGE2 + male_age3*MALE +      
							female_age1*(1-MALE)*AGE1 + female_age2*(1-MALE)*AGE2,					
					a ~ dnorm(0, 100),
					c(fishing, inmigrant, male_age1, male_age2, male_age3, female_age1, female_age2) ~ dnorm(0,10)
			),
			data=as.data.frame(dp),
			control=list(max_treedepth=15),
			start=list(	a=0, fishing=0, inmigrant=0, male_age1=0, male_age2=0, male_age3=0, female_age1=0, female_age2=0),
			warmup=5e2, iter=15e3, chains=1, cores=1
		)
	plot( precis(mp3, depth=2, prob=0.95) )
	save(mp3, dp, des, file=paste0(outfile.base,'190910_participation_differences_mp3.rda'))
	#	posterior check
	sims 	<- sim(mp3, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dp, dpp)		
	dpp[, mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
	#	0.27 
	dpp[, COMM_TYPE_F2:= as.character(factor(COMM_TYPE_F, levels=c(0,1), labels=c('inland','fisherfolk')))]
	dpp[, INMIGRANT2:= as.character(factor(INMIGRANT, levels=c(0,1), labels=c('resident','inmigrant')))]
	dpp[, CLASS:= ifelse(SUC<predicted_obs_l95, 'below 95%', ifelse(SUC>predicted_obs_u95, 'above 95%', 'within 95%'))]
	ggplot(dpp, aes(x=interaction(SEX, AGE_AT_MID_C)) ) +
			geom_bar(aes(fill=CLASS)) +
			facet_grid( INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_participation_differences_mp3_pp.pdf'), w=8, h=8)
	
	
	
	# model with gender:age:mig_status interactions
	mp4 	<- map2stan(
			alist(
					SUC ~ dbinom(TRIAL, p_part),
					logit(p_part) <- a + fishing * COMM_TYPE_F +  
							male_age1_inmigrant*MALE*AGE1*INMIGRANT + male_age2_inmigrant*MALE*AGE2*INMIGRANT + male_age3_inmigrant*MALE*INMIGRANT +      
							female_age1_inmigrant*(1-MALE)*AGE1*INMIGRANT + female_age2_inmigrant*(1-MALE)*AGE2*INMIGRANT + female_age3_inmigrant*(1-MALE)*INMIGRANT +
							male_age1_resident*MALE*AGE1*(1-INMIGRANT) + male_age2_resident*MALE*AGE2*(1-INMIGRANT) + male_age3_resident*MALE*(1-INMIGRANT) +      
							female_age1_resident*(1-MALE)*AGE1*(1-INMIGRANT) + female_age2_resident*(1-MALE)*AGE2*(1-INMIGRANT),					
					a ~ dnorm(0, 100),
					c(fishing, male_age1_inmigrant, male_age2_inmigrant, male_age3_inmigrant, female_age1_inmigrant, 
						female_age2_inmigrant, female_age3_inmigrant,
						male_age1_resident, male_age2_resident, male_age3_resident, female_age1_resident, 
						female_age2_resident) ~ dnorm(0,10)
			),
			data=as.data.frame(dp),
			control=list(max_treedepth=15),
			start=list(	a=0, fishing=0, male_age1_inmigrant=0, male_age2_inmigrant=0, male_age3_inmigrant=0, female_age1_inmigrant=0, 
					female_age2_inmigrant=0, female_age3_inmigrant=0,
					male_age1_resident=0, male_age2_resident=0, male_age3_resident=0, female_age1_resident=0, 
					female_age2_resident=0),
			warmup=5e2, iter=15e3, chains=1, cores=1
		)	
	plot( precis(mp4, depth=2, prob=0.95) )
	save(mp4, dp, des, file=paste0(outfile.base,'190910_participation_differences_mp4.rda'))
	#	posterior check
	sims 	<- sim(mp4, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dp, dpp)		
	dpp[, mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
	#	0.25 
	dpp[, COMM_TYPE_F2:= as.character(factor(COMM_TYPE_F, levels=c(0,1), labels=c('inland','fisherfolk')))]
	dpp[, INMIGRANT2:= as.character(factor(INMIGRANT, levels=c(0,1), labels=c('resident','inmigrant')))]
	dpp[, CLASS:= ifelse(SUC<predicted_obs_l95, 'below 95%', ifelse(SUC>predicted_obs_u95, 'above 95%', 'within 95%'))]
	ggplot(dpp, aes(x=interaction(SEX, AGE_AT_MID_C)) ) +
			geom_bar(aes(fill=CLASS)) +
			facet_grid( INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_participation_differences_mp4_pp.pdf'), w=8, h=8)
	
	
	# model with gender:age:mig_status interactions and dispersion
	mp5 	<- map2stan(
			alist(
					SUC ~ dbetabinom(TRIAL, p_part, dispersion),
					logit(p_part) <- a + fishing * COMM_TYPE_F +  
							male_age1_inmigrant*MALE*AGE1*INMIGRANT + male_age2_inmigrant*MALE*AGE2*INMIGRANT + male_age3_inmigrant*MALE*INMIGRANT +      
							female_age1_inmigrant*(1-MALE)*AGE1*INMIGRANT + female_age2_inmigrant*(1-MALE)*AGE2*INMIGRANT + female_age3_inmigrant*(1-MALE)*INMIGRANT +
							male_age1_resident*MALE*AGE1*(1-INMIGRANT) + male_age2_resident*MALE*AGE2*(1-INMIGRANT) + male_age3_resident*MALE*(1-INMIGRANT) +      
							female_age1_resident*(1-MALE)*AGE1*(1-INMIGRANT) + female_age2_resident*(1-MALE)*AGE2*(1-INMIGRANT),					
					a ~ dnorm(0, 100),
					dispersion ~ dexp(1),
					c(fishing, male_age1_inmigrant, male_age2_inmigrant, male_age3_inmigrant, female_age1_inmigrant, 
							female_age2_inmigrant, female_age3_inmigrant,
							male_age1_resident, male_age2_resident, male_age3_resident, female_age1_resident, 
							female_age2_resident) ~ dnorm(0,10)
			),
			data=as.data.frame(dp),
			control=list(max_treedepth=15),
			start=list(	dispersion=1, a=0, fishing=0, male_age1_inmigrant=0, male_age2_inmigrant=0, male_age3_inmigrant=0, female_age1_inmigrant=0, 
					female_age2_inmigrant=0, female_age3_inmigrant=0,
					male_age1_resident=0, male_age2_resident=0, male_age3_resident=0, female_age1_resident=0, 
					female_age2_resident=0),
			warmup=5e2, iter=15e3, chains=1, cores=1
		)	
	plot( precis(mp5, depth=2, prob=0.95, pars=c( "a","fishing","male_age1_inmigrant","male_age2_inmigrant","male_age3_inmigrant","female_age1_inmigrant","female_age2_inmigrant","female_age3_inmigrant","male_age1_resident","male_age2_resident","male_age3_resident","female_age1_resident","female_age2_resident" )) )
	save(mp5, dp, des, file=paste0(outfile.base,'190910_participation_differences_mp5.rda'))
	#	posterior check
	sims 	<- sim(mp5, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dp, dpp)		
	dpp[, mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
	#	0.034 
	dpp[, COMM_TYPE_F2:= as.character(factor(COMM_TYPE_F, levels=c(0,1), labels=c('inland','fisherfolk')))]
	dpp[, INMIGRANT2:= as.character(factor(INMIGRANT, levels=c(0,1), labels=c('resident','inmigrant')))]
	dpp[, CLASS:= ifelse(SUC<predicted_obs_l95, 'below 95%', ifelse(SUC>predicted_obs_u95, 'above 95%', 'within 95%'))]
	ggplot(dpp, aes(x=interaction(SEX, AGE_AT_MID_C)) ) +
			geom_bar(aes(fill=CLASS)) +
			facet_grid( INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_participation_differences_mp5_pp.pdf'), w=8, h=8)
	
	
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95, colour=interaction(SEX,AGE_AT_MID_C)), position=position_dodge(width=0.9)) +
			geom_point(aes(y=predicted_obs_median, colour=interaction(SEX,AGE_AT_MID_C)), shape=18, position=position_dodge(width=0.9)) +			
			geom_point(aes(y=SUC, colour=interaction(SEX,AGE_AT_MID_C)), position=position_dodge(width=0.9)) +			
			facet_grid(INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_participation_differences_mp5_pred.pdf'), w=30, h=15)
	
	
}

Rakai190327.participitation.differences.betabinomialmodel7<- function()
{
  require(data.table)
  require(rstan)
  require(extraDistr)
  require(bayesplot)
  
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
  
  #	vanilla community IDs
  tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
  des		<- merge(des, tmp, by='COMM_NUM_A')
  
  #	aggregate participation
  dp		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), 
             by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6',
                  'INMIGRANT','MALE','COMM_TYPE_F')]
  dp[, CATEGORY:= paste0(COMM_TYPE_F,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
  
  #	run STAN 
  tmp			<- as.list(subset(dp,select=c('TRIAL','SUC',
                                      'MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6',
                                      'INMIGRANT','COMM_TYPE_F')))
  tmp$N		<- nrow(dp)
  fit.par 	<- stan(	file = infile.participation.stan.model, 
                    data = tmp, 
                    iter = 10e3,
                    warmup = 5e2,
                    cores = 1,
                    chains = 1,
                    control=list(max_treedepth=15),
                    init = list(list(a=0, dispersion=1, fishing=0, inmigrant=0, male=0,age1=0, age2=0,
                                     age3=0, age4=0, age5=0, age6=0, 
                                     age1_female=0, age1_male=0, age4_female=0,
                                     age7_male=0,
                                     male_inland=0, female_inland=0, 
                                     inland_inmig=0,  inland_resid=0,
                                     male_inmig=0, female_inmig=0,female_resid=0)))
  

  # assess convergence
  fit.pars	<- c('a','dispersion','fishing', 'inmigrant', 'male', 'age1', 'age2',
                'age3', 'age4', 'age5', 'age6',
                'age1_female', 'age1_male', 'age4_female', 
                'age7_male',
                'male_inland', 'female_inland', 
                'inland_inmig', 'inland_resid',
                'male_inmig','female_inmig','female_resid')
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
  dps<-dp
  fit.e		<- extract(fit.par)
  set.seed(42)
  tmp			<- sample(length(fit.e$a), nprior)
  dps			<- dps[,	
               {
                 z<- with(fit.e, a + fishing*COMM_TYPE_F + inmigrant*INMIGRANT + male*MALE + age1*AGE1 + age2*AGE2 +
                            age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6 +
                            age1_female * AGE1 * (1-MALE) + age1_male * AGE1 * MALE + age4_female * AGE4 * (1-MALE) +
                            age7_male * (1-AGE1-AGE2-AGE3-AGE4-AGE5-AGE6) * MALE +
                            male_inland * MALE * (1-COMM_TYPE_F) + female_inland * (1-MALE) * (1-COMM_TYPE_F) +
                            inland_inmig * (1-COMM_TYPE_F) * INMIGRANT + 
                            inland_resid * (1-COMM_TYPE_F) * (1-INMIGRANT)+
                            male_inmig * MALE * INMIGRANT + female_inmig * (1-MALE) * INMIGRANT + female_resid * (1-MALE) * (1-INMIGRANT) ,
                 )
                 list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
               },	
               by=c('CATEGORY')]
  dps[, P:= exp(ETA)/(1+exp(ETA))]
  
  require(bde)
  tmp	<- dps[, {
    bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
    list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
  }, by='CATEGORY']
  
  dps	<- merge(dps, tmp, by=c('CATEGORY','SAMPLE'))
  set(dps, NULL, c('ETA'), NULL)		
  
  save(dp, dps, fit.par, file=paste0(outfile.base,'190327_participation_modelage6_samples.rda'))
}

Rakai190910.sequencing.modeldev<- function()
{
	require(data.table)
	require(rstan)
	require(rethinking)
	require(bayesplot)
	
	opt.exclude.onART.from.denominator	<- 1
	indir <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data <- file.path(indir,"190327_sampling_by_gender_age.rda")	
	outfile.base <- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des <- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
	
	#	binarize age, sex
	des[, AGE1:= as.integer(AGE_AT_MID_C=='15-24')]
	des[, AGE2:= as.integer(AGE_AT_MID_C=='25-34')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	
	#	vanilla community IDs
	tmp <- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des <- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate sequencing
	des <- subset(des, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
		stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
		des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
	}
	stopifnot( all(unique(sort(des$COMM_NUM_B))==1:36) )
	ds <-des[,list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))),
			by=c('COMM_NUM_A','COMM_NUM_B','AGE_AT_MID_C','SEX','AGE1','AGE2',
					'INMIGRANT','MALE','COMM_TYPE_F')]
	ds[, CATEGORY:= paste0(COMM_TYPE_F,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	# basic model, no interactions
	ms1 	<- map2stan(
			alist(
					SUC ~ dbinom(TRIAL, p_part),
					logit(p_part) <- a + fishing * COMM_TYPE_F + inmigrant * INMIGRANT + male * MALE +      
							age1 * AGE1 + age2 * AGE2,					
					a ~ dnorm(0, 100),
					c(fishing, inmigrant, male, age1, age2) ~ dnorm(0,10)
			),
			data=as.data.frame(ds),
			control=list(max_treedepth=15),
			start=list(	a=0, fishing=0, inmigrant=0, male=0, age1=0, age2=0),
			warmup=5e2, iter=15e3, chains=1, cores=1
		)
	plot( precis(ms1, depth=2, prob=0.95) )
	save(ms1, ds, des, file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1.rda'))
	#	posterior check
	sims 	<- sim(ms1, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(ds, dpp)		
	dpp[, mean(SUC<predicted_obs_l95 | SUC>predicted_obs_u95)]
	# 0.032
	dpp[, COMM_TYPE_F2:= as.character(factor(COMM_TYPE_F, levels=c(0,1), labels=c('inland','fisherfolk')))]
	dpp[, INMIGRANT2:= as.character(factor(INMIGRANT, levels=c(0,1), labels=c('resident','inmigrant')))]
	dpp[, CLASS:= ifelse(SUC<predicted_obs_l95, 'below 95%', ifelse(SUC>predicted_obs_u95, 'above 95%', 'within 95%'))]
	ggplot(dpp, aes(x=interaction(SEX, AGE_AT_MID_C)) ) +
			geom_bar(aes(fill=CLASS)) +
			facet_grid( INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1_pp.pdf'), w=8, h=8)
	#
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95, colour=interaction(SEX,AGE_AT_MID_C)), position=position_dodge(width=0.9)) +
			geom_point(aes(y=predicted_obs_median, colour=interaction(SEX,AGE_AT_MID_C)), shape=18, position=position_dodge(width=0.9)) +			
			geom_point(aes(y=SUC, colour=interaction(SEX,AGE_AT_MID_C)), position=position_dodge(width=0.9)) +			
			facet_grid(INMIGRANT2~COMM_TYPE_F2, scales='free', space='free') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1_pred.pdf'), w=30, h=15)
	
	
}


Rakai190910.sequencing.age3.ms1<- function()
{
  require(data.table)
  require(rstan)
  require(bayesplot)
  
  opt.exclude.onART.from.denominator	<- 0
  indir <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
  infile.data <- file.path(indir,"190327_sampling_by_gender_age.rda")
  indir.stan <- '~/git/phyloscanner/phyloflows/inst/misc'
  infile.sequencing.stan.model <- file.path(indir.stan,"190910_glm_sequencing_age3_ms1.stan")
  outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
  
  #	set up variables for STAN
  load(infile.data)
  des <- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
  
  # remove individuals in the communities where no sequences were obtained successfully
  tmp <- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
  des <- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')	
  
  # set no INMIGRANT to 0
  set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
  
  # age group
  des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
  
  #	binarize age, sex
  des[, AGE1:= as.integer(AGE_AT_MID_C=='15-24')]
  des[, AGE2:= as.integer(AGE_AT_MID_C=='25-34')]
  des[, MALE:= as.integer(SEX=='M')]
  
  #	binarize community type
  des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
  
  #	vanilla community IDs
  tmp <- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
  des <- merge(des, tmp, by='COMM_NUM_A')
  
  #	aggregate sequencing
  des <- subset(des, HIV_1517==1)
  if(opt.exclude.onART.from.denominator)
  {
    stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
    des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
  }
  stopifnot( all(unique(sort(des$COMM_NUM_B))==1:36) )
  ds <-des[,	list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))),
          		by=c('COMM_NUM_A','COMM_NUM_B','AGE_AT_MID_C','SEX','AGE1','AGE2','INMIGRANT','MALE','COMM_TYPE_F')]
  ds[, CATEGORY:= paste0(COMM_TYPE_F,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
  
  #	run STAN 
  tmp <- as.list(subset(ds,select=c('AGE_AT_MID_C', 'SEX', 'AGE1', 'AGE2', 'INMIGRANT', 'MALE', 'COMM_TYPE_F', 'TRIAL', 'SUC')))
  tmp$N <- nrow(ds)
  fit.seq <- stan(file = infile.sequencing.stan.model, 
                   data = tmp, 
                   iter = 15e3,
                   warmup = 5e2,
                   cores = 1,
                   chains = 1,
                   init = list(list(a=0, fishing=0, inmigrant=0, male=0,age1=0, age2=0)))

  # assess convergence
  fit.pars <- c('a', 'fishing', 'inmigrant', 'male','age1', 'age2')
  any(rhat(fit.seq, pars=fit.pars)>1.02)
  any(neff_ratio(fit.seq, pars=fit.pars) * 9.5e3 < 500)
  po <- as.matrix(fit.seq) 
  po <- po[, colnames(po)[!grepl('p_suc_logit|lp__',colnames(po))]]	
  p <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + geom_vline(xintercept = 0)
  pdf(file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1_marginalposteriors.pdf'), w=7, h=12)
  p
  dev.off()
  p <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
  pdf(file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1_marginaltraces.pdf'), w=7, h=100)
  p
  dev.off()
  
  # summary csv file
  write.csv(  	summary(fit.seq, pars= fit.pars, probs = c(0.025, 0.975))$summary,
              	file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1_summary.csv')
  				)	
  
  # extract samples for unique strata levels
  nprior <- 1e3
  dss <- ds
  set(dss, NULL, c('COMM_NUM_A','COMM_NUM_B','TRIAL','SUC'), NULL)
  dss <- unique(dss)
  fit.e <- extract(fit.seq)
  set.seed(42)
  tmp <- sample(length(fit.e$a), nprior)
  dss <- dss[,	
               {
                 z<- with(fit.e, a + fishing*COMM_TYPE_F + 
                            inmigrant*INMIGRANT + male*MALE + age1*AGE1 + age2*AGE2)
                 list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
               },	
               by=c('CATEGORY')]
  dss[, P:= exp(ETA)/(1+exp(ETA))]
  #	fit kernel density to prior samples with Gourierous and Monfort, 2006 macrobetakernel bounded density estimator, and then estimate log density
  require(bde)
  tmp <- dss[, {
    				bdest <- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
    				list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
  				}, by='CATEGORY']
  dss <- merge(dss, tmp, by=c('CATEGORY','SAMPLE'))
  set(dss, NULL, c('ETA'), NULL)		
  save(ds, dss, fit.seq, file=paste0(outfile.base,'190910_sequencing_excART',opt.exclude.onART.from.denominator,'_age3_ms1_samples.rda'))
}

Rakai190327.sequencing.differences.binomialmodel7<- function()
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
  ds<-des[,list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))),
          by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6',
               'INMIGRANT','MALE','COMM_TYPE_F')]
  ds[, CATEGORY:= paste0(COMM_TYPE_F,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
  
  #	run STAN 
  tmp			<- as.list(subset(ds,select=c('AGE_AT_MID_C', 'SEX', 'AGE1', 'AGE2', 'AGE3', 'AGE4', 'AGE5','AGE6',
                                      'INMIGRANT', 'MALE', 'COMM_TYPE_F', 'TRIAL', 'SUC')))
  tmp$N		<- nrow(ds)
  fit.seq 	<- stan(file = infile.sequencing.stan.model, 
                   data = tmp, 
                   iter = 10e3,
                   warmup = 5e2,
                   cores = 1,
                   chains = 1,
                   init = list(list(a=0, fishing=0, inmigrant=0, male=0,age1=0, age2=0,
                                    age3=0, age4=0, age5=0, age6=0)))
  
  # assess convergence
  fit.pars	<- c('a', 'fishing', 'inmigrant', 'male','age1', 'age2','age3', 'age4', 'age5', 'age6')
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
  dss   <- ds
  fit.e		<- extract(fit.seq)
  set.seed(42)
  tmp			<- sample(length(fit.e$a), nprior)
  dss			<- dss[,	
               {
                 z<- with(fit.e, a + fishing*COMM_TYPE_F + 
                            inmigrant*INMIGRANT + male*MALE + age1*AGE1 + age2*AGE2 +
                            age3 * AGE3 + age4 * AGE4 + age5 * AGE5 + age6 * AGE6)
                 list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
               },	
               by=c('CATEGORY')]
  dss[, P:= exp(ETA)/(1+exp(ETA))]
  
  require(bde)
  tmp	<- dss[, {
    bdest<- bde(P, dataPointsCache=sort(P), b=0.001, estimator='betakernel', lower.limit=0, upper.limit=1, options=list(modified=FALSE, normalization='densitywise', mbc='none', c=0.5))
    list(SAMPLE=SAMPLE, LP=log(density(bdest, P)))
  }, by='CATEGORY']
  
  dss	<- merge(dss, tmp, by=c('CATEGORY','SAMPLE'))
  set(dss, NULL, c('ETA'), NULL)		
  
  save(ds, dss, fit.seq, file=paste0(outfile.base,'190327_sequencing_modelage6_samples.rda'))
}


Rakai190910.analysispipeline.age3model<- function(infile.inference=NULL, infile.participation.prior.samples=NULL, infile.sequencing.prior.samples=NULL, opt=NULL)
{
  require(data.table)	
  require(phyloflows)
  
  #
  #	input args
  #
  if(is.null(opt))
  {
    opt									<- list()
    opt$adjust.sequencing.bias			<- 0
    opt$adjust.participation.bias		<- 0
    opt$migration.def.code				<- '24'
    opt$set.missing.migloc.to.inland	<- 0
    opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
  }
  if(is.null(infile.inference))
  {
    infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
  }
  if(is.null(infile.participation.prior.samples))
  {
	  infile.participation.prior.samples	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190910_participation_age3_mp5_samples.rda"	
  }
  if(is.null(infile.sequencing.prior.samples) & opt$adjust.sequencing.bias<2)
  {
    infile.sequencing.prior.samples		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190910_sequencing_excART1_age3_ms1_samples.rda"
  }
  if(is.null(infile.sequencing.prior.samples) & opt$adjust.sequencing.bias==2)
  {
	infile.sequencing.prior.samples		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190910_sequencing_excART0_age3_ms1_samples.rda"
  }
  if(opt$adjust.sequencing.bias==2 && !grepl('excART0',infile.sequencing.prior.samples))
  {
	stop('unexpected infile.sequencing.prior.samples')
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
  
  rtr	<- subset(rtr, select=c('PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
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
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fish')
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

  # define TR_COMM_TYPE_F, REC_COMM_TYPE_F (i: inland; f: fishing) 
  rtr[,TR_COMM_TYPE_F:=substr(TR_COMM_TYPE,1,1)]
  rtr[substr(TR_COMM_TYPE,1,1)!='f',TR_COMM_TYPE_F:='i']
  rtr[,REC_COMM_TYPE_F:=substr(REC_COMM_TYPE,1,1)]
  rtr[substr(REC_COMM_TYPE,1,1)!='f',REC_COMM_TYPE_F:='i']
  
  # define TR_COMM_TYPE_F_MIG (i: inland; f: fishing; e: external) 
  rtr[,TR_COMM_TYPE_F_MIG:=substr(TR_COMM_NUM_A_MIG,1,1)]
  rtr[substr(TR_COMM_NUM_A_MIG,1,1)=='a' | substr(TR_COMM_NUM_A_MIG,1,1)=='i'|
        substr(TR_COMM_NUM_A_MIG,1,1)=='t',TR_COMM_TYPE_F_MIG:='i']
  
  #	build category to match with sampling data tables 
  rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
  rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_TYPE_F,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
  #	build transmission flow category 
  rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
  rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_TYPE_F_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
  
  # make all combinations of variables
  dac <- expand.grid( 	COMM_TYPE_F= sort(unique(c(rtr$REC_COMM_TYPE_F, rtr$TR_COMM_TYPE_F))),
						SEX=  sort(unique(c(rtr$REC_SEX, rtr$TR_SEX))),
						AGE_AT_MID_C= sort(unique(c(rtr$REC_AGE_AT_MID_C, rtr$TR_AGE_AT_MID_C))),
    					INMIGRANT= sort(unique(c(rtr$REC_INMIGRANT, rtr$TR_INMIGRANT)))
    					)
  dac <- as.data.table(dac)  				
  dac[, CATEGORY:= paste0(COMM_TYPE_F, ':', SEX, ':', AGE_AT_MID_C, ':', INMIGRANT)]  					
  dac <- as.data.table(expand.grid(TR_CATEGORY= dac$CATEGORY, REC_CATEGORY= dac$CATEGORY))
  # ignore Male-Male and Female-Female combinations
  dac <- subset(dac, !(grepl('F',TR_CATEGORY)&grepl('F',REC_CATEGORY)) &
  					 !(grepl('M',TR_CATEGORY)&grepl('M',REC_CATEGORY))  
					  )
  # add transmission categories
  dac[, REC_TRM_CATEGORY:= REC_CATEGORY]
  dac[, TR_TRM_CATEGORY:= TR_CATEGORY]  
  # add inmigrants from external communities
  tmp <- dac[grepl('1$',TR_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^[f|i]','e',TR_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # add inmigrants sampled in inland and migrated from fishing
  tmp <- dac[grepl('1$',TR_CATEGORY) & grepl('^i',TR_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^i','f',TR_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # add inmigrants sampled in fishing and migrated from inland
  tmp <- dac[grepl('1$',TR_CATEGORY) & grepl('^f',TR_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^f','i',TR_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # remove duplicated rows 
  # TR_SAMPLING_CATEGORY f: F: 15-24: 1 and i: F: 15-24: 1 are all set to e: F: 15-24: 1
  dac <- unique(dac)
  setnames(dac, c('TR_CATEGORY', 'REC_CATEGORY'), c('TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'))
  
  #	calculate observed number of transmissions
  #
  dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
  dac[, DUMMY:= 1]
  dobs <- merge(dac, dobs, by=c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'), all=TRUE)
  stopifnot( dobs[, !any(is.na(DUMMY))] )
  set(dobs, NULL, 'DUMMY', NULL)
  set(dobs, dobs[, which(is.na(TRM_OBS))], 'TRM_OBS', 0L)
  
  #	make PAIR_ID
  setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY,TR_SAMPLING_CATEGORY,REC_SAMPLING_CATEGORY)	
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
  # multiply participation and sequencing probabilities
  dprior[, P:= P.x*P.y]		
  # add log posterior densities
  dprior[, LP:= LP.x+LP.y]	
  set(dprior, NULL, c('P.x','LP.x','P.y','LP.y'), NULL)
  setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')
  # consolidate names
  set(dprior, NULL, 'SAMPLING_CATEGORY', dprior[, gsub('^2','e',gsub('^1','f',gsub('^0','i',SAMPLING_CATEGORY)))] )		
  
  #	plot prior densities
  if(0)
  {
	  ggplot(dprior, aes(x=SAMPLING_CATEGORY, y=P)) + 
			  geom_boxplot(bins=50) +
			  scale_y_continuous(labels=scales:::percent, lim=c(0,1), expand=c(0,0)) +
			  labs(x='sampling category', y='sampling probability') + 
			  theme_bw() + 
			  coord_flip()
	  ggsave( file=paste0(outfile.base,"saprior190910_samplingprob.pdf"), w=7, h=7, limitsize=FALSE)    		
  }
  
  #	run MCMC
  mcmc.file	<- paste0(outfile.base,"samcmc190327_nsweep1e5_opt",paste0(opt, collapse=''))
  control		<- list(	seed=42, 
		                    mcmc.n=48 * 1e5,							
		                    verbose=0, 
		                    outfile=mcmc.file,
		                    sweep_group=1e4L)
  source.attribution.mcmc(dobs, dprior, control=control)
  
  #	run diagnostics
  #mcmc.file	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_SAMCMCv190327_mcmc.rda"
  #mcmc.file	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_181006_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_samcmc190327_nsweep1e5_opt112401.rda'  
  mcmc.file.full <- paste0(mcmc.file,'_sweepgrp',1:(1e5%/%1e4),'.rda')
  control		<- list(	burnin.p=0.05, 
                    		regex_pars='*', 
                    		credibility.interval=0.95, 
                    		pdf.plot.all.parameters=FALSE, 
                    		pdf.plot.n.worst.case.parameters=10, 
                    		pdf.height.per.par=1.2, 
                    		outfile.base=mcmc.file)
  source.attribution.mcmc.diagnostics(mcmc.file=mcmc.file.full, control=control)
  
  #	aggregate MCMC output to fish<->inland
  aggregate.file	<- paste0(mcmc.file,'_aggregatedFishInland.csv')
  daggregateTo	<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
  daggregateTo[, TR_TARGETCAT:= gsub('^(.)\\:(.)\\:(.+)\\:(.)$','\\1',TR_TRM_CATEGORY)]
  daggregateTo[, REC_TARGETCAT:= gsub('^(.)\\:(.)\\:(.+)\\:(.)$','\\1',REC_TRM_CATEGORY)]
  set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)	
  control			<- list(	burnin.p=0.05, 
                     thin=NA_integer_, 
                     regex_pars='*', 
                     outfile=aggregate.file)
  source.attribution.mcmc.aggregateToTarget(mcmc.file=mcmc.file.full, 
                                            daggregateTo=daggregateTo, 
                                            control=control)
  
  #	calculate flows sources WAIFM flow_ratio overall		
  control			<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                     			flowratios= list( c('i/f', 'i f', 'f i')),
                     			outfile=gsub('\\.csv','_flowsetc.csv',aggregate.file))
  source.attribution.mcmc.getKeyQuantities(infile=aggregate.file, control=control)
  
  #	aggregate MCMC output to fish<->inland by gender
  aggregate.file	<- paste0(mcmc.file,'_aggregatedFishInlandByGender.csv')
  daggregateTo	<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
  daggregateTo[, TR_TARGETCAT:= paste0(
    gsub('^(.)\\:(.)\\:(.+)\\:(.)$','\\1',TR_TRM_CATEGORY),
    ':',
    gsub('^(.)\\:(.)\\:(.+)\\:(.)$','\\2',TR_TRM_CATEGORY))]
  daggregateTo[, REC_TARGETCAT:= paste0(
    gsub('^(.)\\:(.)\\:(.+)\\:(.)$','\\1',REC_TRM_CATEGORY),
    ':',
    gsub('^(.)\\:(.)\\:(.+)\\:(.)$','\\2',REC_TRM_CATEGORY))]
  set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)	
  control		<- list(	burnin.p=0.05, 
                    thin=NA_integer_, 
                    regex_pars='*', 
                    outfile=aggregate.file)
  source.attribution.mcmc.aggregateToTarget(mcmc.file=mcmc.file.full,
                                            daggregateTo=daggregateTo, 
                                            control=control)
  
  #	calculate flows sources WAIFM flow_ratio by gender		
  control		<- list(	quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                    flowratios= list( c('i:M/f:M', 'i:M f:F', 'f:M i:F'), 
                                      c('i:F/f:F', 'i:F f:M', 'f:F i:M')),
                    outfile=gsub('\\.csv','_flowsetc.csv',aggregate.file))
  source.attribution.mcmc.getKeyQuantities(infile=aggregate.file, control=control)				
}


Rakai190327.analysispipeline.age6model<- function(infile.inference=NULL, infile.participation.prior.samples=NULL, infile.sequencing.prior.samples=NULL, opt=NULL)
{
  require(data.table)	
  require(phyloflows)
  
  #
  #	input args
  #
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
    infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
  }
  if(is.null(infile.participation.prior.samples))
  {
    infile.participation.prior.samples	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_participation_modelage6_samples"	
  }
  if(is.null(infile.sequencing.prior.samples))
  {
    infile.sequencing.prior.samples		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_sequencing_modelage6_samples.rda"
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
  
  rtr	<- subset(rtr, select=c('PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
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
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fish')
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
  rtr[, TR_AGE_AT_MID_C:= cut(TR_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
  rtr[, REC_AGE_AT_MID_C:= cut(REC_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
  stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
  stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
  
  # unique(substr(rtr$TR_COMM_TYPE,1,1))
  # unique(substr(rtr$REC_COMM_TYPE,1,1))
  rtr[,TR_COMM_TYPE_F:=substr(TR_COMM_TYPE,1,1)]
  rtr[substr(TR_COMM_TYPE,1,1)!='f',TR_COMM_TYPE_F:='i']
  rtr[,REC_COMM_TYPE_F:=substr(REC_COMM_TYPE,1,1)]
  rtr[substr(REC_COMM_TYPE,1,1)!='f',REC_COMM_TYPE_F:='i']
  
  # unique(substr(rtr$TR_COMM_NUM_A_MIG,1,1))
  rtr[,TR_COMM_TYPE_F_MIG:=substr(TR_COMM_NUM_A_MIG,1,1)]
  rtr[substr(TR_COMM_NUM_A_MIG,1,1)=='a' | substr(TR_COMM_NUM_A_MIG,1,1)=='i'|
        substr(TR_COMM_NUM_A_MIG,1,1)=='t',TR_COMM_TYPE_F_MIG:='i']
  unique(rtr$TR_COMM_TYPE_F_MIG)
  
  #	build category to match with sampling data tables 
  rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
  rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_TYPE_F,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
  #	build transmission flow category 
  rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
  rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_TYPE_F_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
  
  # make all combinations of variables
  dac <- expand.grid( 	COMM_TYPE_F= sort(unique(c(rtr$REC_COMM_TYPE_F, rtr$TR_COMM_TYPE_F))),
                       SEX=  sort(unique(c(rtr$REC_SEX, rtr$TR_SEX))),
                       AGE_AT_MID_C= unique(rtr$REC_AGE_AT_MID_C),
                       INMIGRANT= sort(unique(c(rtr$REC_INMIGRANT, rtr$TR_INMIGRANT)))
  )
  dac <- as.data.table(dac)  				
  dac[, CATEGORY:= paste0(COMM_TYPE_F, ':', SEX, ':', AGE_AT_MID_C, ':', INMIGRANT)]  					
  dac <- as.data.table(expand.grid(TR_CATEGORY= dac$CATEGORY, REC_CATEGORY= dac$CATEGORY))
  dac <- subset(dac, !(grepl('F',TR_CATEGORY)&grepl('F',REC_CATEGORY)) &
                  !(grepl('M',TR_CATEGORY)&grepl('M',REC_CATEGORY))  
  )
  
  # add transmission categories
  
  # same as sampling categories
  dac[, REC_TRM_CATEGORY:= REC_CATEGORY]
  dac[, TR_TRM_CATEGORY:= TR_CATEGORY]
  
  # inmigrants from external communities
  tmp <- dac[grepl('1$',TR_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^[f|i]','e',TR_CATEGORY)]) 
  dac <- rbind(dac, tmp)
  
  # inmigrants sampled in inland and migrated from fishing
  tmp <- dac[grepl('1$',TR_CATEGORY) & grepl('^i',TR_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^i','f',TR_CATEGORY)]) 
  dac <- rbind(dac, tmp)
  
  # inmigrants sampled in fishing and migrated from inland
  tmp <- dac[grepl('1$',TR_CATEGORY) & grepl('^f',TR_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^f','i',TR_CATEGORY)]) 
  dac <- rbind(dac, tmp)
  
  # remove duplicated rows 
  # TR_SAMPLING_CATEGORY f: F: 15-24: 1 and i: F: 15-24: 1 are all set to e: F: 15-24: 1
  dac <- unique(dac)
  setnames(dac, c('TR_CATEGORY', 'REC_CATEGORY'), c('TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'))
  
  #	calculate observed number of transmissions
  #
  dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
  
  # TODO: can we make adding zeros independent of any other files please
  #	colnames --> setnames
  
  # add TRM_OBS
  # zero TRM_OBS
  tmp1 <- dac[! paste0(TR_SAMPLING_CATEGORY,'-',REC_SAMPLING_CATEGORY, '-', TR_TRM_CATEGORY,'-',REC_TRM_CATEGORY)
              %in% paste0(dobs$TR_SAMPLING_CATEGORY,'-',dobs$REC_SAMPLING_CATEGORY, '-', dobs$TR_TRM_CATEGORY,'-',dobs$REC_TRM_CATEGORY),]
  # non-zero TRM_OBS
  tmp2 <- dac[ paste0(TR_SAMPLING_CATEGORY,'-',REC_SAMPLING_CATEGORY, '-', TR_TRM_CATEGORY,'-',REC_TRM_CATEGORY)
               %in% paste0(dobs$TR_SAMPLING_CATEGORY,'-',dobs$REC_SAMPLING_CATEGORY, '-', dobs$TR_TRM_CATEGORY,'-',dobs$REC_TRM_CATEGORY),]
  
  tmp1[, TRM_OBS:=0]
  dobs <- merge(dobs, tmp2, by = c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY',
                                   'TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'))
  dobs <- rbind(tmp1, dobs)
  
  setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY,TR_SAMPLING_CATEGORY,REC_SAMPLING_CATEGORY)	
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
  mcmc.file	<- paste0(outfile.base,"samcmc190327_nsweep1e5_6age_opt",paste0(opt, collapse=''))
  control		<- list(	seed=42, 
                    mcmc.n=112 * 1e5, 
                    verbose=0, 
                    outfile=mcmc.file,
                    sweep_group=1e4L)
  source.attribution.mcmc(dobs, dprior, control=control)
 				
}
