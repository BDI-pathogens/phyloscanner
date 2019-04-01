lddirichlet_vector	<- function(x, nu){
	ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
	stopifnot(is.finite(ans))
	ans
}


source.attribution.mcmc	<- function(dobs, dprior, control=list(seed=42, mcmc.n=1e3, verbose=1, outfile='SAMCMCv190327.rda')){
  library(data.table)
  library(gtools)
  
  #	
  # basic checks
  #
  if(!all(dobs$TR_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
	  stop('Did not find prior samples of a sampling category for TR_SAMPLING_CATEGORY')
  if(!all(dobs$REC_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
	  stop('Did not find prior samples of a sampling category for REC_SAMPLING_CATEGORY')
  
  ptm	<- Sys.time()
  set.seed(control$seed)
  
  #	
  # set up mcmc 
  #
  mc				<- list()
  mc$time			<- NA_real_
  # construct look-up table so we know which transmission pair categories need to be updated 
  # at every MCMC iteration
  tmp				<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_SAMPLING_CATEGORY, REC_SAMPLING_CATEGORY))
  tmp				<- melt(tmp, id.vars='TRM_CAT_PAIR_ID', value.name='SAMPLING_CATEGORY', variable.name='WHO')
  mc$dl				<- unique(subset(tmp, select=c(WHO, SAMPLING_CATEGORY)))
  mc$dl[, UPDATE_ID:= seq_len(nrow(mc$dl))]
  mc$dl				<- merge(mc$dl, tmp, by=c('WHO','SAMPLING_CATEGORY'))
  setkey(mc$dl, UPDATE_ID)
  mc$dlt			<- dcast.data.table(mc$dl, TRM_CAT_PAIR_ID~WHO, value.var='UPDATE_ID')
  setnames(mc$dlt, c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY'), c('TR_UPDATE_ID','REC_UPDATE_ID'))
  #	every transmission category pair needs to be updated at least one as we sweep through the sampling categories
  # I don t think this is guaranteed, hence the check
  if( !all(dobs$TRM_CAT_PAIR_ID %in% sort(unique(mc$dl$TRM_CAT_PAIR_ID))) )
	  stop('Fatal error. Contact the package maintainer with your input data.')
  mc$nprior			<- max(dprior$SAMPLE)
  mc$sweep			<- max(mc$dl$UPDATE_ID)+1L
  mc$nsweep			<- ceiling( control$mcmc.n/mc$sweep )  
  mc$n				<- mc$nsweep*mc$sweep
  mc$pars			<- list()
  mc$pars$LAMBDA	<- matrix(NA_real_, ncol=nrow(dobs), nrow=1)		#prior for proportions
  mc$pars$XI		<- matrix(NA_real_, ncol=mc$sweep-1L, nrow=mc$nsweep+1L) # prior for sampling in transmitter categories and sampling in recipient categories, concatenated
  mc$pars$XI_LP		<- matrix(NA_real_, ncol=mc$sweep-1L, nrow=mc$nsweep+1L) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated   
  mc$pars$S			<- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$nsweep+1L) # prior probability of sampling transmission pair categories
  mc$pars$S_LP		<- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$nsweep+1L) # log prior density of sampling transmission pair categories
  mc$pars$Z			<- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$nsweep+1L) #augmented data
  mc$pars$NU		<- NA_real_															#prior for N
  mc$pars$N			<- matrix(NA_integer_, ncol=1, nrow=mc$nsweep+1L)							#total number of counts on augmented data
  mc$pars$PI		<- matrix(NA_real_, ncol=nrow(dobs), nrow=mc$nsweep+1L)	#proportions
  mc$it.info		<- data.table(	IT= seq.int(0,mc$n),
									PAR_ID= rep(NA_integer_, mc$n+1L),
									BLOCK= rep(NA_character_, mc$n+1L),
									MHRATIO= rep(NA_real_, mc$n+1L),
									ACCEPT=rep(NA_integer_, mc$n+1L), 
									LOG_LKL=rep(NA_real_, mc$n+1L),
									LOG_PRIOR=rep(NA_real_, mc$n+1L))  
  
  if(1)
  {
	  cat('\nNumber of parameters:\t', ncol(mc$pars$PI)+ncol(mc$pars$N)+ncol(mc$pars$Z)+ncol(mc$pars$S)+ncol(mc$pars$XI) )
	  cat('\nDimension of PI:\t', ncol(mc$pars$PI))
	  cat('\nSweep length:\t', mc$sweep)
	  cat('\nNumber of sweeps:\t', mc$nsweep)
	  cat('\nNumber of iterations:\t', mc$n)
	  tmp				<- mc$dl[, list(N_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
	  cat('\nNumber of transmission pair categories updated per iteration, and their frequencies:\n')
	  print(table(tmp$N_PAIRS))	  
  }
  
	
  #
  # initialise MCMC
  #	  
  mc$curr.it			<- 1L
  set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
  set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
  set(mc$it.info, mc$curr.it, 'MHRATIO', 1)
  set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)  
  #	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
  #	(https://projecteuclid.org/euclid.ba/1422556416)
  mc$pars$LAMBDA[1,]	<- 0.8/nrow(dobs)
  # prior for sampling in transmitter categories
  tmp					<- subset(dprior, SAMPLE==sample(mc$nprior,1))
  tmp					<- merge(unique(subset(mc$dl, WHO=='TR_SAMPLING_CATEGORY', c(SAMPLING_CATEGORY, UPDATE_ID))), tmp, by='SAMPLING_CATEGORY')
  setkey(tmp, UPDATE_ID)
  mc$pars$XI[1, tmp$UPDATE_ID]		<- tmp$P
  mc$pars$XI_LP[1, tmp$UPDATE_ID]	<- tmp$LP
  # prior for sampling in in recipient categories
  tmp					<- subset(dprior, SAMPLE==sample(mc$nprior,1))
  tmp					<- merge(unique(subset(mc$dl, WHO=='REC_SAMPLING_CATEGORY', c(SAMPLING_CATEGORY, UPDATE_ID))), tmp, by='SAMPLING_CATEGORY')
  setkey(tmp, UPDATE_ID)
  mc$pars$XI[1, tmp$UPDATE_ID]		<- tmp$P
  mc$pars$XI_LP[1, tmp$UPDATE_ID]	<- tmp$LP
  # prior for sampling in transmission pair categories
  tmp					<- mc$dl[,  list( 	TR_UPDATE_ID= UPDATE_ID[WHO=='TR_SAMPLING_CATEGORY'][1],
											REC_UPDATE_ID= UPDATE_ID[WHO=='REC_SAMPLING_CATEGORY'][1]), 
										by='TRM_CAT_PAIR_ID']
  setkey(tmp,TRM_CAT_PAIR_ID)
  mc$pars$S[1,]			<- mc$pars$XI[1, tmp$TR_UPDATE_ID] * mc$pars$XI[1, tmp$REC_UPDATE_ID]
  mc$pars$S_LP[1,]		<- mc$pars$XI_LP[1, tmp$TR_UPDATE_ID] + mc$pars$XI_LP[1, tmp$REC_UPDATE_ID]
  #	augmented data: proposal draw under sampling probability
  mc$pars$Z[1,]			<- dobs$TRM_OBS + rnbinom(nrow(dobs),dobs$TRM_OBS,mc$pars$S[1,])
  #	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
  mc$pars$NU			<- sum(dobs$TRM_OBS) / mean(mc$pars$S[1,])
  #	total count: that s just the sum of Z
  mc$pars$N[1,]			<- sum(mc$pars$Z[1,])
  #	proportions: draw from full conditional
  mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])
  #	store log likelihood
  tmp	<- sum( dbinom(dobs$TRM_OBS, size=mc$pars$Z[1,], prob=mc$pars$S[1,], log=TRUE) ) +
    dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
  set(mc$it.info, 1L, 'LOG_LKL', tmp)
  # 	store log prior
  tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
    lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
    sum(mc$pars$S_LP[1,])
  set(mc$it.info, 1L, 'LOG_PRIOR', tmp)

  
  # parameter value at the current step
  XI.curr	<- mc$pars$XI[1,]
  XI_LP.curr<- mc$pars$XI_LP[1,]
  S.curr	<- mc$pars$S[1,]
  S_LP.curr	<- mc$pars$S_LP[1,]
  PI.curr	<- mc$pars$PI[1,]
  N.curr	<- mc$pars$N[1,]
  Z.curr	<- mc$pars$Z[1,]

  # run mcmc
  options(warn=0)  
  for(i in 1L:mc$n){
    mc$curr.it		<- i
    # determine source-recipient combination that will be updated in this iteration
    update.count	<- (i-1L) %% mc$sweep + 1L
    update.round 	<- (i-1L) %/% mc$sweep + 1L
	# update in one go S, Z, N for the ith XI
    if(update.count<mc$sweep)
    {
	  update.info	<- subset(mc$dl, UPDATE_ID==update.count)	
	  # propose single XI
	  XI.prop		<- XI.curr
	  XI_LP.prop	<- XI_LP.curr
	  tmp			<- subset(dprior, SAMPLING_CATEGORY==update.info$SAMPLING_CATEGORY[1] & SAMPLE==sample(mc$nprior,1))
	  XI.prop[ update.info$UPDATE_ID[1] ]	<- tmp$P	  
	  XI_LP.prop[ update.info$UPDATE_ID[1] ]<- tmp$LP
	  # propose all S that involve the one XI from above 
	  S.prop						<- S.curr
	  S_LP.prop						<- S_LP.curr
	  tmp							<- subset(mc$dlt, TR_UPDATE_ID==update.count | REC_UPDATE_ID==update.count)
	  setkey(tmp,TRM_CAT_PAIR_ID)
	  S.prop[tmp$TRM_CAT_PAIR_ID]	<- XI.prop[tmp$TR_UPDATE_ID] * XI.prop[tmp$REC_UPDATE_ID]
	  S_LP.prop[tmp$TRM_CAT_PAIR_ID]<- XI_LP.prop[tmp$TR_UPDATE_ID] + XI_LP.prop[tmp$REC_UPDATE_ID]
	  # propose all Z that involve a new S
	  Z.prop						<- Z.curr
	  trm.obs 						<- merge(tmp, dobs, by='TRM_CAT_PAIR_ID')[, TRM_OBS]	   
	  Z.prop[tmp$TRM_CAT_PAIR_ID]	<- trm.obs + rnbinom(length(trm.obs), trm.obs, S.prop[tmp$TRM_CAT_PAIR_ID])
	  # propose total of Z
	  N.prop						<- sum(Z.prop)
	  #	calculate MH ratio
	  log.prop.ratio	<- sum(dnbinom(Z.curr[tmp$TRM_CAT_PAIR_ID]-trm.obs, size=trm.obs, prob=S.curr[tmp$TRM_CAT_PAIR_ID], log=TRUE)) - 
						   sum(dnbinom(Z.prop[tmp$TRM_CAT_PAIR_ID]-trm.obs, size=trm.obs, prob=S.prop[tmp$TRM_CAT_PAIR_ID], log=TRUE))		
	  log.fc			<- sum(dbinom(trm.obs, size=Z.curr[tmp$TRM_CAT_PAIR_ID], prob=S.curr[tmp$TRM_CAT_PAIR_ID], log=TRUE)) +
						   dmultinom(Z.curr[tmp$TRM_CAT_PAIR_ID], prob=PI.curr[tmp$TRM_CAT_PAIR_ID], log=TRUE) +
						   dpois(N.curr, lambda=mc$pars$NU, log=TRUE)
	  log.fc.prop		<- sum(dbinom(trm.obs, size=Z.prop[tmp$TRM_CAT_PAIR_ID], prob=S.prop[tmp$TRM_CAT_PAIR_ID], log=TRUE)) +
			  			   dmultinom(Z.prop[tmp$TRM_CAT_PAIR_ID], prob=PI.curr[tmp$TRM_CAT_PAIR_ID], log=TRUE) +
			  			   dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
	  log.mh.ratio		<- log.fc.prop - log.fc + log.prop.ratio
	  mh.ratio			<- min(1,exp(log.mh.ratio))	
	  #	update
	  mc$curr.it				<- mc$curr.it+1L
	  set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
	  set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
	  set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
	  set(mc$it.info, mc$curr.it, 'ACCEPT', as.integer(runif(1) < mh.ratio))
	  if(control$verbose & mc$it.info[mc$curr.it, ACCEPT])
	  {
		  cat('\nit ',mc$curr.it,' ACCEPT S-Z-N block ',update.count)
	  }
	  if(mc$it.info[mc$curr.it, ACCEPT])
	  {
		  XI.curr	<- XI.prop
		  XI_LP.curr<- XI_LP.prop
		  S.curr	<- S.prop
		  S_LP.curr	<- S_LP.prop
		  Z.curr	<- Z.prop
		  N.curr	<- N.prop				  				
	  }	  
    }
    # update PI
    if(update.count==mc$sweep)
    {
      #	propose
      PI.prop		<- rdirichlet(1L, Z.curr + mc$pars$LAMBDA[1,])
      #	this is the full conditional of PI given S, N, Z
      #	always accept
      #	update
	  mc$curr.it	<- mc$curr.it+1L
	  set(mc$it.info, mc$curr.it, 'BLOCK', 'PI')
	  set(mc$it.info, mc$curr.it, 'PAR_ID', NA_integer_)
	  set(mc$it.info, mc$curr.it, 'MHRATIO', 1L)
	  set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)		
      PI.curr		<- PI.prop

	  # at the end of sweep, record current parameters
      mc$pars$XI[update.round+1L,]		<- XI.curr
	  mc$pars$XI_LP[update.round+1L,]	<- XI_LP.curr
      mc$pars$S[update.round+1L,]		<- S.curr
	  mc$pars$S_LP[update.round+1L,]	<- S_LP.curr
      mc$pars$Z[update.round+1L,]		<- Z.curr      
      mc$pars$N[update.round+1L,]		<- N.curr
      mc$pars$PI[update.round+1L,]		<- PI.curr
    }
	# 	record log likelihood
	tmp	<- sum(dbinom(dobs$TRM_OBS, size=Z.curr, prob=S.curr, log=TRUE) ) +
			dmultinom(Z.curr, size=N.curr, prob=PI.curr, log=TRUE)
	set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
	# 	record log prior
	tmp	<- dpois(N.curr, lambda=mc$pars$NU, log=TRUE) +
			lddirichlet_vector(PI.curr, nu=mc$pars$LAMBDA[1,]) +
			sum(S_LP.curr)
	set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)
	#
	if(update.count==mc$sweep & update.round %% 100 == 0)
		cat('\nSweeps done:\t',update.round)	
  }
  
  mc$time	<- Sys.time()-ptm
  save(mc,	file=control$outfile)
  NULL
}

source.attribution.mcmc.diagnostics	<- function(mcmc.file, control=list(burnin.p=0.2, regex_pars='*', pdf.height.per.par=1.2, outfile.base='SAMCMCv190327')){
  library(coda)
  library(data.table)
  library(bayesplot)
  diagnostic<-list()

  
  mcmc.file	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_SAMCMCv190327_mcmc.rda"
  control=list(burnin.p=0.2, regex_pars='PI', pdf.height.per.par=1.2, outfile.base=gsub('\\.rda','',mcmc.file))
  
  load(mcmc.file)
  burnin.n	<- floor(control$burnin.p*nrow(mc$pars$S))
  pars		<- matrix(NA,nrow=nrow(mc$pars$S),ncol=0)	
  if(grepl(regex_pars,'S'))
  {
	  tmp	<- mc$pars$S
	  colnames(tmp)	<- paste0('S-',1:ncol(tmp))	   
	  pars	<- cbind(pars, tmp)  
  }	
  if(grepl(regex_pars,'Z'))
  {
	  tmp	<- mc$pars$Z
	  colnames(tmp)	<- paste0('Z-',1:ncol(tmp))	   
	  pars	<- cbind(pars, tmp)  
  }
  if(grepl(regex_pars,'N'))
  {
	  tmp	<- mc$pars$N
	  colnames(tmp)	<- paste0('N-',1:ncol(tmp))	   
	  pars	<- cbind(pars, tmp)  
  }
 if(grepl(regex_pars,'PI'))
 {
	 tmp	<- mc$pars$PI
	 colnames(tmp)	<- paste0('PI-',1:ncol(tmp))	   
	 pars	<- cbind(pars, tmp)  
 }

  #	traces for parameters
  p		<- mcmc_trace(pars, pars=colnames(pars), facet_args = list(ncol = 1), n_warmup=burnin.n)
  pdf(file=paste0(control$outfile.base,'_marginaltraces.pdf'), w=7, h=control$pdf.height.per.par*ncol(pars))
  p
  dev.off()
  
  
  #	acceptance rate per MCMC update ID
  da	<- subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, list(ACC_RATE=mean(ACCEPT)), by='PAR_ID']
  setnames(da, 'PAR_ID', 'UPDATE_ID')
  tmp	<- mc$dl[, list(N_TRM_CAT_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
  da	<- merge(da, tmp, by='UPDATE_ID')
  ggplot(da, aes(x=N_TRM_CAT_PAIRS, y=ACC_RATE)) + 
		  geom_point() +
		  theme_bw() +
		  scale_y_continuous(label=scales:::percent) +
		  labs(	x='\nNumber of transmission pair categories updated per sampling category',
				y='Acceptance rate\n')
  cat('\nAverage acceptance rate= ',subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, round(mean(ACCEPT), d=3)])
  cat('\nUpdate IDs with lowest acceptance rates')
  print( da[order(ACC_RATE)[1:10],] )
  
  # remove burn-in
  cat('\nRemoving burnin in set to ', 100*control$burnin.p,'% of chain, total iterations=',burnin.n)
  it.rm.burnin	<- seq.int(burnin.n,nrow(mc$pars$S))
  pars	<- pars[it.rm.burnin,,drop=FALSE]
  
  
  # effective sampling sizes 
  tmp	<- mcmc(pars)
  ans	<- data.table(VAR= colnames(pars), NEFF=as.numeric(effectiveSize(tmp)))
  ggplot(ans, aes(x=NEFF, y=VAR)) + 
		  geom_point() + 	
		  theme_bw() +

  eff.size.PI<-matrix(NA_real_,nrow=1,ncol=ncol(mc$pars$Z))
  for(i in 1:ncol(mc$pars$Z)){
    eff.size.PI[i]<-effectiveSize(PI[,i])
  }

  sort.effsize.PI<-sort(eff.size.PI,index.return=TRUE)
  if (length(sort.effsize.PI$x)<6){
    worst.effsize.PI.index<-sort.effsize.PI$ix
    worst.effsize.PI<-sort.effsize.PI$x
    print(paste0("The worst effective sample sizes for PI ", worst.effsize.PI.index, " are ", worst.effsize.PI," out of ", length(PI[,1])," samples."))
  }else{
    worst.effsize.PI.index<-sort.effsize.PI$ix[1:6]
    worst.effsize.PI<-sort.effsize.PI$x[1:6]
    print(paste0("The worst effective sample sizes for PI ", worst.effsize.PI.index, " are ", worst.effsize.PI," out of ", length(PI[,1])," samples."))
  }

  diagnostic$PI$eff.size.index<-worst.effsize.PI.index
  diagnostic$PI$eff.size<-worst.effsize.PI

  # traceplot for PI
  png(paste0('traceplot_PI_',mc$method,'.png'), width = 1400, height = 800)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(2,3))
  for (i in worst.effsize.PI.index){
    coda::traceplot(mcmc(PI[,i]),main=paste('Traceplot of Pi',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # autocorrelation plot for PI
  png(paste0('autocorrelation_PI_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.PI.index){
    acf(mcmc(PI[,i]),main='')
    mtext(sprintf(paste('Autocorrelation plot of PI',i)), side=3,cex=2,line = 1,font=2)
  }
  dev.off()

  # density plot for PI
  png(paste0('density_PI_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.PI.index){
    densplot(mcmc(PI[,i]),main=paste('Density plot of Pi',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # effetive sampling size for S
  eff.size.S<-matrix(NA_real_,nrow=1,ncol=ncol(mc$pars$S))
  for(i in 1:ncol(mc$pars$S)){
    eff.size.S[i]<-effectiveSize(S[,i])
  }

  sort.effsize.S<-sort(eff.size.S,index.return=TRUE)
  if (length(sort.effsize.S$x)<6){
    worst.effsize.S.index<-sort.effsize.S$ix
    worst.effsize.S<-sort.effsize.S$x
    print(paste0("The worst effective sample sizes for S ", worst.effsize.S.index, " are ", worst.effsize.S," out of ", length(S[,1])," samples."))
  }else{
    worst.effsize.S.index<-sort.effsize.S$ix[1:6]
    worst.effsize.S<-sort.effsize.S$x[1:6]
    print(paste0("The worst effective sample sizes for S are ", worst.effsize.S.index, " are ", worst.effsize.S," out of ", length(S[,1])," samples."))
  }

  diagnostic$S$eff.size.index<-worst.effsize.S.index
  diagnostic$S$eff.size<-worst.effsize.S

  # traceplot for S
  png(paste0('traceplot_S_',mc$method,'.png'), width = 1400, height = 800)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(2,3))
  for (i in worst.effsize.S.index){
    coda::traceplot(mcmc(S[,i]),main=paste('Traceplot of S',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # autocorrelation plot for S
  png(paste0('autocorrelation_S_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.S.index){
    acf(mcmc(S[,i]),main='')
    mtext(sprintf(paste('Autocorrelation plot of S',i)), side=3,cex=2,line = 1,font=2)
  }
  dev.off()

  # density plot for S
  png(paste0('density_S_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.S.index){
    densplot(mcmc(S[,i]),main=paste('Density plot of S',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # effetive sampling size for Z
  eff.size.Z<-matrix(NA_real_,nrow=1,ncol=ncol(mc$pars$Z))
  for(i in 1:ncol(mc$pars$Z)){
    eff.size.Z[i]<-effectiveSize(Z[,i])
  }

  sort.effsize.Z<-sort(eff.size.Z,index.return=TRUE)
  if (length(sort.effsize.Z$x)<6){
    worst.effsize.Z.index<-sort.effsize.Z$ix
    worst.effsize.Z<-sort.effsize.Z$x
    print(paste0("The worst effective sample sizes for Z are ", worst.effsize.Z.index, " are ", worst.effsize.Z," out of ", length(Z[,1])," samples."))
  }else{
    worst.effsize.Z.index<-sort.effsize.Z$ix[1:6]
    worst.effsize.Z<-sort.effsize.Z$x[1:6]
    print(paste0("The worst effective sample sizes for Z are ", worst.effsize.Z.index, " are ", worst.effsize.Z," out of ", length(Z[,1])," samples."))
  }

  diagnostic$Z$eff.size.index<-worst.effsize.Z.index
  diagnostic$Z$eff.size<-worst.effsize.Z


  # traceplot for Z
  png(paste0('traceplot_Z_',mc$method,'.png'), width = 1400, height = 800)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(2,3))
  for (i in worst.effsize.Z.index){
    coda::traceplot(mcmc(Z[,i]),main=paste('Traceplot of Z',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  # autocorrelation plot for Z
  png(paste0('autocorrelation_Z_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.Z.index){
    acf(mcmc(Z[,i]),main='')
    mtext(sprintf(paste('Autocorrelation plot of Z',i)), side=3,cex=2,line = 1,font=2)
  }
  dev.off()

  # density plot for Z:
  png(paste0('density_Z_',mc$method,'.png'), width = 1400, height = 800)
  par(mfrow=c(2,3))
  for (i in worst.effsize.Z.index){
    densplot(mcmc(Z[,i]),main=paste('Density plot of Z',i),cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  }
  dev.off()

  #########################################N###################################################
  # effetive sampling size for N
  eff.size.N<-effectiveSize(N)
  diagnostic$N$eff.size<-eff.size.N

  print(paste0("The worst effective sample size is ",eff.size.N))


  # traceplot for N
  png(paste0('traceplot_N_',mc$method,'.png'), width = 600, height = 500)
  par(ps = 12, cex = 1, cex.main = 2.5)
  par(mfrow=c(1,1))
  coda::traceplot(mcmc(N),main='Traceplot of N',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  dev.off()

  # autocorrelation plot for N
  png(paste0('autocorrelation_N_',mc$method,'.png'), width = 600, height = 500)
  par(mfrow=c(1,1))
  acf(mcmc(N),main='')
  mtext(sprintf('Autocorrelation plot of N'), side=3,cex=2,line = 1,font=2)
  dev.off()

  # density plot for N
  png(paste0('density_N_',mc$method,'.png'), width = 600, height = 500)
  par(mfrow=c(1,1))
  densplot(mcmc(N),main='Density plot of N',cex.lab=2.5, cex.axis=1.5, cex.main=2.5, cex.sub=2.5)
  dev.off()

  save(diagnostic,file=paste0('core_inference_diagnostic_',mc$method,'.rda'))
}
