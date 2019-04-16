lddirichlet_vector	<- function(x, nu){
	ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
	stopifnot(is.finite(ans))
	ans
}

#' @title Source attribution while adjusting for sampling bias
#' @export
#' @import data.table
#' @importFrom gtools rdirichlet
#' @author Xiaoyue Xi, Oliver Ratmann
#' @param dobs Data.table of observed number of transmission events for each transmission pair category.
#' @param dprior Data.table of draws from the prior distribution of sampling probabilities.
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"seed"}{Random number seed, for reproducibility.}
#'  \item{"mcmc.n"}{Guide on the number of MCMC iterations. The actual number of iterations will be slightly larger, and a multiple of the number of iterations needed to complete one MCMC sweep through all free parameters.}
#'  \item{"verbose"}{Flag to switch on/off verbose output.}
#'  \item{"outfile"}{Full path name to which the MCMC output is saved to in RDA format, in an object called \code{mc}.}
#' }
#' @return NULL. MCMC output is written to an RDA file.
source.attribution.mcmc	<- function(dobs, dprior, control=list(seed=42, mcmc.n=1e3, verbose=1, outfile='SAMCMCv190327.rda')){
  #library(data.table); library(gtools)

  dmode <- function(x) {
    den <- density(x, kernel=c("gaussian"))
    (den$x[den$y==max(den$y)])
  }
  #
  # basic checks
  #
  if(!all(dobs$TR_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
	  stop('Did not find prior samples of a sampling category for TR_SAMPLING_CATEGORY')
  if(!all(dobs$REC_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
	  stop('Did not find prior samples of a sampling category for REC_SAMPLING_CATEGORY')

  ptm	<- Sys.time()
  if('seed'%in%names(control))
  {
	  cat('\nSetting seed to',control$seed)
	  set.seed(control$seed)
  }

  #
  # set up mcmc
  #
  mc				<- list()
  #	determine if the sampling probabilities are <1.
  # If they are NOT, Z will be the same as TRM_OBS, and the algorithm only updates PI
  mc$with.sampling	<- dprior[, list(ALL_ONE=all(P==1)), by='SAMPLING_CATEGORY'][, !all(ALL_ONE)]
  mc$time			<- NA_real_
  # construct look-up table so we know which transmission pair categories need to be updated
  # at every MCMC iteration
  tmp				<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_SAMPLING_CATEGORY, REC_SAMPLING_CATEGORY))
  tmp				<- melt(tmp, id.vars='TRM_CAT_PAIR_ID', value.name='SAMPLING_CATEGORY', variable.name='WHO')
  #	make data.table with unique sampling categories
  mc$dlu			<- unique(subset(tmp, select=c(WHO, SAMPLING_CATEGORY)))
  mc$dlu[, UPDATE_ID:= seq_len(nrow(mc$dlu))]
  setkey(mc$dlu, UPDATE_ID)
  # make data.table that maps UPDATE_IDs to TRM_CAT_PAIR_IDs
  mc$dl				<- merge(mc$dlu, tmp, by=c('WHO','SAMPLING_CATEGORY'))
  setkey(mc$dl, UPDATE_ID)
  #	every transmission category pair needs to be updated at least one as we sweep through the sampling categories
  # I don t think this is guaranteed, hence the check
  if( !all(dobs$TRM_CAT_PAIR_ID %in% sort(unique(mc$dl$TRM_CAT_PAIR_ID))) )
	  stop('Fatal error. Contact the package maintainer with your input data.')
  # make data.table that maps TRM_CAT_PAIR_IDs to transmitter UPDATE_IDs and recipient UPDATE_IDs
  mc$dlt			<- dcast.data.table(mc$dl, TRM_CAT_PAIR_ID~WHO, value.var='UPDATE_ID')
  setnames(mc$dlt, c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY'), c('TR_UPDATE_ID','REC_UPDATE_ID'))
  mc$dlt			<- merge(mc$dlt,subset(dobs,select=c('TRM_CAT_PAIR_ID','TRM_OBS')),by='TRM_CAT_PAIR_ID')
  setkey(mc$dlt,TRM_CAT_PAIR_ID)

  # make indexed lookup table for speed
  update.info	<- vector('list', nrow(mc$dlu))
  for (i in 1:nrow(mc$dlu))
  {
	  update.info[[i]]	<- mc$dl[UPDATE_ID==i,TRM_CAT_PAIR_ID]
  }
  # make indexed prior samples for speed
  setkey(dprior,SAMPLING_CATEGORY)
  dprior2		<- vector('list', nrow(mc$dlu))
  for (i in 1:nrow(mc$dlu))
  {
    tmp				<- mc$dlu[UPDATE_ID==i,SAMPLING_CATEGORY]
    dprior2[[i]]<-dprior[J(tmp),nomatch=0L]
  }

  dprior3   <- dprior[,list(EST_SAMPLING_RATE=dmode(P)),by=SAMPLING_CATEGORY]
  dobs2   <- subset(dobs,select = c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY','TRM_CAT_PAIR_ID'))
  setnames(dprior3,colnames(dprior3),paste0('TR_',colnames(dprior3)))
  dobs2   <- merge(dobs2,dprior3,by='TR_SAMPLING_CATEGORY')
  setnames(dprior3,colnames(dprior3),gsub('TR_','REC_',colnames(dprior3)))
  dobs2   <- merge(dobs2,dprior3,by='REC_SAMPLING_CATEGORY')
  dobs2[,EST_SAMPLING_RATE:=TR_EST_SAMPLING_RATE * REC_EST_SAMPLING_RATE]

  mc$nprior			<- max(dprior$SAMPLE)
  mc$sweep			<- nrow(mc$dlu)+1L
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
  mc$pars$S[1,]			<- mc$pars$XI[1, mc$dlt$TR_UPDATE_ID] * mc$pars$XI[1, mc$dlt$REC_UPDATE_ID]
  mc$pars$S_LP[1,]		<- mc$pars$XI_LP[1, mc$dlt$TR_UPDATE_ID] + mc$pars$XI_LP[1, mc$dlt$REC_UPDATE_ID]
  #	augmented data: proposal draw under sampling probability
  mc$pars$Z[1,]			<- dobs$TRM_OBS + rnbinom(nrow(dobs),dobs$TRM_OBS,mc$pars$S[1,])
  #	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
  mc$pars$NU			<- sum(dobs$TRM_OBS) / mean(dobs2$EST_SAMPLING_RATE)
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
  for(i in 1L:mc$n)
  {
    mc$curr.it		<- i
    # determine source-recipient combination that will be updated in this iteration
    update.count	<- (i-1L) %% mc$sweep + 1L
    update.round 	<- (i-1L) %/% mc$sweep + 1L
	# update in one go S, Z, N for the ith XI
    if(mc$with.sampling & update.count<mc$sweep)
    {
	  # update.info	<- subset(mc$dl, UPDATE_ID==update.count)	# recompute for each sweep + the category and the pair id
      update.cat	<- mc$dlu$SAMPLING_CATEGORY[update.count]
      update.pairs	<- update.info[[update.count]]

      # propose single XI
	  XI.prop		<- XI.curr
	  XI_LP.prop	<- XI_LP.curr
	  tmp			<- dprior2[[update.count]][sample(mc$nprior,1),]
	  if(tmp$SAMPLING_CATEGORY[1]!=update.cat)
		  stop('\nFatal error in dprior2.')
	  XI.prop[ update.count ]	<- tmp$P
	  XI_LP.prop[ update.count ]<- tmp$LP
	  # propose all S that involve the one XI from above
	  S.prop						<- S.curr
	  S_LP.prop						<- S_LP.curr
	  S.prop[update.pairs]			<- XI.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] * XI.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
	  S_LP.prop[update.pairs]		<- XI_LP.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] + XI_LP.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
	  # propose all Z that involve a new S
	  Z.prop						<- Z.curr
	  Z.prop[update.pairs]			<- mc$dlt$TRM_OBS[update.pairs] + rnbinom(length(update.pairs), mc$dlt$TRM_OBS[update.pairs], S.prop[update.pairs])
	  # propose total of Z
	  N.prop						<- sum(Z.prop)
	  #	calculate MH ratio
	  log.prop.ratio	<- sum(dnbinom(Z.curr[update.pairs]-mc$dlt$TRM_OBS[update.pairs], size=mc$dlt$TRM_OBS[update.pairs], prob=S.curr[update.pairs], log=TRUE)) -
						   sum(dnbinom(Z.prop[update.pairs]-mc$dlt$TRM_OBS[update.pairs], size=mc$dlt$TRM_OBS[update.pairs], prob=S.prop[update.pairs], log=TRUE))
	  log.fc			<- sum(dbinom(mc$dlt$TRM_OBS[update.pairs], size=Z.curr[update.pairs], prob=S.curr[update.pairs], log=TRUE)) +
						   dmultinom(Z.curr, prob=PI.curr, log=TRUE) +
						   dpois(N.curr, lambda=mc$pars$NU, log=TRUE)
	  log.fc.prop		<- sum(dbinom(mc$dlt$TRM_OBS[update.pairs], size=Z.prop[update.pairs], prob=S.prop[update.pairs], log=TRUE)) +
			  			   dmultinom(Z.prop, prob=PI.curr, log=TRUE) +
			  			   dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
	  log.mh.ratio		<- log.fc.prop - log.fc + log.prop.ratio
	  mh.ratio			<- min(1,exp(log.mh.ratio))
	  #	update
	  mc$curr.it		<- mc$curr.it+1L
	  set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
	  set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
	  set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
	  accept 			<- as.integer(runif(1) < mh.ratio)
	  set(mc$it.info, mc$curr.it, 'ACCEPT', accept)

	  if(control$verbose & accept)
	  {
		  cat('\nit ',mc$curr.it,' ACCEPT S-Z-N block ',update.count)
	  }
	  if(accept)
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
    if(!mc$with.sampling | update.count==mc$sweep)
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
	if(update.count==mc$sweep & update.round %% 100 == 0){
		cat('\nSweeps done:\t',update.round)
	}
  }

  mc$time	<- Sys.time()-ptm
  if(!'outfile'%in%names(control))
	  return(mc)
  save(mc,	file=control$outfile)
  NULL
}

#' @title Aggregate MCMC output to target parameters
#' @export
#' @import data.table
#' @param mcmc.file Full file name to MCMC output from function \code{source.attribution.mcmc}
#' @param daggregateTo Data.table that maps the categories of transmission pairs used in the MCMC to lower-dimensional categories that are of primary interest.
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"burnin.p"}{Proportion of MCMC iterations that are removed as burn-in period.}
#'  \item{"thin"}{Thin MCMC output to every nth MCMC iteration.}
#'  \item{"regex_pars"}{Regular expression to select the parameters that are to be aggregated.}
#'  \item{"outfile"}{Full file name of the output csv file.}
#' }
#' @return NULL. Monte Carlo samples are written to a csv file.
source.attribution.mcmc.aggregateToTarget	<- function(mcmc.file, daggregateTo, control=list(burnin.p=NA_real_, thin=NA_integer_, regex_pars='*', outfile=gsub('\\.rda','_aggregated.csv',mcmc.file))){
	#	basic checks
	if(!'data.table'%in%class(daggregateTo))
		stop('daggregateTo is not a data.table')
	if(!all(c('TRM_CAT_PAIR_ID','TR_TARGETCAT','REC_TARGETCAT')%in%colnames(daggregateTo)))
		stop('daggregateTo does not contain one of the required columns TRM_CAT_PAIR_ID, TR_TARGETCAT, REC_TARGETCAT')

	#	load MCMC output
	cat('\nLoading MCMC output...')
	load(mcmc.file)

	#	define internal control variables
	burnin.p	<- control$burnin.p
	if(is.na(burnin.p))
		burnin.p<- 0
	burnin.n	<- floor(burnin.p*nrow(mc$pars$S))
	thin		<- control$thin
	if(is.na(thin))
		thin	<- 1

	#	collect parameters
	cat('\nCollecting parameters...')
	pars		<- matrix(NA,nrow=nrow(mc$pars$S),ncol=0)
	if(grepl(control$regex_pars,'Z'))
	{
		tmp	<- mc$pars$Z
		colnames(tmp)	<- paste0('Z-',1:ncol(tmp))
		pars	<- cbind(pars, tmp)
	}
	if(grepl(control$regex_pars,'PI'))
	{
		tmp	<- mc$pars$PI
		colnames(tmp)	<- paste0('PI-',1:ncol(tmp))
		pars	<- cbind(pars, tmp)
	}

	# remove burn-in
	if(burnin.n>0)
	{
		cat('\nRemoving burnin in set to ', 100*burnin.p,'% of chain, total iterations=',burnin.n)
		tmp		<- seq.int(burnin.n,nrow(mc$pars$S))
		pars	<- pars[tmp,,drop=FALSE]
	}

	# thin
	if(thin>1)
	{
		cat('\nThinning to every', thin,'th iteration')
		tmp		<- seq.int(1,nrow(pars),thin)
		pars	<- pars[tmp,,drop=FALSE]
	}

	cat('\nMaking aggregated MCMC output...')
	# make data.table in long format
	pars	<- as.data.table(pars)
	pars[, SAMPLE:= seq_len(nrow(pars))]
	pars	<- melt(pars, id.vars='SAMPLE')
	pars[, VARIABLE:= pars[, gsub('([A-Z]+)-([0-9]+)','\\1',variable)]]
	pars[, TRM_CAT_PAIR_ID:= pars[, as.integer(gsub('([A-Z]+)-([0-9]+)','\\2',variable))]]

	# aggregate MCMC samples
	if(!all(sort(unique(pars$TRM_CAT_PAIR_ID))==sort(unique(daggregateTo$TRM_CAT_PAIR_ID))))
		stop('The transmission count categories in the MCMC output do not match the transmission count categories in the aggregateTo data table.')
	pars	<- merge(pars, daggregateTo, by='TRM_CAT_PAIR_ID')
	pars	<- pars[, list(VALUE=sum(value)), by=c('VARIABLE','TR_TARGETCAT','REC_TARGETCAT','SAMPLE')]

	# save or return
	if(!'outfile'%in%names(control))
		return(pars)
	if(grepl('csv$',control$outfile))
	{
		cat('\nWriting csv file to',control$outfile)
		write.csv(pars, row.names=FALSE, file=control$outfile)
	}
	if(grepl('rda$',control$outfile))
	{
		cat('\nSaving rda file to',control$outfile)
		save(pars, file=control$outfile)
	}
}

#' @title Estimate Flows, Sources, WAIFM, Flow ratios
#' @export
#' @import data.table
#' @author Xiaoyue Xi, Oliver Ratmann
#' @param infile Full file name to aggregated MCMC output from function \code{source.attribution.mcmc}
#' @param control List of input arguments that control the behaviour of the derivation of key quantities:
#' \itemize{
#'  \item{"quantiles"}{Named list of quantiles. Default: c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975)}
#'  \item{"flowratios"}{Cector of length 3. First element: name of flow ratio. Second element: name of transmission pair category for enumerator of flow ratio. Third element: name of transmission pair category for denominator of flow ratio.}
#'  \item{"outfile"}{Full file name for output csv file.}
#' }
#' @return NULL. A csv file is written to disk.
source.attribution.aggmcmc.getKeyQuantities<- function(infile, control)
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


#' @title MCMC diagnostics for the source attribution algorithm
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom bayesplot mcmc_trace mcmc_acf_bar mcmc_hist
#' @importFrom coda mcmc effectiveSize
#' @author Xiaoyue Xi, Oliver Ratmann
#' @param mcmc.file Full file name to MCMC output from function \code{source.attribution.mcmc}
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"burnin.p"}{Proportion of MCMC iterations that are removed as burn-in period.}
#'  \item{"regex_pars"}{Regular expression to select the parameter names for which diagnostics are computed. The default is all parameters, which is '*'.}
#'  \item{"credibility.interval"}{Width of the marginal posterior credibility intervals.}
#'  \item{"pdf.plot.n.worst.case.parameters"}{Integer which specifies the number of parameters with smallest effective sample size that are inspected in detail. If set to 0, worst case analyses are not performed.}
#'  \item{"pdf.plot.all.parameters"}{Flag which specifies if traces shall be plotted for all parameters. If set to TRUE, very large pdfs may be created.}
#'  \item{"pdf.height.per.par"}{Most plots show diagnostics with parameters listed on the y-axis. This value controls the plot height in inches for each free parameter.}
#'  \item{"outfile.base"}{Start of the full file name for all output files.}
#' }
#' @return NULL. Diagnostic plots and csv files are written to disk.
source.attribution.mcmc.diagnostics	<- function(mcmc.file, control=list(burnin.p=0.2, regex_pars='*', credibility.interval=0.95, pdf.plot.n.worst.case.parameters=10, pdf.plot.all.parameters=FALSE, pdf.height.per.par=1.2, outfile.base=gsub('\\.rda','',mcmc.file))){
  #library(coda); library(data.table); library(bayesplot); library(ggplot2)

  cat('\nLoading MCMC output...')
  load(mcmc.file)

  burnin.n	<- floor(control$burnin.p*nrow(mc$pars$S))

  cat('\nCollecting parameters...')
  pars		<- matrix(NA,nrow=nrow(mc$pars$S),ncol=0)
  if(grepl(control$regex_pars,'S'))
  {
	  tmp	<- mc$pars$S
	  colnames(tmp)	<- paste0('S-',1:ncol(tmp))
	  pars	<- cbind(pars, tmp)
  }
  if(grepl(control$regex_pars,'Z'))
  {
	  tmp	<- mc$pars$Z
	  colnames(tmp)	<- paste0('Z-',1:ncol(tmp))
	  pars	<- cbind(pars, tmp)
  }
  if(grepl(control$regex_pars,'N'))
  {
	  tmp	<- mc$pars$N
	  colnames(tmp)	<- paste0('N-',1:ncol(tmp))
	  pars	<- cbind(pars, tmp)
  }
 if(grepl(control$regex_pars,'PI'))
 {
	 tmp	<- mc$pars$PI
	 colnames(tmp)	<- paste0('PI-',1:ncol(tmp))
	 pars	<- cbind(pars, tmp)
 }

  #	traces for parameters
  if(control$pdf.plot.all.parameters)
  {
	  cat('\nPlotting traces for all parameters...')
	  p		<- mcmc_trace(pars, pars=colnames(pars), facet_args = list(ncol = 1), n_warmup=burnin.n)
	  pdf(file=paste0(control$outfile.base,'_marginaltraces.pdf'), w=7, h=control$pdf.height.per.par*ncol(pars))
	  print(p)
	  dev.off()
  }

  #	acceptance rate per MCMC update ID
	cat('\nPlotting acceptance rates...')
  da	<- subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, list(ACC_RATE=mean(ACCEPT)), by='PAR_ID']
  setnames(da, 'PAR_ID', 'UPDATE_ID')
  tmp	<- mc$dl[, list(N_TRM_CAT_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
  da	<- merge(da, tmp, by='UPDATE_ID')
  ggplot(da, aes(x=N_TRM_CAT_PAIRS, y=ACC_RATE)) +
		  geom_point() +
		  theme_bw() +
		  scale_y_continuous(label=scales::percent) +
		  labs(	x='\nNumber of transmission pair categories updated per sampling category',
				y='Acceptance rate\n')
  ggsave(file=paste0(control$outfile.base,'_acceptance_per_updateID.pdf'), w=6, h=6)
  cat('\nAverage acceptance rate= ',subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, round(mean(ACCEPT), d=3)])
  cat('\nUpdate IDs with lowest acceptance rates')
  print( da[order(ACC_RATE)[1:10],] )

  # remove burn-in
  cat('\nRemoving burnin in set to ', 100*control$burnin.p,'% of chain, total iterations=',burnin.n)
  tmp	<- seq.int(burnin.n+1L,nrow(mc$pars$S))
  pars	<- pars[tmp,,drop=FALSE]

  # effective sampling sizes
  cat('\nCalculating effective sample size for all parameters...')
  tmp	<- mcmc(pars)
  ans	<- data.table(ID= seq_len(ncol(pars)), VAR= colnames(pars), NEFF=as.numeric(effectiveSize(tmp)))
  set(ans, NULL, 'ID', ans[, factor(ID, labels=VAR)])
  if(control$pdf.plot.all.parameters)
  {
	  ggplot(ans, aes(x=NEFF, y=ID)) +
			  geom_point() +
			  theme_bw() +
			  labs(x='\neffective sample size', y='')
	  ggsave(file=paste0(control$outfile.base,'_neff.pdf'), w=6, h=control$pdf.height.per.par*ncol(pars)*0.15, limitsize=FALSE)
  }

  # summarise mean, sd, quantiles
  cat('\nCalculating posterior summaries for all parameters...')
  tmp	<- apply(pars, 2, function(x) quantile(x, p=c((1-control$credibility.interval)/2, 0.5, control$credibility.interval+(1-control$credibility.interval)/2)))
  tmp	<- data.table(	VAR= colnames(pars),
		  				MEAN= apply(pars, 2, mean),
						SD= apply(pars, 2, sd),
		  				MEDIAN=tmp[2,],
						CI_L=tmp[1,],
						CI_U=tmp[3,])
  ans	<- merge(tmp, ans, by='VAR')
  cat('\nParameters with lowest effective samples')
  print( ans[order(NEFF)[1:10],] )

  # write to file
  cat('\nWriting summary file to',paste0(control$outfile.base,'_summary.csv'))
  setkey(ans, ID)
  write.csv(ans, file=paste0(control$outfile.base,'_summary.csv'))

  # plots for worst case parameters
  if(control$pdf.plot.n.worst.case.parameters>0)
  {
	  worst.pars	<- pars[, ans[order(NEFF)[1:control$pdf.plot.n.worst.case.parameters], VAR]]
	  #	traces
	  cat('\nPlotting traces for worst parameters...')
	  p				<- mcmc_trace(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol=1))
	  pdf(file=paste0(control$outfile.base,'_worst_traces.pdf'), w=7, h=control$pdf.height.per.par*ncol(worst.pars))
	  print(p)
	  dev.off()

	  #	histograms
	  cat('\nPlotting marginal posterior densities for worst parameters...')
	  p				<- mcmc_hist(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol=4))
	  pdf(file=paste0(control$outfile.base,'_worst_marginalposteriors.pdf'), w=10, h=control$pdf.height.per.par*ncol(worst.pars)/4)
	  print(p)
	  dev.off()

	  #	autocorrelations
	  cat('\nPlotting autocorrelations for worst parameters...')
	  p				<- mcmc_acf_bar(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol = 1))
	  pdf(file=paste0(control$outfile.base,'_worst_acf.pdf'), w=7, h=control$pdf.height.per.par*ncol(worst.pars))
	  print(p)
	  dev.off()
  }
}


