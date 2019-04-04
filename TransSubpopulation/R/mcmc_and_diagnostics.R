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
#' @return NULL. MCMC output is written to file.
source.attribution.mcmc	<- function(dobs, dprior, control=list(seed=42, mcmc.n=1e3, verbose=1, outfile='SAMCMCv190327.rda')){
  #library(data.table); library(gtools)
  
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
  #	determine if the sampling probabilities are <1.  
  # If they are NOT, Z will be the same as TRM_OBS, and the algorithm only updates PI  	
  mc$with.sampling	<- dprior[, list(ALL_ONE=all(P==1)), by='SAMPLING_CATEGORY'][, !all(ALL_ONE)] 
  mc$time			<- NA_real_
  # construct look-up table so we know which transmission pair categories need to be updated 
  # at every MCMC iteration
  tmp				<- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_SAMPLING_CATEGORY, REC_SAMPLING_CATEGORY))
  tmp				<- melt(tmp, id.vars='TRM_CAT_PAIR_ID', value.name='SAMPLING_CATEGORY', variable.name='WHO')
  mc$dlu				<- unique(subset(tmp, select=c(WHO, SAMPLING_CATEGORY)))
  mc$dlu[, UPDATE_ID:= seq_len(nrow(mc$dlu))]
  mc$dl				<- merge(mc$dlu, tmp, by=c('WHO','SAMPLING_CATEGORY'))
  setkey(mc$dl, UPDATE_ID)
  mc$dlt			<- dcast.data.table(mc$dl, TRM_CAT_PAIR_ID~WHO, value.var='UPDATE_ID')
  setnames(mc$dlt, c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY'), c('TR_UPDATE_ID','REC_UPDATE_ID'))
  
  #	every transmission category pair needs to be updated at least one as we sweep through the sampling categories
  # I don t think this is guaranteed, hence the check
  if( !all(dobs$TRM_CAT_PAIR_ID %in% sort(unique(mc$dl$TRM_CAT_PAIR_ID))) )
	  stop('Fatal error. Contact the package maintainer with your input data.')
  
  mc$dlt<-merge(mc$dlt,subset(dobs,select=c('TRM_CAT_PAIR_ID','TRM_OBS')),by='TRM_CAT_PAIR_ID')
  setkey(mc$dlt,TRM_CAT_PAIR_ID)
  
  cat<-mc$dlu$SAMPLING_CATEGORY
  tr.id<-mc$dlt$TR_UPDATE_ID
  rec.id<-mc$dlt$TR_UPDATE_ID
  trm.obs<-mc$dlt$TRM_OBS
  
  update.info<-list()
  for (i in 1:max(mc$dl[,UPDATE_ID])){
    update.info[[i]]<-mc$dl[UPDATE_ID==i,TRM_CAT_PAIR_ID]
  }
  
  
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
#   tmp					<- mc$dl[,  list( 	TR_UPDATE_ID= UPDATE_ID[WHO=='TR_SAMPLING_CATEGORY'][1],
# 											REC_UPDATE_ID= UPDATE_ID[WHO=='REC_SAMPLING_CATEGORY'][1]), 
# 										by='TRM_CAT_PAIR_ID']
  # setkey(tmp,TRM_CAT_PAIR_ID) ## the same as mc$dlt
  mc$pars$S[1,]			<- mc$pars$XI[1, mc$dlt$TR_UPDATE_ID] * mc$pars$XI[1, mc$dlt$REC_UPDATE_ID]
  mc$pars$S_LP[1,]		<- mc$pars$XI_LP[1, mc$dlt$TR_UPDATE_ID] + mc$pars$XI_LP[1, mc$dlt$REC_UPDATE_ID]
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
    if(mc$with.sampling & update.count<mc$sweep)
    {
	  # update.info	<- subset(mc$dl, UPDATE_ID==update.count)	# recompute for each sweep + the category and the pair id
	  
    update.cat<-cat[update.count]  
    update.pair<-update.info[[update.count]]
    
    # propose single XI
	  XI.prop		<- XI.curr
	  XI_LP.prop	<- XI_LP.curr
	  tmp			<- subset(dprior, SAMPLING_CATEGORY==update.cat & SAMPLE==sample(mc$nprior,1))
	  XI.prop[ update.count ]	<- tmp$P	  
	  XI_LP.prop[ update.count ]<- tmp$LP
	  # propose all S that involve the one XI from above 
	  S.prop						<- S.curr
	  S_LP.prop						<- S_LP.curr
	  # tmp							<- subset(mc$dlt, TR_UPDATE_ID==update.count | REC_UPDATE_ID==update.count) # recompute for each sweep + decomposition of the pairs
	  # setkey(tmp,TRM_CAT_PAIR_ID)
	  S.prop[update.pair]	<- XI.prop[tr.id[update.pair]] * XI.prop[rec.id[update.pair]]
	  S_LP.prop[tmp$TRM_CAT_PAIR_ID]<- XI_LP.prop[tr.id[update.pair]] + XI_LP.prop[rec.id[update.pair]]
	  # propose all Z that involve a new S
	  Z.prop						<- Z.curr
#	  trm.obs 						<- merge(tmp, dobs, by='TRM_CAT_PAIR_ID')[, TRM_OBS]	 # do it before the loop
	  Z.prop[update.pair]	<- trm.obs[update.pair] + rnbinom(length(update.pair), trm.obs[update.pair], S.prop[update.pair])
	  # propose total of Z
	  N.prop						<- sum(Z.prop)
	  #	calculate MH ratio
	  log.prop.ratio	<- sum(dnbinom(Z.curr[update.pair]-trm.obs[update.pair], size=trm.obs[update.pair], prob=S.curr[update.pair], log=TRUE)) - 
						   sum(dnbinom(Z.prop[update.pair]-trm.obs[update.pair], size=trm.obs[update.pair], prob=S.prop[update.pair], log=TRUE))		
	  log.fc			<- sum(dbinom(trm.obs[update.pair], size=Z.curr[update.pair], prob=S.curr[update.pair], log=TRUE)) +
						   dmultinom(Z.curr, prob=PI.curr, log=TRUE) +
						   dpois(N.curr, lambda=mc$pars$NU, log=TRUE)
	  log.fc.prop		<- sum(dbinom(trm.obs[update.pair], size=Z.prop[update.pair], prob=S.prop[update.pair], log=TRUE)) +
			  			   dmultinom(Z.prop, prob=PI.curr, log=TRUE) +
			  			   dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
	  log.mh.ratio		<- log.fc.prop - log.fc + log.prop.ratio
	  mh.ratio			<- min(1,exp(log.mh.ratio))	
	  #	update
	  mc$curr.it				<- mc$curr.it+1L
	  set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
	  set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
	  set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
	  accept <-as.integer(runif(1) < mh.ratio)
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
	
	# save
	cat('\nWriting csv file to',control$outfile)
	write.csv(pars, row.names=FALSE, file=control$outfile)
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
#' @return NULL. Diagnostic plots and csv files are written to file.
source.attribution.mcmc.diagnostics	<- function(mcmc.file, control=list(burnin.p=0.2, regex_pars='*', credibility.interval=0.95, pdf.plot.n.worst.case.parameters=10, pdf.plot.all.parameters=FALSE, pdf.height.per.par=1.2, outfile.base=gsub('\\.rda','',mcmc.file))){
  #library(coda); library(data.table); library(bayesplot)
  
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
	  p
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
		  scale_y_continuous(label=scales:::percent) +
		  labs(	x='\nNumber of transmission pair categories updated per sampling category',
				y='Acceptance rate\n')
  ggsave(file=paste0(control$outfile.base,'_acceptance_per_updateID.pdf'), w=6, h=6)
  cat('\nAverage acceptance rate= ',subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, round(mean(ACCEPT), d=3)])
  cat('\nUpdate IDs with lowest acceptance rates')
  print( da[order(ACC_RATE)[1:10],] )
  
  # remove burn-in
  cat('\nRemoving burnin in set to ', 100*control$burnin.p,'% of chain, total iterations=',burnin.n)
  it.rm.burnin	<- seq.int(burnin.n,nrow(mc$pars$S))
  pars	<- pars[it.rm.burnin,,drop=FALSE]
    
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
	  worst.pars	<- pars[, ans[order(NEFF)[1:10], VAR]]
	  #	traces
	  cat('\nPlotting traces for worst parameters...')
	  p				<- mcmc_trace(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol=1))
	  pdf(file=paste0(control$outfile.base,'_worst_traces.pdf'), w=7, h=control$pdf.height.per.par*ncol(worst.pars))
	  p
	  dev.off()
	  
	  #	histograms
	  cat('\nPlotting marginal posterior densities for worst parameters...')
	  p				<- mcmc_hist(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol=4))
	  pdf(file=paste0(control$outfile.base,'_worst_marginalposteriors.pdf'), w=10, h=control$pdf.height.per.par*ncol(worst.pars)/4)
	  p
	  dev.off()
	
	  #	autocorrelations
	  cat('\nPlotting autocorrelations for worst parameters...')
	  p				<- mcmc_acf_bar(worst.pars, pars=colnames(worst.pars), facet_args = list(ncol = 1))
	  pdf(file=paste0(control$outfile.base,'_worst_acf.pdf'), w=7, h=control$pdf.height.per.par*ncol(worst.pars))
	  p
	  dev.off()
  }
}


