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
#' @return If `outfile` is not specified, phyloflows MCMC output is returned, which is a list object. If `outfile` is specified, phyloflows MCMC output is written to an RDA file.
#' @seealso \link{\code{source.attribution.mcmc.diagnostics}}, \link{\code{source.attribution.mcmc.aggregateToTarget}}, \link{\code{source.attribution.mcmc.getKeyQuantities}}
#' @examples
#' require(data.table)
#' require(phyloflows)
#' # load input data set
#' data(twoGroupFlows1, package="phyloflows")
#' dobs <- twoGroupFlows1$dobs
#' dprior <- twoGroupFlows1$dprior
#' 
#' #
#' # run MCMC and return output to R
#' # this is usually fine for simple source attribution models
#' # that don t take long to run
#' #
#' \dontrun{
#' control <- list(seed=42,    # seed for reproducibility
#'                 mcmc.n=5e4, # guide on the number of MCMC iterations, the actual number of MCMC iterations may be a little bit larger
#'                 verbose=0   # switch on for output to each MCMC iteration
#'                 )
#' mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)
#' }
#' 
#' #
#' # run MCMC and save output to file
#' # this is recommended when the source attribution is more complex
#' # so you can deploy the MCMC on a high performance cluster
#' #
#' \dontrun{
#' mcmc.file <- 'mcmcm_output.rda'
#' control <- list(seed=42,    # seed for reproducibility
#'                 mcmc.n=5e4, # guide on the number of MCMC iterations, the actual number of MCMC iterations may be a little bit larger
#'                 verbose=0   # switch on for output to each MCMC iteration
#'                 outfile= mcmc.file # store MCMC output as rda file here
#'                 )
#' mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)
#' }
#' 
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
    dprior2[[i]]	<- dprior[J(tmp),nomatch=0L]
  }
  if(dprior[, max(SAMPLE)]>10)
  {
	  dprior3   <- dprior[,list(EST_SAMPLING_RATE=dmode(P)),by=SAMPLING_CATEGORY]  
  }
  if(dprior[, max(SAMPLE)]<=10)
  {
	  dprior3   <- dprior[,list(EST_SAMPLING_RATE=median(P)),by=SAMPLING_CATEGORY]
  }
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
