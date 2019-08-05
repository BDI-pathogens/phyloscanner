#' @title Source attribution while adjusting for sampling bias
#' @export
#' @import data.table Matrix
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
  library(data.table); library(gtools);library(Matrix)
  
  dmultinom_log  <-  function (x, size = NULL, log_prob, log = FALSE) 
  {
    K <- length(log_prob)
    N <- sum(x)
    if (is.null(size)) {
      size <- N
    } else if(size != N){
      stop("size != sum(x), i.e. one is wrong")}
    r <- lgamma(size + 1) + sum(x * log_prob - lgamma(x + 1))
    if (log) 
      r
    else exp(r)
  }
  
  
  rdirichlet_log<-function(alpha) 
  {
    source('rgamss.R')
    l <- length(alpha)
    log_x <- unlist(lapply(alpha, function(x){rgamss(1,x,do.log=TRUE)}))
    if (any(log_x<-500)){
      z <- 500
      sm <- sum(exp(z+log_x))/exp(z)}
    else{
      sm <- sum(exp(log_x))
    }
    log_x-log(sm)
  }
  
  lddirichlet_log_vector	<- function(log_x, nu){
    ans	<- sum((nu - 1) * log_x) + lgamma(sum(nu)) - sum(lgamma(nu)) 
    stopifnot(is.finite(ans))
    ans
  }
  
  lddirichlet_vector	<- function(x, nu){
    ans	<- sum((nu - 1) * log(x)) + lgamma(sum(nu)) - sum(lgamma(nu)) 
    stopifnot(is.finite(ans))
    ans
  }
  
  zproposal_generator<- function(n,p){
    x <- ceiling( log(1-runif(n))/log(1-p) -1)
    x
  }
  
  zproposal_log_density <- function(x,p){
    d <- log(p) + x*log(1-p)
  }
  
  zproposal_density <- function(x,p){
    d <- p * (1-p)^x
  }
  
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
  mc$dlu			<- unique(subset(tmp, select=c(SAMPLING_CATEGORY)))
  mc$dlu[, UPDATE_ID:= seq_len(nrow(mc$dlu))]
  setkey(mc$dlu, UPDATE_ID)
  # make data.table that maps UPDATE_IDs to TRM_CAT_PAIR_IDs
  mc$dl				<- merge(mc$dlu, tmp, by=c('SAMPLING_CATEGORY'))
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
  mc$sweep			<- nrow(mc$dlu)+sum(lengths(update.info))+1L
  mc$nsweep			<- ceiling( control$mcmc.n/mc$sweep )
  mc$n				<- mc$nsweep*mc$sweep
  # every mc$sweep_group sweep, we reset mc$pars and mc$it.info
  mc$sweep_group <- 1e4L
  mc$pars			<- list()
  mc$pars$LAMBDA	<- matrix(NA_real_, ncol=nrow(dobs), nrow=1)		#prior for proportions
  mc$pars$XI		<- matrix(NA_real_, ncol=length(update.info), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # prior for sampling in transmitter categories and sampling in recipient categories, concatenated
  mc$pars$XI_LP		<- matrix(NA_real_, ncol=length(update.info), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated
  mc$pars$S			<- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # prior probability of sampling transmission pair categories
  mc$pars$S_LP		<- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # log prior density of sampling transmission pair categories
  mc$pars$Z			<- Matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$nsweep+1L,sparse = TRUE) #augmented data
  mc$pars$NU		<- NA_real_															#prior for N
  mc$pars$N			<- matrix(NA_integer_, ncol=1, nrow=min(mc$nsweep+1L,mc$sweep_group+1L))							#total number of counts on augmented data
  mc$pars$LOG_PI		<- matrix(NA_real_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L))	#proportions
  
  mc$it.info		<- data.table(  IT= seq.int(0,min(mc$nsweep*mc$sweep,mc$sweep_group*mc$sweep)),
                             PAR_ID= rep(NA_integer_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                             BLOCK= rep(NA_character_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                             MHRATIO= rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                             ACCEPT=rep(NA_integer_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                             LOG_LKL=rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                             LOG_PRIOR=rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)))
  
  if(1)
  {
    cat('\nNumber of parameters:\t', ncol(mc$pars$LOG_PI)+ncol(mc$pars$N)+ncol(mc$pars$Z)+ncol(mc$pars$XI) )
    cat('\nDimension of PI:\t', ncol(mc$pars$LOG_PI))
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
  update.sweep.group <- 0L
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
  # prior for sampling in recipient categories
  tmp					<- subset(dprior, SAMPLE==sample(mc$nprior,1))
  tmp					<- merge(unique(subset(mc$dl, WHO=='REC_SAMPLING_CATEGORY', c(SAMPLING_CATEGORY, UPDATE_ID))), tmp, by='SAMPLING_CATEGORY')
  setkey(tmp, UPDATE_ID)
  mc$pars$XI[1, tmp$UPDATE_ID]		<- tmp$P
  mc$pars$XI_LP[1, tmp$UPDATE_ID]	<- tmp$LP
  # prior for sampling in transmission pair categories
  mc$pars$S[1,]			<- mc$pars$XI[1, mc$dlt$TR_UPDATE_ID] * mc$pars$XI[1, mc$dlt$REC_UPDATE_ID]
  mc$pars$S_LP[1,]		<- mc$pars$XI_LP[1, mc$dlt$TR_UPDATE_ID] + mc$pars$XI_LP[1, mc$dlt$REC_UPDATE_ID]
  #	augmented data: proposal draw under sampling probability
  dobs.zero  <-  which(dobs$TRM_OBS==0)
  dobs.nonzero  <-  which(dobs$TRM_OBS!=0)
  mc$pars$Z[1,dobs.nonzero]			<-   dobs$TRM_OBS[dobs.nonzero] + 
    rnbinom(length(dobs.nonzero),dobs$TRM_OBS[dobs.nonzero],mc$pars$S[1,dobs.nonzero])
  mc$pars$Z[1,dobs.zero]	 <-   unlist( lapply(mc$pars$S[1,dobs.zero],
                                              function(p){zproposal_generator(1,p)}))
  
  
  #	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
  mc$pars$NU			<- sum(dobs$TRM_OBS) / mean(dobs2$EST_SAMPLING_RATE) #  1276.456
  
  # dobs2<-merge(dobs2,subset(dobs,select = c('TRM_OBS', 'TRM_CAT_PAIR_ID')),by=c('TRM_CAT_PAIR_ID'))
  # dobs2[TRM_OBS!=0,EST_Z:=TRM_OBS/EST_SAMPLING_RATE]
  # dobs2[TRM_OBS==0,EST_Z:=(1-EST_SAMPLING_RATE)/EST_SAMPLING_RATE]
  # sum(dobs2[,TRM_OBS/EST_SAMPLING_RATE]) #1194.338
  # mc$pars$NU			<- sum(dobs2$EST_Z) # 1850.857
  
  #	total count: that s just the sum of Z
  mc$pars$N[1,]			<- sum(mc$pars$Z[1,])
  #	proportions: draw from full conditional
  mc$pars$LOG_PI[1,]		<- rdirichlet_log(mc$pars$Z[1,] + mc$pars$LAMBDA[1,])
  #	store log likelihood
  tmp	<- sum( dbinom(dobs$TRM_OBS, size=mc$pars$Z[1,], prob=mc$pars$S[1,], log=TRUE) ) +
    dmultinom_log(mc$pars$Z[1,], size=mc$pars$N[1,], log_prob=mc$pars$LOG_PI[1,], log=TRUE)
  set(mc$it.info, 1L, 'LOG_LKL', tmp)
  # 	store log prior
  tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
    lddirichlet_log_vector(mc$pars$LOG_PI[1,], nu=mc$pars$LAMBDA[1,]) +
    sum(mc$pars$S_LP[1,])
  set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
  
  
  # parameter value at the current step
  XI.curr	<- mc$pars$XI[1,]
  XI_LP.curr<- mc$pars$XI_LP[1,]
  S.curr	<- mc$pars$S[1,]
  S_LP.curr	<- mc$pars$S_LP[1,]
  LOG_PI.curr	<- mc$pars$LOG_PI[1,]
  N.curr	<- mc$pars$N[1,]
  Z.curr	<- mc$pars$Z[1,]
  
  # determine the iteration in a sweep which updates XI
  update.xi  <- c(1,head(cumsum(1+lengths(update.info))+1,-1))
  
  # run mcmc
  options(warn=0)
  for(i in 1L:mc$n)
  {
    mc$curr.it		<- i
    # determine the current iteration in a sweep
    update.count	<- (i-1L) %% mc$sweep + 1L
    update.round 	<- (i-1L) %/% mc$sweep + 1L
    
    
    # update XI first
    if(mc$with.sampling & update.count %in% update.xi)
    {
      # determine the category to update
      update.group  <-  which(update.xi==update.count)
      update.cat	<- mc$dlu$SAMPLING_CATEGORY[update.group]
      
      # choose the pairs to update
      update.pairs	<- update.info[[update.group]]
      
      # record the current iteration in a sweep which updates XI
      # this help to identiy which pair to update in the next if loop
      update.xi.count <-  update.count
      
      # propose single XI
      XI.prop		<- XI.curr
      XI_LP.prop	<- XI_LP.curr
      tmp			<- dprior2[[update.group]][sample(mc$nprior,1),]
      if(tmp$SAMPLING_CATEGORY[1]!=update.cat)
        stop('\nFatal error in dprior2.')
      XI.prop[ update.group ]	<- tmp$P
      XI_LP.prop[ update.group ]<- tmp$LP
      
      # propose all S that involve the one XI from above
      S.prop						<- S.curr
      S_LP.prop						<- S_LP.curr
      S.prop[update.pairs]			<- XI.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] * XI.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
      S_LP.prop[update.pairs]		<- XI_LP.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] + XI_LP.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
      # calculate MH ratio
      log.fc			<- sum(dbinom(mc$dlt$TRM_OBS[update.pairs], size=Z.curr[update.pairs], prob=S.curr[update.pairs], log=TRUE))
      log.fc.prop		<- sum(dbinom(mc$dlt$TRM_OBS[update.pairs], size=Z.curr[update.pairs], prob=S.prop[update.pairs], log=TRUE)) 
      log.mh.ratio		<- log.fc.prop - log.fc 
      mh.ratio			<- min(1,exp(log.mh.ratio))
      #	update
      mc$curr.it		<- mc$curr.it+1L
      
      # update index for mc$it.info
      if(update.sweep.group==0L){
        mc$curr.it.adj <- mc$curr.it
      }else{
        mc$curr.it.adj <- mc$curr.it-mc$sweep_group*mc$sweep*update.sweep.group-1L
      }
      
      set(mc$it.info, mc$curr.it.adj, 'BLOCK', 'XI')
      set(mc$it.info, mc$curr.it.adj, 'PAR_ID', update.group)
      set(mc$it.info, mc$curr.it.adj, 'MHRATIO', mh.ratio)
      accept 			<- as.integer(runif(1) < mh.ratio)
      set(mc$it.info, mc$curr.it.adj, 'ACCEPT', accept)
      
      if(control$verbose & accept)
      {
        cat('\nit ',mc$curr.it,' ACCEPT XI ',update.group)
      }
      if(accept)
      {
        XI.curr	<- XI.prop
        XI_LP.curr<- XI_LP.prop
        S.curr	<- S.prop
        S_LP.curr	<- S_LP.prop
      }
      
    }
    
    if(mc$with.sampling & !update.count %in% update.xi & update.count<mc$sweep)
    {
      # given the update.group, determine the pair to update in the current iteration
      tmp   <- update.count-update.xi.count 
      update.pairs.count  <-  update.info[[update.group]][tmp]
      
      # propose Z for the pair
      Z.prop						<- Z.curr
      if (mc$dlt$TRM_OBS[update.pairs.count]!=0){
        Z.prop[update.pairs.count]			<- mc$dlt$TRM_OBS[update.pairs.count] + 
          rnbinom(1, mc$dlt$TRM_OBS[update.pairs.count], S.curr[update.pairs.count])
      }else{
        Z.prop[update.pairs.count]	<- zproposal_generator(1,S.curr[update.pairs.count])
      }
      
      # propose total of Z
      N.prop						<- sum(Z.prop)
      #	calculate MH ratio
      if (mc$dlt$TRM_OBS[update.pairs.count]!=0){
        log.prop.ratio  <- sum(dnbinom(Z.curr[update.pairs.count]-mc$dlt$TRM_OBS[update.pairs.count], size=mc$dlt$TRM_OBS[update.pairs.count], prob=S.curr[update.pairs.count], log=TRUE))-
          sum(dnbinom(Z.prop[update.pairs.count]-mc$dlt$TRM_OBS[update.pairs.count], size=mc$dlt$TRM_OBS[update.pairs.count], prob=S.curr[update.pairs.count], log=TRUE))
      }else{
        log.prop.ratio	<- zproposal_log_density(Z.curr[update.pairs.count],S.curr[update.pairs.count]) -
          zproposal_log_density(Z.prop[update.pairs.count],S.curr[update.pairs.count]) 
      }
      log.fc			<- sum(dbinom(mc$dlt$TRM_OBS[update.pairs.count], size=Z.curr[update.pairs.count], prob=S.curr[update.pairs.count], log=TRUE)) +
        dmultinom_log(Z.curr, log_prob=LOG_PI.curr, log=TRUE) +
        dpois(N.curr, lambda=mc$pars$NU, log=TRUE)
      
      log.fc.prop		<- sum(dbinom(mc$dlt$TRM_OBS[update.pairs.count], size=Z.prop[update.pairs.count], prob=S.curr[update.pairs.count], log=TRUE)) +
        dmultinom_log(Z.prop, log_prob=LOG_PI.curr, log=TRUE) +
        dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
      
      log.mh.ratio		<- log.fc.prop - log.fc + log.prop.ratio
      mh.ratio			<- min(1,exp(log.mh.ratio))
      
      #	update
      mc$curr.it		<- mc$curr.it+1L
      
      # set up index for mc$it.info
      if(update.sweep.group==0L){
        mc$curr.it.adj <- mc$curr.it
      }else{
        mc$curr.it.adj <- mc$curr.it-mc$sweep_group*mc$sweep*update.sweep.group-1L
      }
      
      set(mc$it.info, mc$curr.it.adj, 'BLOCK', 'Z-N')
      set(mc$it.info, mc$curr.it.adj, 'PAR_ID', update.count)
      set(mc$it.info, mc$curr.it.adj, 'MHRATIO', mh.ratio)
      accept 			<- as.integer(runif(1) < mh.ratio)
      set(mc$it.info, mc$curr.it.adj, 'ACCEPT', accept)
      
      
      if(control$verbose & accept)
      {
        cat('\nit ',mc$curr.it,' ACCEPT Z-N block ',update.count)
      }
      if(accept)
      {
        Z.curr	<- Z.prop
        N.curr	<- N.prop
      }
    }
    
    # update LOG_PI
    if(!mc$with.sampling | update.count==mc$sweep)
    {
      #	propose
      LOG_PI.prop		<- rdirichlet_log(Z.curr + mc$pars$LAMBDA[1,])
      #	this is the full conditional of LOG_PI given S, N, Z
      #	always accept
      #	update
      mc$curr.it	<- mc$curr.it+1L
      
      # update index for mc$it.info
      if(update.sweep.group==0L){
        mc$curr.it.adj <- mc$curr.it
      }else{
        mc$curr.it.adj <- mc$curr.it-mc$sweep_group*mc$sweep*update.sweep.group-1L
      }
      
      set(mc$it.info, mc$curr.it.adj, 'BLOCK', 'LOG_PI')
      set(mc$it.info, mc$curr.it.adj, 'PAR_ID', NA_integer_)
      set(mc$it.info, mc$curr.it.adj, 'MHRATIO', 1L)
      set(mc$it.info, mc$curr.it.adj, 'ACCEPT', 1L)
      LOG_PI.curr		<- LOG_PI.prop
      
      # at the end of sweep, record current parameters
      if(update.sweep.group==0L){
        update.round.adj <- update.round+1L
      }else{
        update.round.adj <- update.round+1L-mc$sweep_group*update.sweep.group-1L
      }
      
      mc$pars$XI[update.round.adj,]		<- XI.curr
      mc$pars$XI_LP[update.round.adj,]	<- XI_LP.curr
      mc$pars$S[update.round.adj,]		<- S.curr
      mc$pars$S_LP[update.round.adj,]	<- S_LP.curr
      mc$pars$Z[update.round.adj,]		<- Z.curr
      mc$pars$N[update.round.adj,]		<- N.curr
      mc$pars$LOG_PI[update.round.adj,]		<- LOG_PI.curr
    }
    
    # 	record log likelihood
    #   update index for mc$it.info
    if(update.sweep.group==0L){
      mc$curr.it.adj <- mc$curr.it
    }else{
      mc$curr.it.adj <- mc$curr.it-mc$sweep_group*mc$sweep*update.sweep.group-1L
    }
    
    tmp	<- sum(dbinom(dobs$TRM_OBS, size=Z.curr, prob=S.curr, log=TRUE) ) +
      dmultinom_log(Z.curr, size=N.curr, log_prob=LOG_PI.curr, log=TRUE)
    set(mc$it.info, mc$curr.it.adj, 'LOG_LKL', tmp)
    # 	record log prior
    tmp	<- dpois(N.curr, lambda=mc$pars$NU, log=TRUE) +
      lddirichlet_log_vector(LOG_PI.curr, nu=mc$pars$LAMBDA[1,]) +
      sum(S_LP.curr)
    set(mc$it.info, mc$curr.it.adj, 'LOG_PRIOR', tmp)
    
    if(update.count==mc$sweep & update.round %% 100 == 0){
      cat('\nSweeps done:\t',update.round)
    }
    
    if(update.count==mc$sweep & (update.round %% mc$sweep_group == 0|update.round == mc$nsweep)){
      
      mc$time	<- Sys.time()-ptm
      # save into RData
      if(update.sweep.group!=update.round %/% mc$sweep_group){
        save(object = mc, file = paste0(control$outfile,update.round %/% mc$sweep_group,".RData"))
      }else{
        # for the last groups of sweep if update.sweep.group==update.round %/% mc$sweep_group
        save(object = mc, file = paste0(control$outfile,update.round %/% mc$sweep_group + 1,".RData"))
      }
      # update the current group of sweeps
      update.sweep.group <- update.round %/% mc$sweep_group
      mc$pars$XI		<- matrix(NA_real_, ncol=length(update.info), nrow=mc$sweep_group) # prior for sampling in transmitter categories and sampling in recipient categories, concatenated
      mc$pars$XI_LP		<- matrix(NA_real_, ncol=length(update.info), nrow=mc$sweep_group) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated
      mc$pars$S			<- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group) # prior probability of sampling transmission pair categories
      mc$pars$S_LP		<- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group) # log prior density of sampling transmission pair categories
      mc$pars$Z			<- Matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group,sparse = TRUE) #augmented data
      mc$pars$N			<- matrix(NA_integer_, ncol=1, nrow=mc$sweep_group)							#total number of counts on augmented data
      mc$pars$LOG_PI		<- matrix(NA_real_, ncol=nrow(dobs), nrow=mc$sweep_group)	#proportions
      
      mc$it.info		<- data.table(	IT= seq.int(mc$sweep_group*mc$sweep*(update.round%/%mc$sweep_group)+1,mc$sweep_group*mc$sweep*(1+(update.round%/%mc$sweep_group))),
                                 PAR_ID= rep(NA_integer_, mc$sweep_group*mc$sweep),
                                 BLOCK= rep(NA_character_, mc$sweep_group*mc$sweep),
                                 MHRATIO= rep(NA_real_, mc$sweep_group*mc$sweep),
                                 ACCEPT=rep(NA_integer_, mc$sweep_group*mc$sweep),
                                 LOG_LKL=rep(NA_real_, mc$sweep_group*mc$sweep),
                                 LOG_PRIOR=rep(NA_real_, mc$sweep_group*mc$sweep)
      )
      
    }
  }
  
  # if(!'outfile'%in%names(control))
  #   return(mc)
  # save(mc,	file=control$outfile)
  NULL
}
