# log product of poisson densities with parameter lambda_ab
lpois_prod <- function(x, llbd){
    # deal with small lambda whose exp is 0
    c <- min((max(llbd)-min(llbd))/2,700)
    ans <- - exp(-c) * sum(exp(llbd+c)) + sum(x * llbd) - sum(lgamma(x+1))
    ans
}

# log product of gamma densities with parameter alpha_{ab} and beta
lgamma_prod <- function(llbd, a, b){
    # deal with small lambda whose exp is 0
    c <- min((max(llbd)-min(llbd))/2,700)
    ans <- sum( a * log(b) ) - sum( lgamma(a) ) +
    sum( (a-1) * llbd ) - b * exp(-c) * sum(exp(llbd + c))
    ans
}

# log product of gamma densities with parameter alpha_{ab} and beta_{ab}
lgamma_prod_vb <- function(llbd, a, b){
  # deal with small lambda whose exp is 0
  c <- min(
    (max(llbd + log(b)) - min(llbd + log(b)))/2,
    700)
  ans <- sum( a * log(b) ) - sum( lgamma(a) ) +
    sum( (a-1) * llbd ) - exp(-c) * sum(exp(llbd + log(b) + c))
  ans
}

# talen from
# Goal: Sample small-shape gamma random variables via accept-reject
# Paper: arXiv:1302.1884
rgamss <- function(n, shape, scale=1, do.log=TRUE) 
{  
	a <- shape
	if(a > 0.2)
	{
		oo <- ifelse(do.log,log(rgamma(n, a, 1 / scale)),rgamma(n, a, 1 / scale))
	}
	else 
	{
		e1 <- 2.71828182845905
		L <- 1 / a - 1
		w <- a / e1 / (1 - a)
		ww <- 1 / (1 + w)
		eta <- function(z)
		{ 
			ifelse(z >= 0, exp(-z), w * L * exp(L * z)) 
		}
		h <- function(z)
		{
			exp(-z - exp(-z / a))
		} 
		rh <- function(a)
		{            
			repeat 
			{        
				U <- runif(1)
				z <- ifelse(U <= ww, -log(U / ww), log(runif(1)) / L)
				if(h(z) / eta(z) > runif(1)) 
					return(z)        
			}
		}
		Z <- numeric(n)
		for(i in 1:n) 
			Z[i] <- rh(a)
		o <- log(scale) - Z / a
		if(!do.log) 
		{      
			oo <- exp(o)
			if(any(oo == 0)) 
			{        
				oo <- o
				warning("Output given on log-scale since shape is small")        
			}      
		} 
		else 
			oo <- o    
	}
	return(oo)
}

#' @title Source attribution while adjusting for sampling bias
#' @export
#' @import data.table
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
source.attribution.mcmc  <- function(dobs, dprior, control=list(seed=42, mcmc.n=1e3, verbose=1, outfile='SAMCMCv190327.rda'))
{        
    dmode <- function(x) 
	{
        den <- density(x, kernel=c("gaussian"))
        (den$x[den$y==max(den$y)])
    }
    # basic checks on inputs
    if(!all(dobs$TR_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
    	stop('Did not find prior samples of a sampling category for TR_SAMPLING_CATEGORY')
    if(!all(dobs$REC_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
    	stop('Did not find prior samples of a sampling category for REC_SAMPLING_CATEGORY')
	
    # keep track of runtime
    ptm    <- Sys.time()
	
	# set seed if provided
    if('seed'%in%names(control))
    {
        cat('\nSetting seed to',control$seed)
        set.seed(control$seed)
    }    
	
    # set up mcmc
    mc <- list()
    # determine if the sampling probabilities are <1.
    # If they are NOT, Z will be the same as TRM_OBS, and the algorithm only updates LAMBDA
    mc$with.sampling <- dprior[, list(ALL_ONE=all(P==1)), by='SAMPLING_CATEGORY'][, !all(ALL_ONE)]
    mc$time <- NA_real_
	
    # for speed: make look-up table so we know which transmission pair categories need to be updated
    # at every MCMC iteration
    
    # set UPDATE_ID for each sampling category
    mc$dlu <- data.table(
      SAMPLING_CATEGORY = unique(c(dobs$TR_SAMPLING_CATEGORY,dobs$REC_SAMPLING_CATEGORY)))
    setkey(mc$dlu,SAMPLING_CATEGORY)
    mc$dlu[, UPDATE_ID:= seq_len(nrow(mc$dlu))]
    setkey(mc$dlu, UPDATE_ID)
    
    # make data.table that maps TRM_CAT_PAIR_IDs to transmitter UPDATE_IDs and recipient UPDATE_IDs
    mc$dlt <- subset(dobs, select = c('TRM_CAT_PAIR_ID','TR_SAMPLING_CATEGORY',
                                      'REC_SAMPLING_CATEGORY','TRM_OBS'))
    setnames(mc$dlu, colnames(mc$dlu), paste0('TR_',colnames(mc$dlu)))
    mc$dlt <- merge(mc$dlt, mc$dlu, by='TR_SAMPLING_CATEGORY')
    setnames(mc$dlu, colnames(mc$dlu), gsub('TR_','REC_',colnames(mc$dlu)))
    mc$dlt <- merge(mc$dlt, mc$dlu, by='REC_SAMPLING_CATEGORY')
    setnames(mc$dlu, colnames(mc$dlu), gsub('REC_','',colnames(mc$dlu)))
    set(mc$dlt,NULL,c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY'),NULL)
    setkey(mc$dlt,TRM_CAT_PAIR_ID)
    
	
    # for speed: make indexed lookup table 
    update.info    <- vector('list', nrow(mc$dlu))
    for (i in 1:nrow(mc$dlu))
    {
      update.info[[i]] <- sort( unique(c(mc$dlt[TR_UPDATE_ID==i,TRM_CAT_PAIR_ID],
                                         mc$dlt[REC_UPDATE_ID==i,TRM_CAT_PAIR_ID])) )
    }
    # for speed: make indexed prior samples
    setkey(dprior, SAMPLING_CATEGORY)
    dprior2 <- vector('list', nrow(mc$dlu))
    for (i in 1:nrow(mc$dlu))
    {
        tmp <- mc$dlu[UPDATE_ID==i,SAMPLING_CATEGORY]
        dprior2[[i]] <- dprior[J(tmp),nomatch=0L]
    }
    
    # estimate the mean sampling probabilities for each transmission group
    if(dprior[, max(SAMPLE)]>10)
    {
        dprior3 <- dprior[, list(EST_SAMPLING_RATE=dmode(P)), by=SAMPLING_CATEGORY]
    }
    if(dprior[, max(SAMPLE)]<=10)
    {
        dprior3 <- dprior[, list(EST_SAMPLING_RATE=median(P)), by=SAMPLING_CATEGORY]
    }
	
	# add estimated sampling probabilities to dobs
    dobs2 <- subset(dobs, select = c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY','TRM_CAT_PAIR_ID'))
    setnames(dprior3, colnames(dprior3), paste0('TR_',colnames(dprior3)))
    dobs2 <- merge(dobs2, dprior3, by='TR_SAMPLING_CATEGORY')
    setnames(dprior3, colnames(dprior3), gsub('TR_','REC_', colnames(dprior3)))
    dobs2 <- merge(dobs2,dprior3,by='REC_SAMPLING_CATEGORY')
    dobs2[, EST_SAMPLING_RATE:= TR_EST_SAMPLING_RATE * REC_EST_SAMPLING_RATE]
    
    mc$nprior <- max(dprior$SAMPLE)
    mc$sweep <- nrow(mc$dlu)*2L
    mc$nsweep <- ceiling( control$mcmc.n/mc$sweep )
    mc$n <- mc$nsweep*mc$sweep
    # save and reset mc$pars and mc$it.info to file every mc$sweep_group sweeps
	mc$sweep_group <- mc$nsweep
	if('sweep_group'%in%names(control))
	{
		mc$sweep_group <- control$sweep_group	
	}    
    mc$pars <- list()
    mc$pars$ALPHA <- matrix(NA_real_, ncol= nrow(dobs), nrow= 1L)
    mc$pars$BETA <- NA_real_
    mc$pars$XI <- matrix(NA_real_, ncol=length(update.info), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # prior for sampling in categories
    mc$pars$XI_LP <- matrix(NA_real_, ncol=length(update.info), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated
    mc$pars$S <- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # prior probability of sampling transmission pair categories
    mc$pars$S_LP <- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # log prior density of sampling transmission pair categories
    mc$pars$LOG_LAMBDA <- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # mean number of counts on augmented data for transmission pair categories
    
	mc$it.info <- data.table(	IT= seq.int(0,min(mc$nsweep*mc$sweep,mc$sweep_group*mc$sweep)),
   								PAR_ID= rep(NA_integer_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
    							BLOCK= rep(NA_character_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
    							MHRATIO= rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
    							ACCEPT=rep(NA_integer_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
    							LOG_LKL=rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
    							LOG_PRIOR=rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L))
								)
    if(1)
    {
        cat('\nNumber of parameters:\t', ncol(mc$pars$LOG_LAMBDA)+ncol(mc$pars$XI) )
        cat('\nDimension of PI:\t', ncol(mc$pars$LOG_LAMBDA))
        cat('\nSweep length:\t', mc$sweep)
        cat('\nNumber of sweeps:\t', mc$nsweep)
        cat('\nNumber of iterations:\t', mc$n)        
    }
    
    
    #
    # initialise MCMC
    #
    update.sweep.group <- 0L
    mc$curr.it <- 1L
    set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
    set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
    set(mc$it.info, mc$curr.it, 'MHRATIO', 1)
    set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)
	
    # initialise prior on lambda (ALPHA, BETA): 
	# specify prior so that we obtain for PI the Berger objective prior with minimal loss compared to marginal Beta reference prior
    # (https://projecteuclid.org/euclid.ba/1422556416)    
    mc$pars$ALPHA[1,]    <- 0.8/nrow(dobs)
    dobs3 <- merge(	subset(dobs,select = c('TRM_OBS','TRM_CAT_PAIR_ID')),
    				subset(dobs2,select = c('EST_SAMPLING_RATE','TRM_CAT_PAIR_ID')),
    				by=c('TRM_CAT_PAIR_ID'))
    dobs3[, LBD:=TRM_OBS/EST_SAMPLING_RATE]
    dobs3[TRM_OBS==0, LBD:=(1-EST_SAMPLING_RATE)/EST_SAMPLING_RATE]    
    mc$pars$BETA <- 0.8/sum(dobs3$LBD)
    
    # initialise sampling probabilities for pop strata: random draw from input sampling probabilities
    tmp <- subset(dprior, SAMPLE==sample(mc$nprior,1))
    tmp <- merge(mc$dlu, tmp, by='SAMPLING_CATEGORY')
    setkey(tmp, UPDATE_ID)
    mc$pars$XI[1, tmp$UPDATE_ID] <- tmp$P
    mc$pars$XI_LP[1, tmp$UPDATE_ID] <- tmp$LP
    
    # initialise sampling probabilities for transmission groups: multiply
    mc$pars$S[1,] <- mc$pars$XI[1, mc$dlt$TR_UPDATE_ID] * mc$pars$XI[1, mc$dlt$REC_UPDATE_ID]
    mc$pars$S_LP[1,] <- mc$pars$XI_LP[1, mc$dlt$TR_UPDATE_ID] + mc$pars$XI_LP[1, mc$dlt$REC_UPDATE_ID]
    
	# initialise log_lambda: random draw from full conditional
    mc$pars$LOG_LAMBDA[1,] <- mapply( function(x,y){ rgamss(1,shape=x,scale=y) }, dobs$TRM_OBS+mc$pars$ALPHA, 1/(mc$pars$S[1,]+mc$pars$BETA)	)
	
    #    store log likelihood
    tmp <- lpois_prod( dobs$TRM_OBS, mc$pars$LOG_LAMBDA[1,] + log(mc$pars$S[1,]))
    set(mc$it.info, 1L, 'LOG_LKL', tmp)
	
    #     store log prior
    tmp <- lgamma_prod( llbd= mc$pars$LOG_LAMBDA[1,], a = mc$pars$ALPHA , b = mc$pars$BETA) +
    sum(mc$pars$XI_LP[1,])
    set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
    
    
    # parameter value at the current step
    XI.curr <- mc$pars$XI[1,]
    XI_LP.curr <- mc$pars$XI_LP[1,]
    S.curr <- mc$pars$S[1,]
    S_LP.curr <- mc$pars$S_LP[1,]
    LOG_LAMBDA.curr <- mc$pars$LOG_LAMBDA[1,]
    
    # run mcmc
    options(warn=0)
    for(i in 1L:mc$n)
    {
        mc$curr.it <- i
        # determine the current iteration in a sweep
        update.count <- (i-1L) %% mc$sweep + 1L
        update.round <- (i-1L) %/% mc$sweep + 1L
        # update XI first
        if(mc$with.sampling & (update.count %% 2==1))
        {
            # determine the category to update
            update.group <-  update.count %/% 2 + 1
            update.cat <- mc$dlu$SAMPLING_CATEGORY[update.group]
            
            # choose the pairs to update
            update.pairs <- update.info[[update.group]]
            
            # propose single XI
            XI.prop <- XI.curr
            XI_LP.prop <- XI_LP.curr
            tmp <- dprior2[[update.group]][sample(mc$nprior,1),]
            if(tmp$SAMPLING_CATEGORY[1]!=update.cat)
            	stop('\nFatal error in dprior2.')
            XI.prop[ update.group ] <- tmp$P
            XI_LP.prop[ update.group ] <- tmp$LP
            
            # propose all S that involve the one XI from above
            S.prop <- S.curr
            S_LP.prop <- S_LP.curr
            S.prop[update.pairs] <- XI.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] * XI.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
            S_LP.prop[update.pairs] <- XI_LP.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] + XI_LP.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
            # calculate MH ratio
            log.fc <- lpois_prod( mc$dlt$TRM_OBS[update.pairs], LOG_LAMBDA.curr[update.pairs] + log(S.curr[update.pairs]))
            log.fc.prop <- lpois_prod( mc$dlt$TRM_OBS[update.pairs],  LOG_LAMBDA.curr[update.pairs] + log(S.prop[update.pairs]))
            
            log.mh.ratio <- log.fc.prop - log.fc
            mh.ratio <- min(1,exp(log.mh.ratio))
            #    update
            mc$curr.it <- mc$curr.it+1L
            
            if(update.sweep.group==0L)
			{
                mc$curr.it.adj <- mc$curr.it
            }
			if(update.sweep.group!=0L)
			{
                mc$curr.it.adj <- mc$curr.it - mc$sweep_group*mc$sweep*update.sweep.group - 1L
            }            
            set(mc$it.info, mc$curr.it.adj, 'BLOCK', 'XI')
            set(mc$it.info, mc$curr.it.adj, 'PAR_ID', update.group)
            set(mc$it.info, mc$curr.it.adj, 'MHRATIO', mh.ratio)
            accept <- as.integer(runif(1) < mh.ratio)
            set(mc$it.info, mc$curr.it.adj, 'ACCEPT', accept)
            
            if(control$verbose & accept)
            {
                cat('\nit ',mc$curr.it,' ACCEPT XI ',update.group)
            }
            if(accept)
            {
                XI.curr <- XI.prop
                XI_LP.curr <- XI_LP.prop
                S.curr <- S.prop
                S_LP.curr <- S_LP.prop
            }            
        }        
        if(mc$with.sampling & (update.count %% 2==0))
        {            
            # propose
            LOG_LAMBDA.prop <- LOG_LAMBDA.curr
            LOG_LAMBDA.prop[update.pairs] <- mapply( function(x,y){ rgamss(1,shape=x,scale=y) }, mc$dlt$TRM_OBS[update.pairs] + mc$pars$ALPHA[update.pairs], 1/(S.curr[update.pairs] + mc$pars$BETA))
                        
            #    this is the full conditional of LOG_LAMBDA given S
            #    always accept
            #    update            
            mc$curr.it <- mc$curr.it+1L            
            if(update.sweep.group==0L)
			{
                mc$curr.it.adj <- mc$curr.it
            }
			if(update.sweep.group!=0L)
			{
                mc$curr.it.adj <- as.integer(mc$curr.it - mc$sweep_group*mc$sweep*update.sweep.group - 1L)
            }            
            set(mc$it.info, mc$curr.it.adj, 'BLOCK', 'LOG_LAMBDA')
            set(mc$it.info, mc$curr.it.adj, 'PAR_ID', update.group)
            set(mc$it.info, mc$curr.it.adj, 'MHRATIO', 1L)
            set(mc$it.info, mc$curr.it.adj, 'ACCEPT', 1L)            
            LOG_LAMBDA.curr        <- LOG_LAMBDA.prop            
        }
                
        #     record log likelihood
        if(update.sweep.group==0L)
		{
            mc$curr.it.adj <- mc$curr.it
		}
		if(update.sweep.group!=0L)
		{
            mc$curr.it.adj <- as.integer(mc$curr.it - mc$sweep_group*mc$sweep*update.sweep.group - 1L)
        }        
        tmp <- lpois_prod( dobs$TRM_OBS, LOG_LAMBDA.curr + log(S.curr))
        set(mc$it.info, mc$curr.it.adj, 'LOG_LKL', tmp)
        #     record log prior
        tmp <- lgamma_prod( llbd= LOG_LAMBDA.curr, a = mc$pars$ALPHA , b = mc$pars$BETA) + sum(XI_LP.curr)
        set(mc$it.info, mc$curr.it.adj, 'LOG_PRIOR', tmp)

		if(update.count == mc$sweep)
		{
            # at the end of sweep, record current parameters
            if(update.sweep.group==0L)
			{
                update.round.adj <- update.round + 1L
			}
			if(update.sweep.group!=0L)
			{
                update.round.adj <- update.round + 1L - mc$sweep_group*update.sweep.group - 1L
            }            
            mc$pars$XI[update.round.adj,] <- XI.curr
            mc$pars$XI_LP[update.round.adj,] <- XI_LP.curr
            mc$pars$S[update.round.adj,] <- S.curr
            mc$pars$S_LP[update.round.adj,] <- S_LP.curr
            mc$pars$LOG_LAMBDA[update.round.adj,] <- LOG_LAMBDA.curr
        }

        if(update.count==mc$sweep & (update.round %% 100 == 0 | update.round == mc$nsweep))
		{
            cat('\nSweeps done:\t',update.round)
        }
        
        
        if(update.count==mc$sweep & (update.round %% mc$sweep_group == 0 | update.round == mc$nsweep))
		{
            mc$time    <- Sys.time()-ptm
            # save
            if(!'outfile'%in%names(control))
            	return(mc)
            if(update.sweep.group!=update.round %/% mc$sweep_group)
			{
                save(object = mc, file = paste0(control$outfile,'_sweepgrp',update.round %/% mc$sweep_group,".rda"))
            }
			if(update.sweep.group==update.round %/% mc$sweep_group)
			{
                save(object = mc, file = paste0(control$outfile,'_sweepgrp',update.round %/% mc$sweep_group + 1,".rda"))
            }
            update.sweep.group <- update.round %/% mc$sweep_group
            
            # reset mc$pars and mc$it.info
            mc$pars$XI <- matrix(NA_real_, ncol=length(update.info), nrow=mc$sweep_group) # prior for sampling in transmitter categories and sampling in recipient categories, concatenated
            mc$pars$XI_LP <- matrix(NA_real_, ncol=length(update.info), nrow=mc$sweep_group) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated
            mc$pars$S <- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group) # prior probability of sampling transmission pair categories
            mc$pars$S_LP <- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group) # log prior density of sampling transmission pair categories
            mc$pars$LOG_LAMBDA <- matrix(NA_real_, ncol=nrow(dobs), nrow=mc$sweep_group)    #proportions
            
            mc$it.info <- data.table(	IT= seq.int(mc$sweep_group*mc$sweep*(update.round%/%mc$sweep_group)+1,mc$sweep_group*mc$sweep*(1+(update.round%/%mc$sweep_group))),
            							PAR_ID= rep(NA_integer_, mc$sweep_group*mc$sweep),
            							BLOCK= rep(NA_character_, mc$sweep_group*mc$sweep),
            							MHRATIO= rep(NA_real_, mc$sweep_group*mc$sweep),
            							ACCEPT=rep(NA_integer_, mc$sweep_group*mc$sweep),
            							LOG_LKL=rep(NA_real_, mc$sweep_group*mc$sweep),
            							LOG_PRIOR=rep(NA_real_, mc$sweep_group*mc$sweep)
            							)            
        }                
    }
}



source.attribution.mcmc.jointupdate  <- function(dobs, dprior, control=list(seed=42, mcmc.n=1e3, verbose=1, outfile='SAMCMCv190327.rda'))
{        
  dmode <- function(x) 
  {
    den <- density(x, kernel=c("gaussian"))
    (den$x[den$y==max(den$y)])
  }
  # basic checks on inputs
  if(!all(dobs$TR_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
    stop('Did not find prior samples of a sampling category for TR_SAMPLING_CATEGORY')
  if(!all(dobs$REC_SAMPLING_CATEGORY %in% dprior$SAMPLING_CATEGORY))
    stop('Did not find prior samples of a sampling category for REC_SAMPLING_CATEGORY')
  
  # keep track of runtime
  ptm    <- Sys.time()
  
  # set seed if provided
  if('seed'%in%names(control))
  {
    cat('\nSetting seed to',control$seed)
    set.seed(control$seed)
  }    
  
  # set up mcmc
  mc <- list()
  # determine if the sampling probabilities are <1.
  # If they are NOT, Z will be the same as TRM_OBS, and the algorithm only updates LAMBDA
  mc$with.sampling <- dprior[, list(ALL_ONE=all(P==1)), by='SAMPLING_CATEGORY'][, !all(ALL_ONE)]
  mc$time <- NA_real_
  
  # for speed: make look-up table so we know which transmission pair categories need to be updated
  # at every MCMC iteration
  
  # set UPDATE_ID for each sampling category
  mc$dlu <- data.table(
    SAMPLING_CATEGORY = unique(c(dobs$TR_SAMPLING_CATEGORY,dobs$REC_SAMPLING_CATEGORY)))
  setkey(mc$dlu,SAMPLING_CATEGORY)
  mc$dlu[, UPDATE_ID:= seq_len(nrow(mc$dlu))]
  setkey(mc$dlu, UPDATE_ID)
  
  # make data.table that maps TRM_CAT_PAIR_IDs to transmitter UPDATE_IDs and recipient UPDATE_IDs
  mc$dlt <- subset(dobs, select = c('TRM_CAT_PAIR_ID','TR_SAMPLING_CATEGORY',
                                    'REC_SAMPLING_CATEGORY','TRM_OBS'))
  setnames(mc$dlu, colnames(mc$dlu), paste0('TR_',colnames(mc$dlu)))
  mc$dlt <- merge(mc$dlt, mc$dlu, by='TR_SAMPLING_CATEGORY')
  setnames(mc$dlu, colnames(mc$dlu), gsub('TR_','REC_',colnames(mc$dlu)))
  mc$dlt <- merge(mc$dlt, mc$dlu, by='REC_SAMPLING_CATEGORY')
  setnames(mc$dlu, colnames(mc$dlu), gsub('REC_','',colnames(mc$dlu)))
  set(mc$dlt,NULL,c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY'),NULL)
  setkey(mc$dlt,TRM_CAT_PAIR_ID)
  
  
  # for speed: make indexed lookup table 
  update.info    <- vector('list', nrow(mc$dlu))
  for (i in 1:nrow(mc$dlu))
  {
    update.info[[i]] <- sort( unique(c(mc$dlt[TR_UPDATE_ID==i,TRM_CAT_PAIR_ID],
                                       mc$dlt[REC_UPDATE_ID==i,TRM_CAT_PAIR_ID])) )
  }
  # for speed: make indexed prior samples
  setkey(dprior, SAMPLING_CATEGORY)
  dprior2 <- vector('list', nrow(mc$dlu))
  for (i in 1:nrow(mc$dlu))
  {
    tmp <- mc$dlu[UPDATE_ID==i,SAMPLING_CATEGORY]
    dprior2[[i]] <- dprior[J(tmp),nomatch=0L]
  }
  
  # estimate the mean sampling probabilities for each transmission group
  if(dprior[, max(SAMPLE)]>10)
  {
    dprior3 <- dprior[, list(EST_SAMPLING_RATE=dmode(P)), by=SAMPLING_CATEGORY]
  }
  if(dprior[, max(SAMPLE)]<=10)
  {
    dprior3 <- dprior[, list(EST_SAMPLING_RATE=median(P)), by=SAMPLING_CATEGORY]
  }
  
  # add estimated sampling probabilities to dobs
  dobs2 <- subset(dobs, select = c('TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY','TRM_CAT_PAIR_ID'))
  setnames(dprior3, colnames(dprior3), paste0('TR_',colnames(dprior3)))
  dobs2 <- merge(dobs2, dprior3, by='TR_SAMPLING_CATEGORY')
  setnames(dprior3, colnames(dprior3), gsub('TR_','REC_', colnames(dprior3)))
  dobs2 <- merge(dobs2,dprior3,by='REC_SAMPLING_CATEGORY')
  dobs2[, EST_SAMPLING_RATE:= TR_EST_SAMPLING_RATE * REC_EST_SAMPLING_RATE]
  
  mc$nprior <- max(dprior$SAMPLE)
  mc$sweep <- nrow(mc$dlu)
  mc$nsweep <- ceiling( control$mcmc.n/mc$sweep )
  mc$n <- mc$nsweep*mc$sweep
  # save and reset mc$pars and mc$it.info to file every mc$sweep_group sweeps
  mc$sweep_group <- mc$nsweep
  if('sweep_group'%in%names(control))
  {
    mc$sweep_group <- control$sweep_group	
  }    
  mc$pars <- list()
  mc$pars$ALPHA <- matrix(NA_real_, ncol= nrow(dobs), nrow= 1L)
  mc$pars$BETA <- NA_real_
  mc$pars$XI <- matrix(NA_real_, ncol=length(update.info), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # prior for sampling in categories
  mc$pars$XI_LP <- matrix(NA_real_, ncol=length(update.info), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated
  mc$pars$S <- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # prior probability of sampling transmission pair categories
  mc$pars$S_LP <- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # log prior density of sampling transmission pair categories
  mc$pars$LOG_LAMBDA <- matrix(NA_integer_, ncol=nrow(dobs), nrow=min(mc$nsweep+1L,mc$sweep_group+1L)) # mean number of counts on augmented data for transmission pair categories
  
  mc$it.info <- data.table(	IT= seq.int(0,min(mc$nsweep*mc$sweep,mc$sweep_group*mc$sweep)),
                            PAR_ID= rep(NA_integer_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                            BLOCK= rep(NA_character_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                            MHRATIO= rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                            ACCEPT=rep(NA_integer_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                            LOG_LKL=rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L)),
                            LOG_PRIOR=rep(NA_real_, min(mc$nsweep*mc$sweep+1L,mc$sweep_group*mc$sweep+1L))
  )
  if(1)
  {
    cat('\nNumber of parameters:\t', ncol(mc$pars$LOG_LAMBDA)+ncol(mc$pars$XI) )
    cat('\nDimension of PI:\t', ncol(mc$pars$LOG_LAMBDA))
    cat('\nSweep length:\t', mc$sweep)
    cat('\nNumber of sweeps:\t', mc$nsweep)
    cat('\nNumber of iterations:\t', mc$n)        
  }
  
  
  #
  # initialise MCMC
  #
  update.sweep.group <- 0L
  mc$curr.it <- 1L
  set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
  set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
  set(mc$it.info, mc$curr.it, 'MHRATIO', 1)
  set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)
  
  # initialise prior on lambda (ALPHA, BETA): 
  # specify prior so that we obtain for PI the Berger objective prior with minimal loss compared to marginal Beta reference prior
  # (https://projecteuclid.org/euclid.ba/1422556416)    
  mc$pars$ALPHA[1,]    <- 0.8/nrow(dobs)
  dobs3 <- merge(	subset(dobs,select = c('TRM_OBS','TRM_CAT_PAIR_ID')),
                  subset(dobs2,select = c('EST_SAMPLING_RATE','TRM_CAT_PAIR_ID')),
                  by=c('TRM_CAT_PAIR_ID'))
  dobs3[, LBD:=TRM_OBS/EST_SAMPLING_RATE]
  dobs3[TRM_OBS==0, LBD:=(1-EST_SAMPLING_RATE)/EST_SAMPLING_RATE]    
  mc$pars$BETA <- 0.8/sum(dobs3$LBD)
  
  # initialise sampling probabilities for pop strata: random draw from input sampling probabilities
  tmp <- subset(dprior, SAMPLE==sample(mc$nprior,1))
  tmp <- merge(mc$dlu, tmp, by='SAMPLING_CATEGORY')
  setkey(tmp, UPDATE_ID)
  mc$pars$XI[1, tmp$UPDATE_ID] <- tmp$P
  mc$pars$XI_LP[1, tmp$UPDATE_ID] <- tmp$LP
  
  # initialise sampling probabilities for transmission groups: multiply
  mc$pars$S[1,] <- mc$pars$XI[1, mc$dlt$TR_UPDATE_ID] * mc$pars$XI[1, mc$dlt$REC_UPDATE_ID]
  mc$pars$S_LP[1,] <- mc$pars$XI_LP[1, mc$dlt$TR_UPDATE_ID] + mc$pars$XI_LP[1, mc$dlt$REC_UPDATE_ID]
  
  # initialise log_lambda: random draw from full conditional
  mc$pars$LOG_LAMBDA[1,] <- mapply( function(x,y){ rgamss(1,shape=x,scale=y) }, dobs$TRM_OBS+mc$pars$ALPHA, 1/(mc$pars$S[1,]+mc$pars$BETA)	)
  
  #    store log likelihood
  tmp <- lpois_prod( dobs$TRM_OBS, mc$pars$LOG_LAMBDA[1,] + log(mc$pars$S[1,]))
  set(mc$it.info, 1L, 'LOG_LKL', tmp)
  
  #     store log prior
  tmp <- lgamma_prod( llbd= mc$pars$LOG_LAMBDA[1,], a = mc$pars$ALPHA , b = mc$pars$BETA) +
    sum(mc$pars$XI_LP[1,])
  set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
  
  
  # parameter value at the current step
  XI.curr <- mc$pars$XI[1,]
  XI_LP.curr <- mc$pars$XI_LP[1,]
  S.curr <- mc$pars$S[1,]
  S_LP.curr <- mc$pars$S_LP[1,]
  LOG_LAMBDA.curr <- mc$pars$LOG_LAMBDA[1,]
  
  # run mcmc
  options(warn=0)
  for(i in 1L:mc$n)
  {
    mc$curr.it <- i
    # determine the current iteration in a sweep
    update.count <- (i-1L) %% mc$sweep + 1L
    update.round <- (i-1L) %/% mc$sweep + 1L

    if(mc$with.sampling)
    {
      # determine the category to update
      update.cat <- mc$dlu$SAMPLING_CATEGORY[update.count]
      
      # choose the pairs to update
      update.pairs <- update.info[[update.count]]
      
      # propose single XI
      XI.prop <- XI.curr
      XI_LP.prop <- XI_LP.curr
      tmp <- dprior2[[update.count]][sample(mc$nprior,1),]
      if(tmp$SAMPLING_CATEGORY[1]!=update.cat)
        stop('\nFatal error in dprior2.')
      XI.prop[ update.count ] <- tmp$P
      XI_LP.prop[ update.count ] <- tmp$LP
      
      # propose all S that involve the one XI from above
      S.prop <- S.curr
      S_LP.prop <- S_LP.curr
      S.prop[update.pairs] <- XI.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] * XI.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
      S_LP.prop[update.pairs] <- XI_LP.prop[ mc$dlt$TR_UPDATE_ID[update.pairs] ] + XI_LP.prop[ mc$dlt$REC_UPDATE_ID[update.pairs] ]
      
      # propose lambda
      LOG_LAMBDA.prop <- LOG_LAMBDA.curr
      LOG_LAMBDA.prop[update.pairs] <- mapply( function(x,y){ rgamss(1,shape=x,scale=y) }, mc$dlt$TRM_OBS[update.pairs] + mc$pars$ALPHA[update.pairs], 1/(S.prop[update.pairs] + mc$pars$BETA))
      
      
      # calculate MH ratio
      log.fc <- lpois_prod( mc$dlt$TRM_OBS[update.pairs], 
                            LOG_LAMBDA.curr[update.pairs] + log(S.curr[update.pairs]))+
        lgamma_prod(LOG_LAMBDA.curr[update.pairs], 
                    mc$pars$ALPHA[update.pairs],
                    mc$pars$BETA) 
        
      
      log.fc.prop <- lpois_prod( mc$dlt$TRM_OBS[update.pairs],  
                                 LOG_LAMBDA.prop[update.pairs] + log(S.prop[update.pairs]))+
        lgamma_prod(LOG_LAMBDA.prop[update.pairs], 
                    mc$pars$ALPHA[update.pairs],
                    mc$pars$BETA) 
      
      log.prop.ratio <- lgamma_prod_vb(LOG_LAMBDA.curr[update.pairs],
                                       mc$dlt$TRM_OBS[update.pairs] + mc$pars$ALPHA[update.pairs],
                                       mc$pars$BETA + S.curr[update.pairs])-
        lgamma_prod_vb(LOG_LAMBDA.prop[update.pairs],
                       mc$dlt$TRM_OBS[update.pairs] + mc$pars$ALPHA[update.pairs],
                       mc$pars$BETA + S.prop[update.pairs])
      
      log.mh.ratio <- log.fc.prop - log.fc + log.prop.ratio
      mh.ratio <- min(1,exp(log.mh.ratio))
      #    update
      mc$curr.it <- mc$curr.it+1L
      
      if(update.sweep.group==0L)
      {
        mc$curr.it.adj <- mc$curr.it
      }
      if(update.sweep.group!=0L)
      {
        mc$curr.it.adj <- mc$curr.it - mc$sweep_group*mc$sweep*update.sweep.group - 1L
      }            
      set(mc$it.info, mc$curr.it.adj, 'BLOCK', 'XI-LAMBDA')
      set(mc$it.info, mc$curr.it.adj, 'PAR_ID', update.count)
      set(mc$it.info, mc$curr.it.adj, 'MHRATIO', mh.ratio)
      accept <- as.integer(runif(1) < mh.ratio)
      set(mc$it.info, mc$curr.it.adj, 'ACCEPT', accept)
      
      if(control$verbose & accept)
      {
        cat('\nit ',mc$curr.it,' ACCEPT XI-LAMBDA ',update.count)
      }
      if(accept)
      {
        XI.curr <- XI.prop
        XI_LP.curr <- XI_LP.prop
        S.curr <- S.prop
        S_LP.curr <- S_LP.prop
        LOG_LAMBDA.curr        <- LOG_LAMBDA.prop  
      }            
    }        

    #     record log likelihood
    if(update.sweep.group==0L)
    {
      mc$curr.it.adj <- mc$curr.it
    }
    if(update.sweep.group!=0L)
    {
      mc$curr.it.adj <- as.integer(mc$curr.it - mc$sweep_group*mc$sweep*update.sweep.group - 1L)
    }        
    tmp <- lpois_prod( dobs$TRM_OBS, LOG_LAMBDA.curr + log(S.curr))
    set(mc$it.info, mc$curr.it.adj, 'LOG_LKL', tmp)
    #     record log prior
    tmp <- lgamma_prod( llbd= LOG_LAMBDA.curr, a = mc$pars$ALPHA , b = mc$pars$BETA) + 
      sum(XI_LP.curr)
    set(mc$it.info, mc$curr.it.adj, 'LOG_PRIOR', tmp)
    
    if(update.count == mc$sweep)
    {
      # at the end of sweep, record current parameters
      if(update.sweep.group==0L)
      {
        update.round.adj <- update.round + 1L
      }
      if(update.sweep.group!=0L)
      {
        update.round.adj <- update.round + 1L - mc$sweep_group*update.sweep.group - 1L
      }            
      mc$pars$XI[update.round.adj,] <- XI.curr
      mc$pars$XI_LP[update.round.adj,] <- XI_LP.curr
      mc$pars$S[update.round.adj,] <- S.curr
      mc$pars$S_LP[update.round.adj,] <- S_LP.curr
      mc$pars$LOG_LAMBDA[update.round.adj,] <- LOG_LAMBDA.curr
    }
    
    if(update.count==mc$sweep & (update.round %% 100 == 0 | update.round == mc$nsweep))
    {
      cat('\nSweeps done:\t',update.round)
    }
    
    
    if(update.count==mc$sweep & (update.round %% mc$sweep_group == 0 | update.round == mc$nsweep))
    {
      mc$time    <- Sys.time()-ptm
      # save
      if(!'outfile'%in%names(control))
        return(mc)
      if(update.sweep.group!=update.round %/% mc$sweep_group)
      {
        save(object = mc, file = paste0(control$outfile,'_sweepgrp',update.round %/% mc$sweep_group,".rda"))
      }
      if(update.sweep.group==update.round %/% mc$sweep_group)
      {
        save(object = mc, file = paste0(control$outfile,'_sweepgrp',update.round %/% mc$sweep_group + 1,".rda"))
      }
      update.sweep.group <- update.round %/% mc$sweep_group
      
      # reset mc$pars and mc$it.info
      mc$pars$XI <- matrix(NA_real_, ncol=length(update.info), nrow=mc$sweep_group) # prior for sampling in transmitter categories and sampling in recipient categories, concatenated
      mc$pars$XI_LP <- matrix(NA_real_, ncol=length(update.info), nrow=mc$sweep_group) # log prior density for sampling in transmitter categories and sampling in recipient categories, concatenated
      mc$pars$S <- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group) # prior probability of sampling transmission pair categories
      mc$pars$S_LP <- matrix(NA_integer_, ncol=nrow(dobs), nrow=mc$sweep_group) # log prior density of sampling transmission pair categories
      mc$pars$LOG_LAMBDA <- matrix(NA_real_, ncol=nrow(dobs), nrow=mc$sweep_group)    #proportions
      
      mc$it.info <- data.table(	IT= seq.int(mc$sweep_group*mc$sweep*(update.round%/%mc$sweep_group)+1,mc$sweep_group*mc$sweep*(1+(update.round%/%mc$sweep_group))),
                                PAR_ID= rep(NA_integer_, mc$sweep_group*mc$sweep),
                                BLOCK= rep(NA_character_, mc$sweep_group*mc$sweep),
                                MHRATIO= rep(NA_real_, mc$sweep_group*mc$sweep),
                                ACCEPT=rep(NA_integer_, mc$sweep_group*mc$sweep),
                                LOG_LKL=rep(NA_real_, mc$sweep_group*mc$sweep),
                                LOG_PRIOR=rep(NA_real_, mc$sweep_group*mc$sweep)
      )            
    }                
  }
}
