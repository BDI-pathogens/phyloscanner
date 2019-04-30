#' @title Aggregate MCMC output to target parameters
#' @export
#' @import data.table
#' @param mcmc.file Full file name to MCMC output from function \code{source.attribution.mcmc}
#' @param mc MCMC object from function \code{source.attribution.mcmc}; this can be specified alternatively to `mcmc.file`.
#' @param daggregateTo Data.table that maps the categories of transmission pairs used in the MCMC to lower-dimensional categories that are of primary interest.
#' @param control List of input arguments that control the behaviour of the source attribution algorithm:
#' \itemize{
#'  \item{"burnin.p"}{Proportion of MCMC iterations that are removed as burn-in period.}
#'  \item{"thin"}{Thin MCMC output to every nth MCMC iteration.}
#'  \item{"regex_pars"}{Regular expression to select the parameters that are to be aggregated.}
#'  \item{"outfile"}{Full file name of the output csv file.}
#' }
#' @return If `outfile` is not specified, aggregated Monte Carlo samples are returned, as a data.table. If `outfile` is specified and ends in `csv`, the output is written to a csv file. If `outfile` is specified and ends in `rda`, the output is written to a rda file.
source.attribution.mcmc.aggregateToTarget	<- function(mcmc.file=NULL, mc=NULL, daggregateTo=NULL, control=list(burnin.p=NA_real_, thin=NA_integer_, regex_pars='*', outfile=gsub('\\.rda','_aggregated.csv',mcmc.file))){
	#	basic checks
	if(!'data.table'%in%class(daggregateTo))
		stop('daggregateTo is not a data.table')
	if(!all(c('TRM_CAT_PAIR_ID','TR_TARGETCAT','REC_TARGETCAT')%in%colnames(daggregateTo)))
		stop('daggregateTo does not contain one of the required columns TRM_CAT_PAIR_ID, TR_TARGETCAT, REC_TARGETCAT')
	
	if(!is.null(mc))
	{
		cat('\nUsing MCMC output specified as input...') 
	}
	if(is.null(mc))
	{
		#	load MCMC output
		cat('\nLoading MCMC output from file...')
		load(mcmc.file)		
	}
	
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
