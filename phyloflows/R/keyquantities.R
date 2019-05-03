#' @title Estimate Flows, Sources, WAIFM, Flow ratios
#' @export
#' @import data.table
#' @author Xiaoyue Xi, Oliver Ratmann
#' @param infile Full file name to MCMC output or aggregated MCMC output from function \code{source.attribution.mcmc}
#' @param control List of input arguments that control the behaviour of the derivation of key quantities:
#' \itemize{
#'  \item{"quantiles"}{Named list of quantiles. Default: c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975)}
#'  \item{"flowratios"}{Cector of length 3. First element: name of flow ratio. Second element: name of transmission pair category for enumerator of flow ratio. Third element: name of transmission pair category for denominator of flow ratio.}
#'  \item{"outfile"}{Full file name for output csv file.}
#' }
#' @return If `outfile` is not specified, the output is returned as a data.table. If `outfile` is specified, the output is written to a csv file.
#' @seealso \link{\code{source.attribution.mcmc}}, \link{\code{source.attribution.mcmc.aggregateToTarget}}
#' @examples 
#' require(data.table)
#' require(phyloflows)
#' # load input data set
#' data(twoGroupFlows1, package="phyloflows")
#' dobs <- twoGroupFlows1$dobs
#' dprior <- twoGroupFlows1$dprior
#' # load MCMC output obtained with
#' # control <- list(seed=42, mcmc.n=5e4, verbose=0)
#' # mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)
#' data(twoGroupFlows1_mcmc, package="phyloflows")
#' 
#' #
#' # calculate sources, onward flows & flow ratios from MCMC output 
#' #
#' control <- list( quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
#'                  flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
#'                  )
#' ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(mc=mc, dobs=dobs, control=control)
#' 
#' #
#' # calculate sources, onward flows & flow ratios from aggregated MCMC output 
#' #
#' daggregateTo <- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
#' daggregateTo[, TR_TARGETCAT:= TR_TRM_CATEGORY]
#' daggregateTo[, REC_TARGETCAT:= 'Any']
#' set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
#' control <- list( burnin.p=0.05, thin=NA_integer_, regex_pars='PI')
#' pars <- phyloflows:::source.attribution.mcmc.aggregateToTarget(mc=mc, daggregateTo=daggregateTo, control=control)
#' \dontrun{
#' control <- list( quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
#'                  flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
#'                  )
#' ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(pars=pars, control=control)
#' }
#' 
#' #
#' # calculate sources, onward flows & flow ratios from file 
#' #
#' \dontrun{
#' infile<- 'blah.rda'     # either MCMC output or aggregated MCMC output
#' control <- list( quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
#'                  flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
#'                  )
#' ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(infile=infile, dobs=dobs, control=control)
#' }
#' 
source.attribution.mcmc.getKeyQuantities<- function(infile=NULL, mc=NULL, pars=NULL, dobs=NULL, control=NULL)
{
	if(is.null(mc) & is.null(pars) & !is.null(infile) && grepl('rda$',infile))
	{
		tmp			<- load(infile)
		if(tmp=='pars')
		{
			cat('\nReading aggregated MCMC output...')			
		}
		if(tmp=='mc')
		{
			cat('\nReading MCMC output...')			
		}				
	}
	if(is.null(mc) & is.null(pars) & !is.null(infile) && grepl('csv$',infile))
	{
		cat('\nReading aggregated MCMC output...')
		pars		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))			
	}
	if(!is.null(mc))
	{
		if(is.null(dobs))
			stop('Found unaggregated MCMC output, need dobs to continue')
		# make data.table in long format
		pars	<- mc[['pars']][['PI']]
		colnames(pars)	<- paste0('PI-',seq_len(ncol(mc[['pars']][['PI']])))
		pars	<- as.data.table(pars)
		pars[, SAMPLE:= seq_len(nrow(pars))]
		pars	<- melt(pars, id.vars='SAMPLE')
		pars[, VARIABLE:= pars[, gsub('([A-Z]+)-([0-9]+)','\\1',variable)]]
		pars[, TRM_CAT_PAIR_ID:= pars[, as.integer(gsub('([A-Z]+)-([0-9]+)','\\2',variable))]]
		if(!all(sort(unique(pars$TRM_CAT_PAIR_ID))==sort(unique(dobs$TRM_CAT_PAIR_ID))))
			stop('The transmission count categories in the MCMC output do not match the transmission count categories in the aggregateTo data table.')		
		pars	<- merge(pars, dobs, by='TRM_CAT_PAIR_ID')
		setnames(pars, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), c('TR_TARGETCAT','REC_TARGETCAT'))
		pars	<- pars[, list(VALUE=sum(value)), by=c('VARIABLE','TR_TARGETCAT','REC_TARGETCAT','SAMPLE')]		
	}
	pars		<- subset(pars, VARIABLE=='PI')
	if(any(is.na(pars$VALUE)))
	{
		cat('\nRemoving NA output for samples n=', nrow(subset(pars, is.na(VALUE))))
		pars		<- subset(pars, !is.na(VALUE))
	}
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
			if(!control$flowratios[[ii]][2]%in%colnames(z))
			{
				cat('\nColumn names in MCMC output are', colnames(z))
				warning('\nColumn name ',control$flowratios[[ii]][2],' not in MCMC output. Setting to 0.')
				set(z, NULL, control$flowratios[[ii]][2], 0)
			}
			if(!control$flowratios[[ii]][3]%in%colnames(z))
			{
				cat('\nColumn names in MCMC output are', colnames(z))
				warning('\nColumn name ',control$flowratios[[ii]][3],' not in MCMC output. Setting to 0.')
				set(z, NULL, control$flowratios[[ii]][3], 0)
			}
			set(z, NULL, control$flowratios[[ii]][1], z[[ control$flowratios[[ii]][2] ]] /  z[[ control$flowratios[[ii]][3] ]])			
		}
		for(ii in seq_along(control$flowratios))
		{
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
	
	if(!'outfile'%in%names(control))
		return(ans)	
	cat('\nWriting output to',control$outfile)
	write.csv(ans, row.names=FALSE,file=control$outfile)
}
