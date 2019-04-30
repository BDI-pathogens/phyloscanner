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
				warning('\nColumn name ',control$flowratios[[ii]][2],' not in MCMC output. Setting to 0.')
				set(z, NULL, control$flowratios[[ii]][2], 0)
			}
			if(!control$flowratios[[ii]][3]%in%colnames(z))
			{
				warning('\nColumn name ',control$flowratios[[ii]][3],' not in MCMC output. Setting to 0.')
				set(z, NULL, control$flowratios[[ii]][3], 0)
			}
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
