% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keyquantities_zeros.R
\name{source.attribution.mcmc.getKeyQuantities}
\alias{source.attribution.mcmc.getKeyQuantities}
\title{Estimate Flows, Sources, WAIFM, Flow ratios}
\usage{
source.attribution.mcmc.getKeyQuantities(
  infile = NULL,
  mc = NULL,
  pars = NULL,
  dobs = NULL,
  control = NULL
)
}
\arguments{
\item{infile}{Full file name to MCMC output or aggregated MCMC output from function \code{source.attribution.mcmc}}

\item{control}{List of input arguments that control the behaviour of the derivation of key quantities:
\itemize{
 \item{"quantiles"}{Named list of quantiles. Default: c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975)}
 \item{"flowratios"}{Cector of length 3. First element: name of flow ratio. Second element: name of transmission pair category for enumerator of flow ratio. Third element: name of transmission pair category for denominator of flow ratio.}
 \item{"outfile"}{Full file name for output csv file.}
}}
}
\value{
If `outfile` is not specified, the output is returned as a data.table. If `outfile` is specified, the output is written to a csv file.
}
\description{
Estimate Flows, Sources, WAIFM, Flow ratios
}
\examples{
require(data.table)
require(phyloflows)
# load input data set
data(twoGroupFlows1, package="phyloflows")
dobs <- twoGroupFlows1$dobs
dprior <- twoGroupFlows1$dprior
# load MCMC output obtained with
# control <- list(seed=42, mcmc.n=5e4, verbose=0)
# mc <- phyloflows:::source.attribution.mcmc(dobs, dprior, control)
data(twoGroupFlows1_mcmc, package="phyloflows")

#
# calculate sources, onward flows & flow ratios from MCMC output
#
control <- list( quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                 flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
                 )
ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(mc=mc, dobs=dobs, control=control)

#
# calculate sources, onward flows & flow ratios from aggregated MCMC output
#
daggregateTo <- subset(dobs, select=c(TRM_CAT_PAIR_ID, TR_TRM_CATEGORY, REC_TRM_CATEGORY))
daggregateTo[, TR_TARGETCAT:= TR_TRM_CATEGORY]
daggregateTo[, REC_TARGETCAT:= 'Any']
set(daggregateTo, NULL, c('TR_TRM_CATEGORY','REC_TRM_CATEGORY'), NULL)
control <- list( burnin.p=0.05, thin=NA_integer_, regex_pars='PI')
pars <- phyloflows:::source.attribution.mcmc.aggregateToTarget(mc=mc, daggregateTo=daggregateTo, control=control)
\dontrun{
control <- list( quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                 flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
                 )
ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(pars=pars, control=control)
}

#
# calculate sources, onward flows & flow ratios from file
#
\dontrun{
infile<- 'blah.rda'     # either MCMC output or aggregated MCMC output
control <- list( quantiles= c('CL'=0.025,'IL'=0.25,'M'=0.5,'IU'=0.75,'CU'=0.975),
                 flowratios= list( c('1/2', '1 2', '2 1'), c('2/1', '2 1', '1 2'))
                 )
ans <- phyloflows:::source.attribution.mcmc.getKeyQuantities(infile=infile, dobs=dobs, control=control)
}

}
\seealso{
[source.attribution.mcmc()], [source.attribution.mcmc.aggregateToTarget()]
}
\author{
Xiaoyue Xi, Oliver Ratmann
}
