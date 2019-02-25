\name{mcmc.core.inference}
\alias{mcmc.core.inference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC sampler to correct for sampling bias.
}
\description{
This function implements MCMC sampler in order to estimate the transmission flows between subpopulations when different sampling rates among subpopulations introduce bias.
}
\usage{
mcmc.core.inference(pair,S.DTL.prior,TR_OBS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pair}{
a matrix whose first and second column are the indices of subpopulations for transmitters and recipients.
}
\item{S.DTL.prior}{
a list which contains samples of sampling rates from prior, corresponding densities and empirical sampling rates.
}
\item{TR_OBS}{
the number of observed transmission counts between subpopulations.
}
}
\examples{
# run example in samples.from.GLM.prior or samples.from.SARWS.prior first
# load the data saved from running these functions
mcmc.core.inference(pair=unique.pairs.index,S.DTL.prior=S.DTL.prior,TR_OBS=TR_OBS)
}
