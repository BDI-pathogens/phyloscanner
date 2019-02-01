\name{mcmc.sequence.approximation}
\alias{mcmc.sequence.approximation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC sampler to correct for bias introduced by varying sequencing rates.
}
\description{
This function implement MCMC sampler in order to estimate the transmission flows between subpopulations when different sequencing rates among subpopulations introduce bias.
}
\usage{
mcmc.sequence.approximation(S.DTL.prior,TR_OBS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{S.DTL.prior}{
a list which contains samples of sampling rates from prior, corresponding densities and empirical sampling rates.
}
\item{TR_OBS}{
the number of observed transmission counts between subpopulations.
}
}
\examples{
# TODO: simulated data
mcmc.sequence.approximation(S.DTL.prior,TR_OBS)
}
