\name{mcmc.sarws.biass}
\alias{mcmc.sarws.biass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC sampler under SARWS (bias: sequencing rate).
}
\description{
This function implement MCMC sampler in order to estimate the transmission flows between subpopulations under the assumption that individuals are sampled at random within each subpopulation. Bias are introduced by different sequencing rates among subpopulations.
}
\usage{
mcmc.sarws.biass(dc,iteration,seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dc}{
 data table includes the observed transmissions and sampling information for each subpopulation pair.
}
\item{iteration}{
the number of iterations.
}
\item{seed}{
the seed of random number generator.
}
}

\examples{
mcmc.sarws.biass(dc,1e6,20)
}
