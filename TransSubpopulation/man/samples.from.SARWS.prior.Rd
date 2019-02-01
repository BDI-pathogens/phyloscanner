\name{samples.from.SARWS.prior}
\alias{samples.from.SARWS.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drawing samples from the prior distribution under sampling at random within strata (SARWS) assumption.
}
\description{
This function provides the inputs for mcmc sampler. Samples are from the prior distribution of the sampling rate and evaluate corresponding density. The prior distribution is informative and given by a posterior distribution inferred from Bayesian hierarchical model using empirical sampling rate under the assumption of missing at random within strata (SARWS).
}
\usage{
samples.from.SARWS.prior(dc,iteration,biastype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{iteration}{
the number of iterations for mcmc sampler, which decides the number of prior samples needed.
}
\item{dc}{
 a data table includes transmission counts between subpopulations and covariates of transmitters' and recipients' subpopulations.
}
\item{biastype}{
the type of bias considered. biastype='sequence' indicates only the bias introduced by sequencing rate is considered, biastype='participation' indicates only the bias introduced by sequencing rate is considered
}
}

\examples{
# need to create a simulated dc
samples.from.SARWS.prior(dc,iteration=1e6,biastype='sequence')
}
