\name{samples.from.GLM.prior}
\alias{samples.from.GLM.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drawing samples from the prior distribution given by generalised linear model.
}
\description{
This function provides the inputs for mcmc sampler. Samples are from the prior distribution of the sampling rate and evaluate corresponding density. The prior distribution is informative and given by a posterior distribution inferred from Bayesian hierarchical model using empirical sampling rate and group covariates.
}
\usage{
samples.from.GLM.prior(dp,dg,dc,iteration,method,biastype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dp}{
 a data table includes data which infer the distribution of participation rates through hierarchical model (i.e. the number of eligible individuals and individual participating in each subpopulation, the covariates of individuals in a subpopulation). It is not required if merely bias from sequencing rate is accounted for.
}
\item{dg}{
 a data table includes data which infer the distribution of sequencing rates through hierarchical model (i.e. the number of infected (and ART naive) individuals and individual sequenced in each subpopulation, the covariates of individuals in a subpopulation).
}
\item{dc}{
 a data table includes transmission counts between subpopulations and covariates of transmitters' and recipients' subpopulations.
}
\item{iteration}{
the number of iterations for mcmc sampler, which decides the number of prior samples needed.
}
\item{method}{
the methods for evaluating density of each samples. method='empirical' indicates density of samples come from an empirical distribution, method='betaapprox' indicates a beta distribution is used to approximate the distribution of samples.
}
\item{biastype}{
the type of bias considered. biastype='sequence' indicates only the bias introduced by sequencing rate is considered, biastype='participation' indicates only the bias introduced by sequencing rate is considered
}
}

\examples{
# need to create a simulated dp, dg, dc
samples.from.GLM.prior(dp,dg,dc,iteration=1e6,method='empirical',biastype='sequence')
}
