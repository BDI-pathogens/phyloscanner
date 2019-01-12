\name{mcmc.sarws.diagnostics}
\alias{mcmc.sarws.diagnostics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic of MCMC under the assumption of SARWS.
}
\description{
This function makes convergence diagnostics for MCMC sampler.
}
\usage{
mcmc.sarws.diagnostics(mc, iteration, burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mc}{
   a list containing information on MCMC sampler.
}
\item{iteration}{
  the number of iterations for MCMC sampler.
}
\item{burnin}{
  the number of steps used for burn-in.
}
}

\examples{
mcmc.sarws.diagnostics(mc, 1e6, 2e3)
}
