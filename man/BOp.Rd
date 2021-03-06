% Generated by roxygen2 (4.0.1): do not edit by hand
\name{BOp}
\alias{BOp}
\title{Compute the posterior of the fitted values from a Bayesian
ordination on the probability scale}
\usage{
BOp(BO, eta, .seed = 1)
}
\arguments{
\item{BO}{Results of \code{\link{getBO}} (not required if
\code{eta} is present)}

\item{eta}{Results of \code{\link{BOeta}} (not required if
\code{BO} is present)}

\item{.seed}{Seed for random number generator associated with the
standard normal residual error term}
}
\value{
An array containing the posterior of the fitted values on
the probability scale with dimensions relating to (1) MCMC
iterations, (2) sites, and (3) species
}
\description{
Compute the posterior of the fitted values from a Bayesian
ordination on the probability scale
}

