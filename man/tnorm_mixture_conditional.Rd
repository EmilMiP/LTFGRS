% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helper_Functions.R
\name{tnorm_mixture_conditional}
\alias{tnorm_mixture_conditional}
\title{Title: Calculates mean and variance of mixture of two truncated normal distributions}
\usage{
tnorm_mixture_conditional(mu, var, lower, upper, K_i, K_pop)
}
\arguments{
\item{mu}{mean value of normal distribution}

\item{var}{variance of normal distribution}

\item{lower}{lower threshold (can be -Inf)}

\item{upper}{upper threshold (can be Inf)}

\item{Kp}{(stratified) cumulative incidence proportion for the individual}
}
\value{
mean and variance of mixture distribution between two truncated normal distributions
}
\description{
Title: Calculates mean and variance of mixture of two truncated normal distributions
}
\examples{
tnorm_mixture_conditional(mu = 0, var = 1, lower = -Inf, upper = Inf, Kp = 0)
tnorm_mixture_conditional(mu = 0, var = 1, lower = -Inf, upper = 2, Kp = .01)
}
