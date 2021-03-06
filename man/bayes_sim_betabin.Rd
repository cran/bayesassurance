% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_sim_betabinomial.R
\name{bayes_sim_betabin}
\alias{bayes_sim_betabin}
\title{Bayesian Assurance Computation in the Beta-Binomial Setting}
\usage{
bayes_sim_betabin(
  n1,
  n2,
  p1,
  p2,
  alpha_1,
  alpha_2,
  beta_1,
  beta_2,
  sig_level,
  alt,
  mc_iter
)
}
\arguments{
\item{n1}{sample size of first group}

\item{n2}{sample size of second group}

\item{p1}{proportion of successes in first group. Takes on a NULL (default)
assignment if unknown.}

\item{p2}{proportion of successes in second group. Takes on a NULL
(default) assignment if unknown.}

\item{alpha_1, beta_1}{shape parameters for the distribution of \code{p1}
if \code{p1} is unknown: \eqn{p1 ~ Beta(\alpha_1, \beta_1)}}

\item{alpha_2, beta_2}{shape parameters for the distribution of p2 if
p2 is unknown: \eqn{p2 ~ Beta(\alpha_2, \beta_2)}}

\item{sig_level}{significance level}

\item{alt}{a character string specifying the alternative hypothesis,
must select one of following choices: \code{"two.sided"} (default),
\code{"greater"} or \code{"less"}.}

\item{mc_iter}{number of MC samples evaluated under the analysis objective}
}
\value{
approximate Bayesian assurance of independent two-sample proportion
test
}
\description{
Returns the Bayesian assurance corresponding to a hypothesis test for
difference in two independent proportions.
}
\examples{

#########################################################
# alpha1 = 0.5, beta1 = 0.5, alpha2 = 0.5, beta2 = 0.5 ##
#########################################################
n <- seq(200, 1000, 10)
assur_vals <- bayesassurance::bayes_sim_betabin(n1 = n, n2 = n, 
p1 = 0.25, p2 = 0.2, alpha_1 = 0.5, beta_1 = 0.5, alpha_2 = 0.5, 
beta_2 = 0.5, sig_level = 0.05, alt = "greater", mc_iter = 1000)

assur_vals$assurance_table
assur_vals$assurance_plot

}
