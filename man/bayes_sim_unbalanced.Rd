% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_sim_unbalanced.R
\name{bayes_sim_unbalanced}
\alias{bayes_sim_unbalanced}
\title{Unbalanced Bayesian Simulation in Conjugate Linear Model Framework}
\usage{
bayes_sim_unbalanced(
  n1,
  n2,
  repeats = 1,
  u,
  C,
  Xn = NULL,
  Vn = NULL,
  Vbeta_d,
  Vbeta_a_inv,
  sigsq,
  mu_beta_d,
  mu_beta_a,
  alt,
  alpha,
  mc_iter,
  surface_plot = TRUE
)
}
\arguments{
\item{n1}{first sample size (vector or scalar).}

\item{n2}{second sample size (vector or scalar).}

\item{repeats}{an positive integer specifying number of times to repeat
\code{c(n1, n2)}. Applicable for studies that consider multiple measures
within each group. Default setting is \code{repeats = 1} if not applicable.}

\item{u}{a scalar or vector to evaluate \deqn{u'\beta > C,}
where \eqn{\beta} is an unknown parameter that is to be estimated.
Default setting is \code{u = 1}.}

\item{C}{constant value to be compared to when evaluating \eqn{u'\beta > C}}

\item{Xn}{design matrix that characterizes where the data is to be
generated from. This is specifically designed under the normal linear
regression model \deqn{yn = Xn\beta + \epsilon,}
\deqn{\epsilon ~ N(0, \sigma^2 Vn).} When set to \code{NULL},
\code{Xn} is generated in-function using \code{bayesassurance::gen_Xn()}.
Note that setting \code{Xn = NULL} also enables user to pass in a
vector of sample sizes to undergo evaluation.}

\item{Vn}{a correlation matrix for the marginal distribution of the
sample data \code{yn}. Takes on an identity matrix when set to \code{NULL}.}

\item{Vbeta_d}{correlation matrix that helps describe the prior
information on \eqn{\beta} in the design stage}

\item{Vbeta_a_inv}{inverse-correlation matrix that helps describe the prior
information on \eqn{\beta} in the analysis stage}

\item{sigsq}{a known and fixed constant preceding all correlation matrices
\code{Vn}, \code{Vbeta_d} and \code{Vbeta_a_inv}.}

\item{mu_beta_d}{design stage mean}

\item{mu_beta_a}{analysis stage mean}

\item{alt}{specifies alternative test case, where alt = "greater" tests if
\eqn{u'\beta > C},
alt = "less" tests if \eqn{u'\beta < C}, and alt = "two.sided"
performs a two-sided test. By default, alt = "greater".}

\item{alpha}{significance level}

\item{mc_iter}{number of MC samples evaluated under the analysis objective}

\item{surface_plot}{when set to \code{TRUE} and \eqn{n1} and \eqn{n2} are vectors,
a contour map showcasing various assurance values corresponding to
different combinations of \eqn{n1} and \eqn{n2} is produced.}
}
\value{
a list of objects corresponding to the assurance approximations
\itemize{
\item{assurance_table:} table of sample size and corresponding assurance
values
\item{contourplot:} contour map of assurance values
\item{mc_samples:} number of Monte Carlo samples that were generated
and evaluated
}
}
\description{
Approximates the Bayesian assurance of attaining \eqn{u'\beta > C}
for unbalanced study designs through Monte Carlo sampling.
See Argument descriptions for more detail.
}
\examples{
## Example 1
## Sample size vectors are passed in for n1 and n2 to evaluate
## assurance. 
n1 <- seq(20, 75, 5)
n2 <- seq(50, 160, 10)

assur_out <- bayes_sim_unbalanced(n1 = n1, n2 = n2, repeats = 1, u = c(1, -1),
C = 0, Xn = NULL, Vbeta_d = matrix(c(50, 0, 0, 10),nrow = 2, ncol = 2),
Vbeta_a_inv = matrix(rep(0, 4), nrow = 2, ncol = 2),
Vn = NULL, sigsq = 100,  mu_beta_d = c(1.17, 1.25),
mu_beta_a = c(0, 0), alt = "two.sided", alpha = 0.05, mc_iter = 1000,
surface_plot = FALSE)

assur_out$assurance_table


## Example 2
## We can produce a contour plot that evaluates unique combinations of n1
## and n2 simply by setting `surfaceplot = TRUE`.
\donttest{
n1 <- seq(20, 75, 5)
n2 <- seq(50, 160, 10)
assur_out <- bayes_sim_unbalanced(n1 = n1, n2 = n2, repeats = 1, 
u = c(1, -1), C = 0, Xn = NULL, Vbeta_d = matrix(c(50, 0, 0, 10),
nrow = 2, ncol = 2), Vbeta_a_inv = matrix(rep(0, 4), nrow = 2, ncol = 2),
Vn = NULL, sigsq = 100,  mu_beta_d = c(1.17, 1.25),
mu_beta_a = c(0, 0), alt = "two.sided", alpha = 0.05, mc_iter = 1000,
surface_plot = TRUE)

assur_out$assurance_table
assur_out$contourplot
}
}
