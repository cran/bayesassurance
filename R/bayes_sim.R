#' Bayesian Simulation in Conjugate Linear Model Framework
#'
#' Approximates the Bayesian assurance of attaining either \eqn{u'\beta > C},
#' \eqn{u'\beta < C}, or \eqn{u'\beta \neq C},
#' for equal-sized samples through Monte Carlo sampling.
#' The function also carries the capability to process longitudinal data.
#' See Argument descriptions for more detail.
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline aes theme element_text 
#' @importFrom ggplot2 xlab ylab ggtitle scale_x_continuous scale_y_continuous
#' @importFrom rlang .data
#' @importFrom pbapply pbsapply
#' @importFrom stats qnorm rnorm
#' @import mathjaxr
#' @param n sample size (either scalar or vector). When `longitudinal = TRUE`, 
#' `n` denotes the number of observations per subject.
#' @param p column dimension of design matrix `Xn`. If `Xn = NULL`,
#' `p` must be specified to denote the column dimension of the default 
#' design matrix generated by the function.
#' @param u a scalar or vector included in the expression to be evaluated, 
#' e.g. \deqn{u'\beta > C,} where \eqn{\beta} is an unknown parameter 
#' that is to be estimated.
#' @param C constant to be compared
#' @param Xn design matrix that characterizes where the data is to be 
#' generated from. This is specifically given by the normal linear 
#' regression model \deqn{yn = Xn\beta + \epsilon,} 
#' \deqn{\epsilon ~ N(0, \sigma^2 Vn).} When set to `NULL`, `Xn` is generated 
#' in-function using either `bayesassurance::gen_Xn()` or 
#' `bayesassurance::gen_Xn_longitudinal()`. Note that setting `Xn = NULL` 
#' also enables user to pass in a vector of sample sizes to undergo evaluation 
#' as the function will automatically adjust `Xn` accordingly based on the 
#' sample size.
#' @param Vn a correlation matrix for the marginal distribution of the sample 
#' data \eqn{yn}. Takes on an identity matrix when set to `NULL`.
#' @param mu_beta_d design stage mean
#' @param mu_beta_a analysis stage mean
#' @param Vbeta_d correlation matrix that helps describe the prior information 
#' on \eqn{\beta} in the design stage
#' @param Vbeta_a_inv inverse-correlation matrix that helps describe the prior 
#' information on \eqn{\beta} in the analysis stage
#' @param sigsq a known and fixed constant preceding all correlation matrices 
#' `Vn`, `Vbeta_d`, and `Vbeta_a_inv`.
#' @param alt specifies alternative test case, where alt = "greater" tests if 
#' \eqn{u'\beta > C}, alt = "less" tests if \eqn{u'\beta < C}, and 
#' alt = "two.sided" performs a two-sided test. By default, alt = "greater".
#' @param alpha significance level
#' @param mc_iter number of MC samples evaluated under the analysis objective
#' @param longitudinal when set to `TRUE`, constructs design matrix using 
#' inputs that correspond to a balanced longitudinal study design.
#' @param ids vector of unique subject ids, usually of length 2
#' for study design purposes.
#' @param from start time of repeated measures for
#' each subject
#' @param to end time of repeated measures for each
#' subject
#' @param poly_degree only needed if `longitudinal = TRUE`, specifies 
#' highest degree taken in the longitudinal model. 
#' @return a list of objects corresponding to the assurance approximations
#' \itemize{
#'      \item{assurance_table:} table of sample size and corresponding assurance
#'      values
#'      \item{assur_plot:} plot of assurance values
#'      \item{mc_samples:} number of Monte Carlo samples that were generated
#'      and evaluated
#' }
#' @examples 
#' 
#' ## Example 1
#' ## A single Bayesian assurance value obtained given a scalar sample size
#' ## n and p=1. Note that setting p=1 implies that
#' ## beta is a scalar parameter.
#'
#' bayesassurance::bayes_sim(n=100, p = 1, u = 1, C = 0.15, Xn = NULL, 
#' Vbeta_d = 1e-8, Vbeta_a_inv = 0, Vn = NULL, sigsq = 0.265, mu_beta_d = 0.3, 
#' mu_beta_a = 0, alt = "two.sided", alpha = 0.05, mc_iter = 5000)
#'
#'
#' ## Example 2
#' ## Illustrates a scenario in which weak analysis priors and strong 
#' ## design priors are assigned to enable overlap between the frequentist 
#' ## power and Bayesian assurance.
#'
#' \donttest{
#' library(ggplot2)
#' n <- seq(100, 250, 5)
#'
#'  ## Frequentist Power
#'  power <- bayesassurance::pwr_freq(n, sigsq = 0.265, theta_0 = 0.15,
#'  theta_1 = 0.25, alt = "greater", alpha = 0.05)
#'
#'  ## Bayesian simulation values with specified values from the n vector
#'  assurance <- bayesassurance::bayes_sim(n, p = 1, u = 1, C = 0.15, Xn = NULL,
#'  Vbeta_d = 1e-8, Vbeta_a_inv = 0, Vn = NULL, sigsq = 0.265, mu_beta_d = 0.25,
#'  mu_beta_a = 0, alt = "greater", alpha = 0.05, mc_iter = 1000)
#'
#' ## Visual representation of plots overlayed on top of one another
#' df1 <- as.data.frame(cbind(n, power = power$pwr_table$Power))
#' df2 <- as.data.frame(cbind(n, assurance = 
#' assurance$assurance_table$Assurance))
#'
#' plot_curves <- ggplot2::ggplot(df1, alpha = 0.5, ggplot2::aes(x = n, y = power,
#' color="Frequentist")) + ggplot2::geom_line(lwd=1.2)
#' plot_curves <- plot_curves + ggplot2::geom_point(data = df2, alpha = 0.5,
#' aes(x = n, y = assurance, color="Bayesian"),lwd=1.2) +
#' ggplot2::ggtitle("Bayesian Simulation vs. Frequentist Power Computation")
#' plot_curves
#'}
#'
#' ## Example 3
#' ## Longitudinal example where n now denotes the number of repeated measures 
#' ## per subject and design matrix is determined accordingly.
#'
#' ## subject ids
#' n <- seq(10, 100, 5)
#' ids <- c(1,2)
#' sigsq <- 100
#' Vbeta_a_inv <- matrix(rep(0, 16), nrow = 4, ncol = 4)
#' Vbeta_d <- (1 / sigsq) * 
#' matrix(c(4, 0, 3, 0, 0, 6, 0, 0, 3, 0, 4, 0, 0, 0, 0, 6), 
#' nrow = 4, ncol = 4)
#'
#' assur_out <- bayes_sim(n = n, p = NULL, u = c(1, -1, 1, -1), C = 0, 
#'                       Xn = NULL, Vbeta_d = Vbeta_d, 
#'                       Vbeta_a_inv = Vbeta_a_inv,
#'                       Vn = NULL, sigsq = 100,
#'                       mu_beta_d = as.matrix(c(5, 6.5, 62, 84)),
#'                       mu_beta_a = as.matrix(rep(0, 4)), mc_iter = 1000,
#'                       alt = "two.sided", alpha = 0.05, 
#'                       longitudinal = TRUE, ids = ids,
#'                       from = 10, to = 120)
#' assur_out$assurance_plot
#' 
#' @seealso \code{\link{pwr_freq}} for frequentist power function,
#' \code{\link{assurance_nd_na}} for a closed form assurance function, and
#' \code{\link{bayes_sim_unknownvar}} for a Bayesian assurance function
#' assumes unvariance.
#' @export

bayes_sim <- function(n, p = NULL, u, C, Xn = NULL, Vn = NULL, Vbeta_d, 
                      Vbeta_a_inv, sigsq, mu_beta_d, mu_beta_a, 
                      alt = "two.sided", alpha, mc_iter, 
                      longitudinal = FALSE, ids = NULL, from = NULL,
                      to = NULL, poly_degree = NULL){

  is.scalar <- function(x) is.atomic(x) && length(x) == 1L
  
  # test for parameters passed in
  if(is.null(Xn) == FALSE){
    if(dim(Xn)[2] != dim(as.matrix(u))[1]){
    stop("Column dimension of Xn must be equal to row dimension of u.")
    }
  }
  
  
  if(is.null(Xn) == FALSE){
    if(dim(Xn)[2] != dim(as.matrix(mu_beta_d))[1]){
    stop("Column dimension of Xn must be equal to row dimension of mu_beta_d.")
    }
  }
  
  if(is.null(Xn) == FALSE){
    if(dim(Xn)[2] != dim(as.matrix(mu_beta_a))[1]){
    stop("Column dimension of Xn must be equal to row dimension of mu_beta_a.")
    }
  }
  
  if(is.scalar(C) == FALSE){
    stop("C must be a scalar value.")
  }
  
  if(is.scalar(sigsq) == FALSE){
    stop("sigsq must be scalar.")
  }
  
  if(alpha < 0 | alpha > 1){
    stop("Not a valid significance level, alpha must be between 0 and 1.")
  }
  
  if((alt %in% c("greater", "less", "two.sided")) == FALSE | is.null(alt)){
    stop("Please specify one of the three options for alternative test case: 
         greater, less, two.sided.")
  }
  
  # this embedded function is used to see if the values of n meet the
  # analysis stage objective
  MC_sample <- function(n = n){

    # initial checks
    if(is.null(Xn) == FALSE & longitudinal == FALSE){
      p <- ncol(Xn)
      Xn_t <- t(Xn)
    }

    if(is.null(Xn) & longitudinal == FALSE){
      if(is.null(p)){
        stop("Need to specify column dimension of design matrix if design
           matrix wasn't specified in function call.")
      }else{
        Xn <- bayesassurance::gen_Xn(n = rep(n, p))
        Xn_t <- t(Xn)
      }
    }

    if(is.null(Xn) & longitudinal == TRUE){
      if(is.null(ids) | is.null(from) | is.null(to)){
        stop("Error: At least one longitudinal parameter is not specified.")
      }
      Xn <- bayesassurance::gen_Xn_longitudinal(ids, from, to, n)
      Xn_t <- t(Xn)
    }

    if(is.null(Vn) == TRUE){
      Vn <- diag(x = 1, nrow = nrow(Xn), ncol = nrow(Xn))
    }

    Vn_inv <- chol2inv(chol(Vn))

    count <- 0 # counter for iterations satisfying the analysis objective

    # Analysis Stage Begins
    L_tran <- chol(Vbeta_a_inv + t(Xn) %*% Vn_inv %*% Xn)
    v <- backsolve(L_tran, Vn)
    M <- v %*% t(v)

    # Components that make up the marginal distribution from which
    # y is to be generated from
    V_star <- sigsq * (Xn %*% Vbeta_d %*% Xn_t + Vn)
    L <- t(chol(V_star))

    # posterior_vals <- c()
    for(i in 1:mc_iter){

      # Design Stage Begins
      z <- matrix(stats::rnorm(nrow(Xn), 0, 1), nrow(Xn), 1)
      yi <- Xn %*% mu_beta_d + L %*% z
      # Design Stage Ends

      # Components that make up the posterior distribution
      m <- Vbeta_a_inv %*% mu_beta_a + Xn_t %*% Vn_inv %*% yi
      Mm <- M %*% m

      # Checks the analysis stage objective based on the specified alternative
      if(alt == "greater"){
        Zi <- ifelse((C - t(u) %*% Mm) / (sqrt(sigsq) * sqrt(t(u) %*% M %*% u))
                     < stats::qnorm(alpha), 1, 0)
        count <- ifelse(Zi == 1, count <- count + 1, count <- count)
      }else if(alt == "less"){
        Zi <- ifelse((C - t(u) %*% Mm) / (sqrt(sigsq) * sqrt(t(u) %*% M %*% u))
                     > stats::qnorm(1-alpha), 1, 0)
        count <- ifelse(Zi == 1, count <- count + 1, count <- count)
      }else if(alt == "two.sided"){
        Zi <- ifelse((C - t(u) %*% Mm) / (sqrt(sigsq) * sqrt(t(u) %*% M %*% u))
                     > stats::qnorm(1-alpha/2) | (C - t(u) %*% Mm) / 
                       (sqrt(sigsq) * sqrt(t(u) %*% M %*% u))
                     < stats::qnorm(alpha/2), 1, 0)
        count <- ifelse(Zi == 1, count <- count + 1, count <- count)
      }

      # Analysis Stage Ends
    }

    assurance <- count / mc_iter
    return(assurance)
  }

  # evaluates the analysis objective for all n by applying
  # MC_sample()
  assurance <- pbapply::pbsapply(n, function(i) MC_sample(n=i))

  # Assurance table
    assur_tab <- as.data.frame(cbind(n, assurance))
    colnames(assur_tab) <- c("Observations per Group (n)", "Assurance")
    assur_tab <- structure(assur_tab, class = "data.frame")


  # ggplot
  if(length(n) > 1){
    assur_plot <- ggplot2::ggplot(assur_tab, alpha = 0.5, 
                                  aes(x = .data$`Observations per Group (n)`,
                                  y = .data$Assurance)) +
      ggplot2::geom_point(aes(x = .data$`Observations per Group (n)`, 
                                  y = .data$Assurance), lwd = 1.2) +
      ggplot2::ggtitle("Assurance Values") +
      ggplot2::xlab("Sample Size n") + ggplot2::ylab("Assurance")
    assur_plot <- structure(assur_plot, class = "ggplot")
  }

  if(length(n) > 1){
    return(list(assurance_table = assur_tab, assurance_plot = assur_plot, 
                mc_samples = mc_iter))
  }else{
    return(list(assur_val = paste0("Assurance: ", round(assurance, 3)), 
                mc_samples = mc_iter))
  }

}
