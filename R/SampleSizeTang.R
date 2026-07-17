#' Sample size for a two-group negative binomial rate comparison (Tang, 2015)
#'
#' Computes the total sample size required to detect a rate ratio between
#' two groups of negative binomial event counts, using the Wald test on the
#' log rate ratio with the variance evaluated at the alternative. Follow-up
#' time is either uniform on \code{[0, tau]} with a fraction of dropouts or
#' exponentially shortened by dropout. The parameterization is
#' \code{var(Y) = mu + mu^2 / nu}, matching \code{\link{GenSimData}} and the
#' estimators in this package.
#'
#' @param rate_T,rate_C Positive numerics. True event rates per unit
#'   follow-up time in the treatment and control groups.
#' @param tau Positive numeric. Common planned follow-up period.
#' @param nu Positive numeric dispersion parameter (the
#'   \code{rnbinom(size = nu)} convention). Tang (2015) uses
#'   \code{kappa = 1 / nu}, so a dispersion reported as \code{kappa}
#'   corresponds to \code{nu = 1 / kappa}.
#' @param dropout_T,dropout_C Numerics in \code{[0, 1)}. Dropout proportions
#'   in the treatment and control groups.
#' @param alloc Positive numeric. Allocation ratio \code{n_T : n_C}
#'   (\code{1} for 1:1, \code{2} for 2:1). Defaults to 1.
#' @param alpha Numeric in \code{(0, 1)}. Two-sided significance level.
#'   Defaults to 0.05.
#' @param beta Numeric in \code{(0, 1)}. Type II error. Defaults to 0.20.
#' @param dist Character. Follow-up mechanism, \code{"uniform"} (default) or
#'   \code{"exponential"}.
#'
#' @return A list with components \code{N} (total sample size, rounded up),
#'   \code{n_T}, \code{n_C}, the per-subject expected information
#'   \code{k_T} and \code{k_C}, \code{log_RR}, and the allocation
#'   proportions \code{C_T} and \code{C_C}.
#'
#' @details
#' The per-subject expected information follows Tang (2015). Under uniform
#' follow-up,
#' \code{k_j = (1 - dropout_j) nu rate_j tau / (nu + rate_j tau)
#'   + dropout_j nu (1 - (nu / (rate_j tau)) log(1 + rate_j tau / nu))}.
#' Under exponential dropout with rate
#' \code{lambda_j = -log(1 - dropout_j) / tau},
#' \code{k_j = nu rate_j tau exp(-lambda_j tau) / (nu + rate_j tau)}
#' plus the integral over \code{[0, tau]} of
#' \code{nu rate_j lambda_j t exp(-lambda_j t) / (nu + rate_j t)}.
#' The total sample size is
#' \code{N = (z_{1 - alpha/2} + z_{1 - beta})^2 / (log RR)^2
#'   (1 / (C_T k_T) + 1 / (C_C k_C))} with \code{C_T = alloc / (1 + alloc)}
#' and \code{C_C = 1 / (1 + alloc)}.
#'
#' @examples
#' # Tang (2015) Table 2: kappa = 1, control rate 0.4, rate ratio 0.5,
#' # tau = 2, 20 percent dropout, 90 percent power gives 138 per group.
#' SampleSizeTang(rate_T = 0.2, rate_C = 0.4, tau = 2, nu = 1,
#'                dropout_T = 0.2, dropout_C = 0.2, beta = 0.10)$n_C
#'
#' @importFrom stats qnorm integrate
#' @export
SampleSizeTang <- function(rate_T, rate_C, tau, nu,
                           dropout_T = 0, dropout_C = 0, alloc = 1,
                           alpha = 0.05, beta = 0.20,
                           dist = c("uniform", "exponential")) {

  dist <- match.arg(dist)
  stopifnot(rate_T > 0, rate_C > 0, tau > 0, nu > 0,
            dropout_T >= 0, dropout_T < 1, dropout_C >= 0, dropout_C < 1,
            alloc > 0, alpha > 0, alpha < 1, beta > 0, beta < 1)

  C_T <- alloc / (1 + alloc)
  C_C <- 1 / (1 + alloc)

  # Per-subject expected negative binomial information
  k_per_subject <- function(rate_j, dropout_j) {
    if (dist == "uniform") {
      term1 <- (1 - dropout_j) * nu * rate_j * tau / (nu + rate_j * tau)
      if (dropout_j == 0) return(term1)
      term2 <- dropout_j * nu *
        (1 - (nu / (rate_j * tau)) * log(1 + rate_j * tau / nu))
      term1 + term2
    } else {
      if (dropout_j == 0) {
        return(nu * rate_j * tau / (nu + rate_j * tau))
      }
      lambda_j <- -log(1 - dropout_j) / tau
      term1 <- nu * rate_j * tau * exp(-lambda_j * tau) / (nu + rate_j * tau)
      integrand <- function(t) {
        nu * rate_j * lambda_j * t * exp(-lambda_j * t) / (nu + rate_j * t)
      }
      term2 <- stats::integrate(integrand, lower = 0, upper = tau,
                                rel.tol = 1e-10)$value
      term1 + term2
    }
  }

  k_T <- k_per_subject(rate_T, dropout_T)
  k_C <- k_per_subject(rate_C, dropout_C)

  log_RR <- log(rate_T / rate_C)
  z_a <- stats::qnorm(1 - alpha / 2)
  z_b <- stats::qnorm(1 - beta)
  N_raw <- ((z_a + z_b)^2 / log_RR^2) * (1 / (C_T * k_T) + 1 / (C_C * k_C))
  N   <- ceiling(N_raw)
  n_T <- ceiling(C_T * N)
  n_C <- ceiling(C_C * N)

  list(N = N, n_T = n_T, n_C = n_C, k_T = k_T, k_C = k_C,
       log_RR = log_RR, C_T = C_T, C_C = C_C)
}
