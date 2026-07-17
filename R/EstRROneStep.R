#' One-step closed-form rate ratio estimator for two-group negative binomial data
#'
#' Estimates the rate ratio for a two-group recurrent event trial under a
#' negative binomial model without iterative maximum likelihood fitting.
#' The estimator anchors at the Poisson group rates \code{S_j / T_j},
#' plugs in a Pearson moment estimator of the dispersion parameter, and
#' applies a single group-wise Newton (Fisher scoring) step to the
#' negative binomial score for the group rates. The construction is
#' non-iterative and reproduces the negative binomial maximum likelihood
#' estimator exactly when the follow-up time is constant within each
#' group. Returns an object of class \code{"rrfit"} carrying the point
#' estimate, the standard error of the log rate ratio, the two-sided Wald
#' confidence interval, the Wald test statistic, and the p-value.
#'
#' @param count A non-negative integer vector of event counts for all
#'   subjects (pooled over both groups).
#' @param followup A positive numeric vector of subject-level follow-up
#'   times, aligned with \code{count}.
#' @param group A vector of group labels aligned with \code{count}. Any
#'   type that supports equality comparison is accepted.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group (j = 0). Subjects with
#'   \code{group != control} are treated as the experimental treatment
#'   group (j = 1).
#' @param nu Optional positive scalar dispersion parameter, parameterized
#'   so that \code{var(Y) = mu + mu^2 / nu} (the \code{rnbinom(size = nu)}
#'   convention). When \code{NULL} (default), \code{nu} is estimated in
#'   closed form by the Pearson moment plug-in at the Poisson anchor. Pass
#'   a value when a dispersion estimate is available externally, for
#'   example when reconstructing an analysis from reported summaries.
#' @param conf_level Numeric in \code{(0, 1)}. Confidence level for the
#'   two-sided Wald confidence interval. Defaults to 0.95.
#' @param alternative A character string specifying the alternative
#'   hypothesis for the Wald test of \code{log_RR = 0}, one of
#'   \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#'   The choice sets the p-value only. The confidence interval is always
#'   the two-sided interval at \code{conf_level}, following the usual
#'   clinical trial reporting convention.
#'
#' @return An object of class \code{"rrfit"}, a list with elements:
#'   \describe{
#'     \item{RR}{Point estimate of the rate ratio.}
#'     \item{log_RR}{Point estimate of the log rate ratio.}
#'     \item{SE_log_RR}{Standard error of \code{log_RR}.}
#'     \item{CI_RR}{Numeric vector of length 2 with the lower and upper
#'       two-sided Wald confidence limits for the rate ratio.}
#'     \item{Z}{Wald test statistic for testing \code{log_RR = 0}.}
#'     \item{p_value}{P-value for the Wald test (two-sided or one-sided
#'       according to \code{alternative}).}
#'     \item{nu}{Dispersion parameter used (supplied or plug-in estimate).}
#'     \item{n_T, n_C}{Treatment and control group sizes.}
#'     \item{alternative}{The alternative hypothesis used.}
#'     \item{conf_level}{The confidence level used.}
#'     \item{method}{A label describing the estimator.}
#'   }
#'   All numeric elements are \code{NA_real_} if the estimate cannot be
#'   computed (e.g. all events in one group, zero events overall, or a
#'   numerical failure in the closed-form step).
#'
#' @details
#' Let \code{S_j} and \code{T_j} denote the total event count and total
#' follow-up time in group j, so that \code{theta_j = S_j / T_j} is the
#' Poisson anchor rate. Conditional on the dispersion \code{nu}, the
#' negative binomial score for the two group rates separates by group,
#' and \code{theta_j} solves
#' \code{sum_i (Y_ij - t_ij theta_j) / (nu + t_ij theta_j) = 0}.
#' A single Newton step from the Poisson anchor \code{S_j / T_j} yields the
#' proposed group rate. When the follow-up time is constant within a
#' group the anchor already solves this equation, so the step vanishes and
#' the estimator equals the negative binomial maximum likelihood estimator.
#' The standard error uses the negative binomial information at the fitted
#' rates, \code{sqrt(1 / W_T + 1 / W_C)} with
#' \code{W_j = sum_i nu mu_ij / (nu + mu_ij)}. In the Poisson limit
#' (\code{nu = Inf}) the step vanishes and \code{W_j} reduces to the total
#' fitted mean in group j.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' fu <- runif(n, 0.5, 1)
#' grp <- rep(c("T", "C"), each = n / 2)
#' theta <- ifelse(grp == "T", 0.5, 1)
#' y <- rnbinom(n, mu = fu * theta, size = 2)
#' EstRROneStep(y, fu, grp, "C")
#'
#' @importFrom stats pnorm qnorm
#' @export
EstRROneStep <- function(count, followup, group, control,
                         nu = NULL, conf_level = 0.95,
                         alternative = c("two.sided", "less", "greater")) {

  # Input validation
  if (!is.null(nu)) {
    if (!is.numeric(nu) || length(nu) != 1L || nu <= 0) {
      stop("'nu' must be NULL or a single positive numeric value")
    }
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single numeric value in (0, 1)")
  }
  alternative <- match.arg(alternative)
  method_label <- "One-step negative binomial rate ratio"
  n <- length(count)
  if (length(followup) != n || length(group) != n) {
    stop("'count', 'followup', and 'group' must have the same length")
  }
  na_result <- structure(
    list(RR = NA_real_, log_RR = NA_real_, SE_log_RR = NA_real_,
         CI_RR = c(NA_real_, NA_real_), Z = NA_real_, p_value = NA_real_,
         nu = NA_real_, n_T = NA_integer_, n_C = NA_integer_,
         alternative = alternative, conf_level = conf_level,
         method = method_label),
    class = "rrfit")
  if (n == 0L) return(na_result)
  if (any(followup <= 0)) {
    stop("'followup' must be strictly positive for all subjects")
  }
  if (any(count < 0) || any(count != as.integer(count))) {
    stop("'count' must contain non-negative integers")
  }

  # Treatment indicator: 1 = treatment, 0 = control
  j <- as.integer(group != control)

  # Group-level totals and Poisson anchor rates
  S_0 <- sum(count[j == 0L])
  S_1 <- sum(count[j == 1L])
  T_0 <- sum(followup[j == 0L])
  T_1 <- sum(followup[j == 1L])
  n_0 <- sum(j == 0L)
  n_1 <- sum(j == 1L)

  if (n_0 == 0L || n_1 == 0L) return(na_result)
  if (S_0 == 0 || S_1 == 0) return(na_result)

  theta_0 <- S_0 / T_0
  theta_1 <- S_1 / T_1

  # Dispersion: supplied value or Pearson moment plug-in at the anchor.
  # E(phi_hat) ~ 1 + mu_bar / nu, so nu_plug = mu_bar / (phi_hat - 1).
  # phi_hat <= 1 + eps signals no overdispersion, giving the Poisson
  # limit nu = Inf under which the estimator reduces to the anchor.
  mu_anchor <- followup * ifelse(j == 1L, theta_1, theta_0)
  if (is.null(nu)) {
    df_resid <- n - 2L
    phi_hat <- if (df_resid > 0L) {
      sum((count - mu_anchor)^2 / mu_anchor) / df_resid
    } else NA_real_
    mu_bar  <- sum(mu_anchor) / n
    eps_phi <- 1e-8
    nu_use <- if (is.finite(phi_hat) && phi_hat > 1 + eps_phi) {
      mu_bar / (phi_hat - 1)
    } else Inf
  } else {
    nu_use <- nu
  }

  # Single group-wise Newton step on the negative binomial rate score.
  # g(theta) = sum_i (Y_i - t_i theta) / (nu + t_i theta) = 0.
  newton_rate <- function(theta0, y_g, t_g) {
    if (!is.finite(nu_use)) return(theta0)
    denom <- nu_use + t_g * theta0
    g     <- sum((y_g - t_g * theta0) / denom)
    gp    <- -sum(t_g * (nu_use + y_g) / denom^2)
    if (!is.finite(gp) || gp == 0) return(NA_real_)
    theta0 - g / gp
  }

  theta_1_hat <- newton_rate(theta_1, count[j == 1L], followup[j == 1L])
  theta_0_hat <- newton_rate(theta_0, count[j == 0L], followup[j == 0L])
  if (!is.finite(theta_1_hat) || !is.finite(theta_0_hat) ||
      theta_1_hat <= 0 || theta_0_hat <= 0) {
    return(na_result)
  }

  log_RR <- log(theta_1_hat / theta_0_hat)

  # Standard error from the negative binomial information at the fit.
  mu_hat <- followup * ifelse(j == 1L, theta_1_hat, theta_0_hat)
  w_hat  <- if (is.finite(nu_use)) nu_use * mu_hat / (nu_use + mu_hat) else mu_hat
  W_T <- sum(w_hat[j == 1L])
  W_C <- sum(w_hat[j == 0L])
  if (W_T <= 0 || W_C <= 0) return(na_result)
  SE_log_RR <- sqrt(1 / W_T + 1 / W_C)

  if (!is.finite(log_RR) || !is.finite(SE_log_RR) || SE_log_RR <= 0) {
    return(na_result)
  }

  # Two-sided Wald confidence interval and p-value keyed to alternative.
  Z      <- log_RR / SE_log_RR
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  CI_log <- log_RR + c(-1, 1) * z_crit * SE_log_RR
  p_value <- switch(alternative,
                    two.sided = 2 * stats::pnorm(-abs(Z)),
                    less      = stats::pnorm(Z),
                    greater   = stats::pnorm(-Z))
  structure(
    list(RR = exp(log_RR), log_RR = log_RR, SE_log_RR = SE_log_RR,
         CI_RR = exp(CI_log), Z = Z, p_value = p_value,
         nu = nu_use, n_T = n_1, n_C = n_0,
         alternative = alternative, conf_level = conf_level,
         method = method_label),
    class = "rrfit")
}
