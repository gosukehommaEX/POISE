#' Reconstruct negative binomial rate ratio inference from group summaries
#'
#' Reconstructs the two-group negative binomial rate ratio estimate, its
#' standard error, the Wald confidence interval, and the Wald test from
#' group-level summary statistics, without access to subject-level data.
#' When the follow-up time is constant within each group the returned
#' quantities equal the negative binomial maximum likelihood estimates
#' that would be obtained from the full subject-level data. The intended
#' use is meta-analysis and reanalysis of published trials that report
#' group event rates, follow-up, and a dispersion parameter but do not
#' release individual patient data.
#'
#' @param rate_T Positive numeric. Estimated event rate per unit time in
#'   the treatment group, that is total events divided by total follow-up.
#' @param rate_C Positive numeric. Estimated event rate per unit time in
#'   the control group.
#' @param followup_T Positive numeric. Common follow-up time per subject
#'   in the treatment group, that is total follow-up divided by
#'   \code{n_T}.
#' @param followup_C Positive numeric. Common follow-up time per subject
#'   in the control group.
#' @param n_T Positive integer. Number of subjects in the treatment group.
#' @param n_C Positive integer. Number of subjects in the control group.
#' @param nu Positive numeric dispersion parameter, parameterized so that
#'   \code{var(Y) = mu + mu^2 / nu} (the \code{rnbinom(size = nu)}
#'   convention). Use \code{Inf} for the Poisson limit. When a trial
#'   reports the dispersion as \code{k = 1 / nu}, pass \code{nu = 1 / k}.
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
#'     \item{nu}{Dispersion parameter used.}
#'     \item{n_T, n_C}{Treatment and control group sizes.}
#'     \item{alternative}{The alternative hypothesis used.}
#'     \item{conf_level}{The confidence level used.}
#'     \item{method}{A label describing the estimator.}
#'   }
#'
#' @details
#' Let \code{mu_j = followup_j rate_j} denote the mean event count per
#' subject in group j, so that the total event count is
#' \code{S_j = n_j mu_j}. The log rate ratio is
#' \code{log(rate_T / rate_C)}, and the standard error is
#' \code{sqrt(1 / W_T + 1 / W_C)} with
#' \code{W_j = n_j nu mu_j / (nu + mu_j)} the negative binomial
#' information contributed by group j. In the Poisson limit
#' (\code{nu = Inf}) each \code{W_j} reduces to the total event count
#' \code{S_j}. When the follow-up time varies within a group the same
#' expressions provide the quasi-Poisson point estimate and a
#' dispersion-adjusted standard error, which approximate rather than
#' reproduce the subject-level maximum likelihood fit.
#'
#' @examples
#' # Reconstruct an analysis from reported group summaries
#' EstRRSummary(rate_T = 0.9, rate_C = 1.5,
#'              followup_T = 1, followup_C = 1,
#'              n_T = 200, n_C = 200, nu = 1.5)
#'
#' @importFrom stats pnorm qnorm
#' @export
EstRRSummary <- function(rate_T, rate_C, followup_T, followup_C,
                         n_T, n_C, nu, conf_level = 0.95,
                         alternative = c("two.sided", "less", "greater")) {

  # Input validation
  scalars <- list(rate_T = rate_T, rate_C = rate_C,
                  followup_T = followup_T, followup_C = followup_C,
                  n_T = n_T, n_C = n_C, nu = nu)
  for (nm in names(scalars)) {
    v <- scalars[[nm]]
    if (!is.numeric(v) || length(v) != 1L) {
      stop(sprintf("'%s' must be a single numeric value", nm))
    }
  }
  if (rate_T <= 0 || rate_C <= 0) {
    stop("'rate_T' and 'rate_C' must be strictly positive")
  }
  if (followup_T <= 0 || followup_C <= 0) {
    stop("'followup_T' and 'followup_C' must be strictly positive")
  }
  if (n_T < 1 || n_C < 1) {
    stop("'n_T' and 'n_C' must be at least 1")
  }
  if (nu <= 0) {
    stop("'nu' must be strictly positive (use Inf for the Poisson limit)")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single numeric value in (0, 1)")
  }
  alternative <- match.arg(alternative)

  # Mean event count per subject and negative binomial group information
  mu_T <- followup_T * rate_T
  mu_C <- followup_C * rate_C
  W_T <- if (is.finite(nu)) n_T * nu * mu_T / (nu + mu_T) else n_T * mu_T
  W_C <- if (is.finite(nu)) n_C * nu * mu_C / (nu + mu_C) else n_C * mu_C

  log_RR    <- log(rate_T / rate_C)
  SE_log_RR <- sqrt(1 / W_T + 1 / W_C)

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
         nu = nu, n_T = n_T, n_C = n_C,
         alternative = alternative, conf_level = conf_level,
         method = "Negative binomial rate ratio from group summaries"),
    class = "rrfit")
}
