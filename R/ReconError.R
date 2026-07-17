#' Reconstruction error of the summary-based rate ratio estimator
#'
#' Quantifies how closely the summary-based reconstruction
#' (\code{\link{EstRRSummary}}) reproduces the subject-level negative
#' binomial maximum likelihood fit (\code{MASS::glm.nb}) when only
#' group-level summaries are available. For each Monte Carlo replicate the
#' function fits the negative binomial model to the full data, forms the
#' group summaries a published report would provide (group event rate,
#' mean follow-up, group size, and the estimated dispersion), reconstructs
#' the rate ratio and its standard error from those summaries, and records
#' the reconstruction error. When the follow-up time is constant within
#' each group the reconstruction is exact up to numerical precision. When
#' the follow-up time varies within a group the residual error measures the
#' information lost in aggregation, reported both in absolute terms and
#' relative to the maximum likelihood standard error.
#'
#' @param nsim Integer. Number of Monte Carlo replicates. Defaults to 1000.
#' @param n_T,n_C Positive integers. Treatment and control group sizes.
#' @param rate_T,rate_C Positive numerics. Event rates per unit time.
#' @param nu Positive numeric dispersion parameter (the
#'   \code{rnbinom(size = nu)} convention, \code{var(Y) = mu + mu^2 / nu}).
#' @param tau Positive numeric. Planned follow-up duration. Defaults to 1.
#' @param dropout_T,dropout_C Numerics in \code{[0, 1)}. Group dropout
#'   proportions passed to \code{\link{GenSimData}}. Zero (default) gives
#'   constant follow-up, under which the reconstruction is exact.
#' @param dist Character. Dropout mechanism, \code{"uniform"} (default) or
#'   \code{"exponential"}, passed to \code{\link{GenSimData}}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{estimates}{A \code{data.frame} with one row per replicate and
#'       columns \code{simID}, \code{log_RR_NB}, \code{SE_NB},
#'       \code{nu_NB}, \code{log_RR_rec}, \code{SE_rec},
#'       \code{err_log_RR} (reconstruction minus reference log rate ratio),
#'       \code{err_SE}, and \code{rel_err} (absolute log rate ratio error
#'       divided by the reference standard error).}
#'     \item{summary}{A one-row \code{data.frame} with \code{n_valid},
#'       \code{na_rate}, and the median and maximum of the absolute
#'       reconstruction errors and of the relative error.}
#'     \item{truth}{The true log rate ratio \code{log(rate_T / rate_C)}.}
#'     \item{scenario}{A one-row \code{data.frame} of the input parameters.}
#'   }
#'
#' @details
#' The reference is \code{MASS::glm.nb} fitted with a treatment indicator
#' and offset \code{log(followup)}. The reported dispersion is taken to be
#' the fitted \code{theta}. The reconstruction uses the group event rates
#' \code{S_j / T_j}, the mean follow-up \code{T_j / n_j}, the group sizes,
#' and the fitted dispersion, which are the quantities a published trial
#' provides. Replicates in which one group has no events, or in which
#' \code{glm.nb} does not return a fit, are recorded as \code{NA}.
#'
#' @examples
#' \donttest{
#' re <- ReconError(nsim = 200, n_T = 150, n_C = 150,
#'                  rate_T = 0.5, rate_C = 1.0, nu = 1.5,
#'                  dropout_T = 0.2, dropout_C = 0.2, seed = 1)
#' re$summary
#' }
#'
#' @importFrom MASS glm.nb
#' @importFrom stats coef vcov median
#' @export
ReconError <- function(nsim = 1000L, n_T, n_C, rate_T, rate_C, nu, tau = 1,
                       dropout_T = 0, dropout_C = 0,
                       dist = c("uniform", "exponential"), seed = NULL) {

  dist <- match.arg(dist)
  dat <- GenSimData(nsim = nsim, n_T = n_T, n_C = n_C,
                    rate_T = rate_T, rate_C = rate_C, nu = nu, tau = tau,
                    dropout_T = dropout_T, dropout_C = dropout_C,
                    dist = dist, seed = seed)

  truth <- log(rate_T / rate_C)
  cnt <- dat$count
  fu  <- dat$followup
  j   <- as.integer(dat$group != "C")
  idx <- split(seq_len(nrow(dat)), dat$simID)
  nsim <- length(idx)

  one_rep <- function(ii) {
    c_s <- cnt[ii]; f_s <- fu[ii]; j_s <- j[ii]
    S_1 <- sum(c_s[j_s == 1L]); S_0 <- sum(c_s[j_s == 0L])
    if (S_1 == 0 || S_0 == 0) {
      return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
    }
    # Negative binomial maximum likelihood reference on the full data
    d_s <- data.frame(c_s = c_s, j_s = j_s, f_s = f_s)
    fit <- suppressWarnings(tryCatch(
      MASS::glm.nb(c_s ~ j_s + offset(log(f_s)), data = d_s),
      error = function(e) NULL))
    if (is.null(fit)) {
      return(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
    }
    log_RR_NB <- unname(stats::coef(fit)[2L])
    SE_NB <- tryCatch(sqrt(stats::vcov(fit)[2L, 2L]),
                      error = function(e) NA_real_)
    nu_NB <- fit$theta

    # Group summaries a published report would provide, then reconstruct
    n_1 <- sum(j_s == 1L); n_0 <- sum(j_s == 0L)
    T_1 <- sum(f_s[j_s == 1L]); T_0 <- sum(f_s[j_s == 0L])
    rec <- EstRRSummary(rate_T = S_1 / T_1, rate_C = S_0 / T_0,
                        followup_T = T_1 / n_1, followup_C = T_0 / n_0,
                        n_T = n_1, n_C = n_0, nu = nu_NB)
    c(log_RR_NB, SE_NB, nu_NB, rec$log_RR, rec$SE_log_RR)
  }

  mat <- do.call(rbind, lapply(idx, one_rep))
  est <- data.frame(
    simID      = seq_len(nsim),
    log_RR_NB  = mat[, 1L],
    SE_NB      = mat[, 2L],
    nu_NB      = mat[, 3L],
    log_RR_rec = mat[, 4L],
    SE_rec     = mat[, 5L],
    stringsAsFactors = FALSE)
  est$err_log_RR <- est$log_RR_rec - est$log_RR_NB
  est$err_SE     <- est$SE_rec - est$SE_NB
  est$rel_err    <- abs(est$err_log_RR) / est$SE_NB

  ok <- is.finite(est$err_log_RR)
  summ <- data.frame(
    n_valid            = sum(ok),
    na_rate            = mean(!ok),
    med_abs_err_log_RR = stats::median(abs(est$err_log_RR[ok])),
    max_abs_err_log_RR = max(abs(est$err_log_RR[ok])),
    med_abs_err_SE     = stats::median(abs(est$err_SE[ok])),
    med_rel_err        = stats::median(est$rel_err[ok]),
    max_rel_err        = max(est$rel_err[ok]),
    stringsAsFactors = FALSE)

  scenario <- data.frame(nsim = nsim, n_T = n_T, n_C = n_C,
                         rate_T = rate_T, rate_C = rate_C, nu = nu, tau = tau,
                         dropout_T = dropout_T, dropout_C = dropout_C,
                         dist = dist, stringsAsFactors = FALSE)

  list(estimates = est, summary = summ, truth = truth, scenario = scenario)
}
