#' Compare operating characteristics of rate ratio estimators
#'
#' Runs a Monte Carlo comparison of three rate ratio estimators for a
#' single two-group negative binomial scenario: the negative binomial
#' maximum likelihood estimator (\code{MASS::glm.nb}) used as the
#' reference, the closed-form quasi-Poisson benchmark, and the proposed
#' one-step estimator (\code{\link{EstRROneStep}}). Returns per-replicate
#' estimates and a summary of accuracy, agreement with the reference,
#' confidence interval coverage, and computation time.
#'
#' one-step estimator (\code{\link{EstRROneStep}}). Returns per-replicate
#' estimates and a summary of accuracy, agreement with the reference,
#' confidence interval coverage, and computation time.
#'
#' @param nsim Integer. Number of Monte Carlo replicates. Defaults to 1000.
#' @param n_T,n_C Positive integers. Treatment and control group sizes.
#' @param rate_T,rate_C Positive numerics. Event rates per unit time.
#' @param nu Positive numeric dispersion parameter (the
#'   \code{rnbinom(size = nu)} convention, \code{var(Y) = mu + mu^2 / nu}).
#'   Use \code{Inf} for the Poisson limit.
#' @param tau Positive numeric. Planned follow-up duration. Defaults to 1.
#' @param dropout_T,dropout_C Numerics in \code{[0, 1)}. Group dropout
#'   proportions passed to \code{\link{GenSimData}}. Zero (default) gives
#'   constant follow-up.
#' @param dist Character. Dropout mechanism, \code{"uniform"} (default) or
#'   \code{"exponential"}, passed to \code{\link{GenSimData}}.
#' @param conf_level Numeric in \code{(0, 1)}. Confidence level for the
#'   two-sided Wald interval used in coverage. Defaults to 0.95.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{estimates}{A \code{data.frame} with one row per replicate per
#'       method and columns \code{simID}, \code{method}, \code{log_RR},
#'       \code{SE_log_RR}, \code{covered}.}
#'     \item{summary}{A \code{data.frame} with one row per method and
#'       columns \code{method}, \code{n_valid}, \code{na_rate},
#'       \code{bias}, \code{rmse}, \code{mad_vs_NB}, \code{rmse_vs_NB},
#'       \code{coverage}, and \code{time_per_fit} (elapsed seconds per
#'       replicate). Accuracy columns are relative to the true log rate
#'       ratio, and \code{mad_vs_NB} and \code{rmse_vs_NB} measure
#'       agreement with the negative binomial maximum likelihood estimator.}
#'     \item{truth}{The true log rate ratio \code{log(rate_T / rate_C)}.}
#'     \item{scenario}{A one-row \code{data.frame} of the input parameters.}
#'   }
#'
#' @details
#' The negative binomial reference is fitted by \code{MASS::glm.nb} with a
#' treatment indicator and offset \code{log(followup)}. The quasi-Poisson
#' benchmark uses the Poisson anchor \code{log((S_T / T_T) / (S_C / T_C))}
#' with a negative binomial standard error obtained from the Pearson
#' moment dispersion plug-in evaluated at the anchor. Replicates in which
#' one group has no events are recorded as \code{NA}. Computation time is
#' measured in bulk per method and reported per replicate, which is more
#' stable than timing individual closed-form calls.
#'
#' @examples
#' \donttest{
#' oc <- CompareOC(nsim = 200, n_T = 150, n_C = 150,
#'                 rate_T = 0.5, rate_C = 1.0, nu = 1.5,
#'                 dropout_T = 0.2, dropout_C = 0.2, seed = 1)
#' oc$summary
#' }
#'
#' @importFrom MASS glm.nb
#' @importFrom stats coef vcov qnorm median
#' @export
CompareOC <- function(nsim = 1000L, n_T, n_C, rate_T, rate_C, nu, tau = 1,
                      dropout_T = 0, dropout_C = 0,
                      dist = c("uniform", "exponential"),
                      conf_level = 0.95, seed = NULL) {

  dist <- match.arg(dist)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single numeric value in (0, 1)")
  }

  dat <- GenSimData(nsim = nsim, n_T = n_T, n_C = n_C,
                    rate_T = rate_T, rate_C = rate_C, nu = nu, tau = tau,
                    dropout_T = dropout_T, dropout_C = dropout_C,
                    dist = dist, seed = seed)

  truth  <- log(rate_T / rate_C)
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  cnt <- dat$count
  fu  <- dat$followup
  j   <- as.integer(dat$group != "C")
  idx <- split(seq_len(nrow(dat)), dat$simID)
  nsim <- length(idx)

  # Negative binomial maximum likelihood reference via MASS::glm.nb
  fit_nb <- function(ii) {
    c_s <- cnt[ii]; f_s <- fu[ii]; j_s <- j[ii]
    if (sum(c_s[j_s == 1L]) == 0 || sum(c_s[j_s == 0L]) == 0) {
      return(c(NA_real_, NA_real_))
    }
    d_s <- data.frame(c_s = c_s, j_s = j_s, f_s = f_s)
    fit <- suppressWarnings(tryCatch(
      MASS::glm.nb(c_s ~ j_s + offset(log(f_s)), data = d_s),
      error = function(e) NULL))
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    vc <- tryCatch(stats::vcov(fit)[2L, 2L], error = function(e) NA_real_)
    c(unname(stats::coef(fit)[2L]), sqrt(vc))
  }

  # Closed-form quasi-Poisson benchmark with negative binomial SE
  fit_qp <- function(ii) {
    c_s <- cnt[ii]; f_s <- fu[ii]; j_s <- j[ii]
    S_1 <- sum(c_s[j_s == 1L]); S_0 <- sum(c_s[j_s == 0L])
    T_1 <- sum(f_s[j_s == 1L]); T_0 <- sum(f_s[j_s == 0L])
    if (S_1 == 0 || S_0 == 0) return(c(NA_real_, NA_real_))
    theta_1 <- S_1 / T_1; theta_0 <- S_0 / T_0
    b  <- log(theta_1 / theta_0)
    mu <- f_s * ifelse(j_s == 1L, theta_1, theta_0)
    n  <- length(c_s)
    phi <- if (n > 2L) sum((c_s - mu)^2 / mu) / (n - 2L) else NA_real_
    nu_p <- if (is.finite(phi) && phi > 1 + 1e-8) mean(mu) / (phi - 1) else Inf
    w  <- if (is.finite(nu_p)) nu_p * mu / (nu_p + mu) else mu
    W_T <- sum(w[j_s == 1L]); W_C <- sum(w[j_s == 0L])
    c(b, sqrt(1 / W_T + 1 / W_C))
  }

  # Proposed one-step estimator
  fit_os <- function(ii) {
    r <- EstRROneStep(cnt[ii], fu[ii], dat$group[ii], "C",
                      conf_level = conf_level)
    c(r$log_RR, r$SE_log_RR)
  }

  run_method <- function(fn) {
    t_elapsed <- system.time(res <- lapply(idx, fn))[["elapsed"]]
    mat <- do.call(rbind, res)
    list(log_RR = mat[, 1L], SE = mat[, 2L], time = t_elapsed)
  }

  methods <- c(NB = fit_nb, QP = fit_qp, OneStep = fit_os)
  out <- lapply(methods, run_method)

  # Per-replicate estimates in long form
  est <- do.call(rbind, lapply(names(out), function(m) {
    data.frame(simID     = seq_len(nsim),
               method    = m,
               log_RR    = out[[m]]$log_RR,
               SE_log_RR = out[[m]]$SE,
               stringsAsFactors = FALSE)
  }))
  est$covered <- abs(est$log_RR - truth) <= z_crit * est$SE_log_RR

  # Summary, including agreement with the NB reference aligned by simID
  nb_logRR <- out[["NB"]]$log_RR
  summ <- do.call(rbind, lapply(names(out), function(m) {
    b  <- out[[m]]$log_RR
    dev_nb <- b - nb_logRR
    ok <- is.finite(b)
    data.frame(
      method       = m,
      n_valid      = sum(ok),
      na_rate      = mean(!ok),
      bias         = mean(b[ok] - truth),
      rmse         = sqrt(mean((b[ok] - truth)^2)),
      mad_vs_NB    = stats::median(abs(dev_nb), na.rm = TRUE),
      rmse_vs_NB   = sqrt(mean(dev_nb^2, na.rm = TRUE)),
      coverage     = mean(abs(b[ok] - truth) <= z_crit * out[[m]]$SE[ok]),
      time_per_fit = out[[m]]$time / nsim,
      stringsAsFactors = FALSE)
  }))

  scenario <- data.frame(nsim = nsim, n_T = n_T, n_C = n_C,
                         rate_T = rate_T, rate_C = rate_C, nu = nu, tau = tau,
                         dropout_T = dropout_T, dropout_C = dropout_C,
                         dist = dist, stringsAsFactors = FALSE)

  list(estimates = est, summary = summ, truth = truth, scenario = scenario)
}
