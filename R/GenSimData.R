#' Generate two-group negative binomial recurrent event data
#'
#' Generates Monte Carlo replicates of a two-group parallel trial with
#' subject-specific follow-up times and negative binomial event counts.
#' The output is arranged for direct use by the rate ratio estimators in
#' this package. Follow-up may be constant within group (no dropout, the
#' default) or heterogeneous through a uniform or exponential dropout
#' mechanism, which is the regime in which the closed-form one-step
#' correction departs from the Poisson anchor.
#'
#' @param nsim Integer. Number of Monte Carlo replicates. Defaults to 1000.
#' @param n_T Positive integer. Number of subjects in the treatment group
#'   (group label \code{"T"}).
#' @param n_C Positive integer. Number of subjects in the control group
#'   (group label \code{"C"}).
#' @param rate_T Positive numeric. Event rate per unit time in the
#'   treatment group.
#' @param rate_C Positive numeric. Event rate per unit time in the control
#'   group.
#' @param nu Positive numeric dispersion parameter, parameterized so that
#'   \code{var(Y) = mu + mu^2 / nu} (the \code{rnbinom(size = nu)}
#'   convention). Use \code{Inf} for the Poisson limit, in which counts are
#'   drawn from a Poisson distribution.
#' @param tau Positive numeric. Planned follow-up duration. Defaults to 1.
#' @param dropout_T Numeric in \code{[0, 1)}. Dropout proportion in the
#'   treatment group. \code{0} (default) gives constant follow-up
#'   \code{tau} for every subject in the group.
#' @param dropout_C Numeric in \code{[0, 1)}. Dropout proportion in the
#'   control group.
#' @param dist Character. Dropout mechanism for subjects who do not
#'   complete \code{tau}. \code{"uniform"} draws the follow-up of a
#'   \code{dropout_j} fraction of subjects uniformly on \code{(0, tau]} and
#'   sets the remainder to \code{tau}. \code{"exponential"} draws each
#'   follow-up as \code{min(rexp(rate), tau)} with the rate chosen so that
#'   the probability of dropping out before \code{tau} equals
#'   \code{dropout_j}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A \code{data.frame} with one row per subject per replicate and
#'   columns \code{simID} (replicate index), \code{group} (factor with
#'   levels \code{"T"} and \code{"C"}), \code{followup}, and \code{count}.
#'   The \code{count}, \code{followup}, and \code{group} columns are the
#'   inputs expected by \code{\link{EstRROneStep}} with \code{control = "C"}.
#'
#' @examples
#' d <- GenSimData(nsim = 50, n_T = 100, n_C = 100,
#'                 rate_T = 0.5, rate_C = 1.0, nu = 2, tau = 1,
#'                 dropout_T = 0.1, dropout_C = 0.1, dist = "uniform",
#'                 seed = 1)
#' head(d)
#'
#' @importFrom stats runif rexp rnbinom rpois
#' @export
GenSimData <- function(nsim = 1000L, n_T, n_C, rate_T, rate_C, nu, tau = 1,
                       dropout_T = 0, dropout_C = 0,
                       dist = c("uniform", "exponential"), seed = NULL) {

  # Input validation
  if (!is.null(seed)) set.seed(seed)
  dist <- match.arg(dist)
  if (nsim < 1 || n_T < 1 || n_C < 1) {
    stop("'nsim', 'n_T', and 'n_C' must be at least 1")
  }
  if (rate_T <= 0 || rate_C <= 0) {
    stop("'rate_T' and 'rate_C' must be strictly positive")
  }
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0) {
    stop("'nu' must be a single positive numeric value (Inf for Poisson)")
  }
  if (tau <= 0) {
    stop("'tau' must be strictly positive")
  }
  if (dropout_T < 0 || dropout_T >= 1 || dropout_C < 0 || dropout_C >= 1) {
    stop("'dropout_T' and 'dropout_C' must be in [0, 1)")
  }

  nsim <- as.integer(nsim)
  n_T  <- as.integer(n_T)
  n_C  <- as.integer(n_C)
  N_per_sim <- n_T + n_C
  total_n   <- N_per_sim * nsim

  # Subject-level group labels and rates, blocked by replicate
  group_one <- rep(c("T", "C"), c(n_T, n_C))
  group_all <- rep(group_one, nsim)
  rate_all  <- ifelse(group_all == "T", rate_T, rate_C)

  # Follow-up times
  if (dist == "uniform") {
    # A fixed fraction of subjects per group drop out with follow-up
    # uniform on (0, tau]; the remainder complete tau.
    n_drop_T <- ceiling(dropout_T * n_T)
    n_drop_C <- ceiling(dropout_C * n_C)
    drop_one <- c(rep(c(1L, 0L), c(n_drop_T, n_T - n_drop_T)),
                  rep(c(1L, 0L), c(n_drop_C, n_C - n_drop_C)))
    drop_all <- rep(drop_one, nsim)
    followup <- ifelse(drop_all == 1L,
                       stats::runif(total_n, min = 0, max = tau),
                       tau)
  } else {
    # Exponential dropout with administrative censoring at tau, rate set
    # so that P(drop out before tau) equals the group dropout proportion.
    dropout_all <- ifelse(group_all == "T", dropout_T, dropout_C)
    lambda_all  <- ifelse(dropout_all > 0, -log(1 - dropout_all) / tau, 0)
    raw_fu <- ifelse(lambda_all > 0,
                     stats::rexp(total_n, rate = pmax(lambda_all, 1e-12)),
                     Inf)
    followup <- pmin(raw_fu, tau)
  }

  # Negative binomial (or Poisson limit) counts
  mu <- followup * rate_all
  count <- if (is.finite(nu)) {
    stats::rnbinom(total_n, mu = mu, size = nu)
  } else {
    stats::rpois(total_n, lambda = mu)
  }

  data.frame(
    simID    = rep(seq_len(nsim), each = N_per_sim),
    group    = factor(group_all, levels = c("T", "C")),
    followup = followup,
    count    = as.integer(count)
  )
}
