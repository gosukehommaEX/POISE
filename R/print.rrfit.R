#' Print a rate ratio fit
#'
#' Formats and prints an object of class \code{"rrfit"} returned by the
#' rate ratio estimators in this package. Displays the rate ratio, the log
#' rate ratio with its standard error, the two-sided confidence interval,
#' the Wald test result, the dispersion parameter, and the group sizes,
#' with the value column aligned across rows.
#'
#' @param x An object of class \code{"rrfit"}.
#' @param digits Integer. Number of significant digits for the printed
#'   estimates. Defaults to 3.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return The input \code{x}, returned invisibly.
#'
#' @examples
#' fit <- EstRRSummary(rate_T = 0.9, rate_C = 1.5,
#'                     followup_T = 1, followup_C = 1,
#'                     n_T = 200, n_C = 200, nu = 1.5)
#' print(fit)
#'
#' @export
print.rrfit <- function(x, digits = 3, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  if (!is.finite(x$RR)) {
    cat("  Estimate not available (NA).\n\n")
    return(invisible(x))
  }

  fnum <- function(v) trimws(formatC(v, format = "g", digits = digits))
  alt_lab <- switch(x$alternative,
                    two.sided = "two-sided",
                    less      = "one-sided, less",
                    greater   = "one-sided, greater")
  ci_pct <- format(100 * x$conf_level, trim = TRUE)
  nu_lab <- if (is.finite(x$nu)) fnum(x$nu) else "Inf (Poisson limit)"

  # Assemble label-value rows, then left-justify labels to a common width
  # so that the value column is aligned across all rows.
  rows <- list(
    c("Rate ratio", fnum(x$RR)),
    c("log rate ratio",
      paste0(fnum(x$log_RR), "  (SE ", fnum(x$SE_log_RR), ")")),
    c(paste0(ci_pct, "% CI (rate ratio)"),
      paste0(fnum(x$CI_RR[1]), " to ", fnum(x$CI_RR[2]))),
    c(paste0("Wald test (", alt_lab, ")"),
      paste0("Z = ", fnum(x$Z), ", p = ", fnum(x$p_value))),
    c("Dispersion (nu)", nu_lab)
  )
  if (!is.na(x$n_T) && !is.na(x$n_C)) {
    rows[[length(rows) + 1L]] <-
      c("Group sizes", paste0("n_T = ", x$n_T, ", n_C = ", x$n_C))
  }

  labels <- vapply(rows, `[`, character(1), 1L)
  width  <- max(nchar(labels))
  for (r in rows) {
    cat("  ", formatC(r[1L], width = width, flag = "-"), " : ", r[2L],
        "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}
