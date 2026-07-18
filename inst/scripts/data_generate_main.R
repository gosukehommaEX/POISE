# Generate and cache the Monte Carlo results for the POISE manuscript.
#
# This script runs the simulation once and saves the numeric results to
# inst/scripts/data as .rds files. The companion script
# table_and_figure_manuscript.R reads those files and produces the figures
# and tables, so plotting can be revised without rerunning the simulation.
# All paths are resolved relative to this script's folder.
#
# Approximate runtime scales with NSIM; NSIM = 2000 takes several minutes.
# Lower NSIM for a quick test, raise to 10000 for the final manuscript.

library(POISE)

# ---- Resolve paths relative to this script (inst/scripts) ----
get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(sub("^--file=", "", file_arg[1])))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(p))
  }
  getwd()
}
script_dir <- get_script_dir()
data_dir   <- file.path(script_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Global settings ----
NSIM  <- 2000L      # raise to 10000 for the final manuscript
TAU   <- 1
ALPHA <- 0.05
BETA  <- 0.10       # 90 percent power
DROP_DESIGN <- 0.2  # dropout assumed at the design stage
BASE_SEED   <- 20260718

# Phase 2 uses a larger effect (RR = 0.5); Phase 3 a smaller one (RR = 0.7).
phase_RR <- c("Phase 2" = 0.5, "Phase 3" = 0.7)
rate_C   <- 1.0
nu_levels <- c(0.5, 1, 2)

# ============================================================
# Part 1: operating characteristics (accuracy, coverage, time)
# ============================================================
oc_grid <- expand.grid(phase = names(phase_RR),
                       nu    = nu_levels,
                       alloc = c(1, 2),
                       stringsAsFactors = FALSE)

oc_summary   <- list()
oc_estimates <- list()
m <- 0L
for (i in seq_len(nrow(oc_grid))) {
  ph <- oc_grid$phase[i]; nu <- oc_grid$nu[i]; alloc <- oc_grid$alloc[i]
  RR_alt <- phase_RR[[ph]]
  # Size the trial at the alternative, then reuse the same n for the
  # matched null (rate ratio 1) so the null run yields the type I error.
  ss <- SampleSizeTang(rate_T = RR_alt * rate_C, rate_C = rate_C, tau = TAU,
                       nu = nu, dropout_T = DROP_DESIGN, dropout_C = DROP_DESIGN,
                       alloc = alloc, alpha = ALPHA, beta = BETA,
                       dist = "uniform")
  for (effect in c("Alternative", "Null")) {
    m <- m + 1L
    rate_T <- if (effect == "Alternative") RR_alt * rate_C else rate_C
    oc <- CompareOC(nsim = NSIM, n_T = ss$n_T, n_C = ss$n_C,
                    rate_T = rate_T, rate_C = rate_C, nu = nu, tau = TAU,
                    dropout_T = DROP_DESIGN, dropout_C = DROP_DESIGN,
                    dist = "uniform", conf_level = 0.95,
                    seed = BASE_SEED + m)
    lab <- data.frame(phase = ph, nu = nu, alloc = alloc, effect = effect,
                      RR = rate_T / rate_C, n_T = ss$n_T, n_C = ss$n_C,
                      stringsAsFactors = FALSE)
    oc_summary[[m]] <- cbind(lab[rep(1L, nrow(oc$summary)), ], oc$summary)
    est <- oc$estimates
    est$phase <- ph; est$nu <- nu; est$alloc <- alloc; est$effect <- effect
    oc_estimates[[m]] <- est
    cat(sprintf("OC %d: %s %s nu=%.1f alloc=%d n=(%d,%d)\n",
                m, ph, effect, nu, alloc, ss$n_T, ss$n_C))
  }
}
oc_summary   <- do.call(rbind, oc_summary)
oc_estimates <- do.call(rbind, oc_estimates)
rownames(oc_summary) <- NULL
saveRDS(list(summary = oc_summary, estimates = oc_estimates, grid = oc_grid),
        file.path(data_dir, "sim_oc.rds"))

# ============================================================
# Part 2: reconstruction error vs follow-up heterogeneity
# ============================================================
drop_grid <- c(0, 0.1, 0.2, 0.4)
re_design <- expand.grid(phase = names(phase_RR), nu = nu_levels,
                         stringsAsFactors = FALSE)

recon_summary   <- list()
recon_estimates <- list()
k <- 0L
for (r in seq_len(nrow(re_design))) {
  ph <- re_design$phase[r]; nu <- re_design$nu[r]
  RR <- phase_RR[[ph]]; rate_T <- RR * rate_C
  # Fix the design size at the design dropout, then vary the actual dropout
  # so the figure isolates the effect of follow-up heterogeneity.
  ss <- SampleSizeTang(rate_T = rate_T, rate_C = rate_C, tau = TAU, nu = nu,
                       dropout_T = DROP_DESIGN, dropout_C = DROP_DESIGN,
                       alloc = 1, alpha = ALPHA, beta = BETA, dist = "uniform")
  for (d in drop_grid) {
    k <- k + 1L
    re <- ReconError(nsim = NSIM, n_T = ss$n_T, n_C = ss$n_C,
                     rate_T = rate_T, rate_C = rate_C, nu = nu, tau = TAU,
                     dropout_T = d, dropout_C = d, dist = "uniform",
                     seed = BASE_SEED + 1000L + k)
    lab <- data.frame(phase = ph, nu = nu, RR = RR, dropout = d,
                      n_T = ss$n_T, n_C = ss$n_C, stringsAsFactors = FALSE)
    recon_summary[[k]] <- cbind(lab, re$summary)
    est <- re$estimates
    est$phase <- ph; est$nu <- nu; est$dropout <- d
    recon_estimates[[k]] <- est
    cat(sprintf("Recon %d/%d: %s nu=%.3f drop=%.1f\n",
                k, nrow(re_design) * length(drop_grid), ph, nu, d))
  }
}
recon_summary   <- do.call(rbind, recon_summary)
recon_estimates <- do.call(rbind, recon_estimates)
rownames(recon_summary) <- NULL
saveRDS(list(summary = recon_summary, estimates = recon_estimates,
             drop_grid = drop_grid),
        file.path(data_dir, "sim_recon.rds"))

cat("\nDone. Saved sim_oc.rds and sim_recon.rds to:\n  ", data_dir, "\n")
