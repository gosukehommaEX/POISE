# Build the manuscript figures and tables for POISE.
#
# Reads the cached simulation results from inst/scripts/data and writes
# figures (PDF and EPS) and LaTeX tables to inst/scripts/results. Run
# data_generate_main.R first to produce the .rds inputs. Plotting can be
# revised and rerun here without repeating the simulation. All paths are
# resolved relative to this script's folder.

library(ggplot2)

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
script_dir  <- get_script_dir()
data_dir    <- file.path(script_dir, "data")
results_dir <- file.path(script_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

oc    <- readRDS(file.path(data_dir, "sim_oc.rds"))
recon <- readRDS(file.path(data_dir, "sim_recon.rds"))

# ---- Shared display settings ----
METHOD_LEVELS <- c("NB", "QP", "OneStep")
METHOD_LABELS <- c("NB MLE", "QP", "POISE")
METHOD_COLORS <- c("NB MLE" = "#000000", "QP" = "#D55E00", "POISE" = "#0072B2")

as_method  <- function(m) factor(METHOD_LABELS[match(m, METHOD_LEVELS)],
                                 levels = METHOD_LABELS)
as_nu_lab  <- function(nu) factor(paste0("nu == ", nu),
                                  levels = paste0("nu == ", sort(unique(nu))))
as_alloc   <- function(a) factor(ifelse(a == 1, "1:1", "2:1"),
                                 levels = c("1:1", "2:1"))
as_nu_f    <- function(nu) factor(nu, levels = sort(unique(nu)))

base_theme <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = grid::unit(2, "cm"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

save_fig <- function(p, name, width = 7.2, height = 5.6) {
  ggsave(file.path(results_dir, paste0(name, ".pdf")), p,
         width = width, height = height)
  ggsave(file.path(results_dir, paste0(name, ".eps")), p,
         width = width, height = height, device = "eps")
}

nu_labeller <- labeller(nu_lab = label_parsed)

# ============================================================
# Figure 1: reconstruction error vs follow-up heterogeneity
# ============================================================
rs <- recon$summary
rs$nu_lab <- as_nu_lab(rs$nu)
fig1_df <- rbind(
  data.frame(dropout = rs$dropout, phase = rs$phase, nu_lab = rs$nu_lab,
             stat = "Median", value = rs$med_rel_err),
  data.frame(dropout = rs$dropout, phase = rs$phase, nu_lab = rs$nu_lab,
             stat = "Maximum", value = rs$max_rel_err))
fig1_df$stat  <- factor(fig1_df$stat, levels = c("Maximum", "Median"))
fig1_df$value <- pmax(fig1_df$value, 1e-16)

fig1 <- ggplot(fig1_df, aes(dropout, value, colour = stat, group = stat)) +
  geom_line() + geom_point(size = 2) +
  facet_grid(nu_lab ~ phase, labeller = nu_labeller) +
  scale_y_log10() +
  scale_colour_manual(values = c("Median" = "#0072B2",
                                 "Maximum" = "#D55E00")) +
  labs(x = "Dropout proportion",
       y = "Reconstruction error relative to SE") +
  base_theme
save_fig(fig1, "fig1_reconstruction_error", height = 7.2)

# ============================================================
# Figure 2: agreement with the NB MLE (QP vs POISE)
# ============================================================
os <- oc$summary
os$nu_f      <- as_nu_f(os$nu)
os$method_f  <- as_method(os$method)
os$alloc_f   <- as_alloc(os$alloc)
os$time_plot <- pmax(os$time_per_fit, 1e-7)
os_alt  <- os[os$effect == "Alternative", ]
os_null <- os[os$effect == "Null", ]

fig2_df <- os_alt[os_alt$method %in% c("QP", "OneStep"), ]
fig2 <- ggplot(fig2_df, aes(nu_f, mad_vs_NB,
                            colour = method_f, group = method_f)) +
  geom_line() + geom_point(size = 2.6) +
  facet_grid(alloc_f ~ phase) +
  scale_y_log10() +
  scale_colour_manual(values = METHOD_COLORS) +
  labs(x = expression(nu),
       y = "Median |estimate - NB MLE| (log rate ratio)") +
  base_theme
save_fig(fig2, "fig2_agreement_with_nb")

# ============================================================
# Figure 3: confidence interval coverage
# ============================================================
fig3 <- ggplot(os_alt, aes(nu_f, coverage, colour = method_f, group = method_f)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
  geom_line() + geom_point(size = 2.6) +
  facet_grid(alloc_f ~ phase) +
  scale_colour_manual(values = METHOD_COLORS) +
  coord_cartesian(ylim = c(0.90, 0.99)) +
  labs(x = expression(nu),
       y = "95% Wald interval coverage") +
  base_theme
save_fig(fig3, "fig3_coverage")

# ============================================================
# Figure 4: computation time per fit
# ============================================================
fig4 <- ggplot(os_alt, aes(nu_f, time_plot, colour = method_f, group = method_f)) +
  geom_line() + geom_point(size = 2.6) +
  facet_grid(alloc_f ~ phase) +
  scale_y_log10() +
  scale_colour_manual(values = METHOD_COLORS) +
  labs(x = expression(nu),
       y = "Elapsed time per fit (seconds)") +
  base_theme
save_fig(fig4, "fig4_computation_time")

# ============================================================
# Figure 5: type I error under the null (rate ratio = 1)
# ============================================================
os_null$type_I <- 1 - os_null$coverage
fig5 <- ggplot(os_null, aes(nu_f, type_I, colour = method_f, group = method_f)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey50") +
  geom_line() + geom_point(size = 2.6) +
  facet_grid(alloc_f ~ phase) +
  scale_colour_manual(values = METHOD_COLORS) +
  coord_cartesian(ylim = c(0, 0.10)) +
  labs(x = expression(nu),
       y = "Type I error (two-sided, nominal 0.05)") +
  base_theme
save_fig(fig5, "fig5_type_I_error")

# ============================================================
# LaTeX table helpers
# ============================================================
fmt_sci <- function(x) {
  vapply(x, function(v) {
    if (!is.finite(v)) return("---")
    if (v == 0) return("$0$")
    e <- floor(log10(abs(v)))
    m <- v / 10^e
    if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("$%.1f \\times 10^{%d}$", m, e)
  }, character(1))
}
fmt_f <- function(x, d) formatC(x, format = "f", digits = d)

write_latex_table <- function(cells, header, colspec, caption, label, file) {
  body <- apply(cells, 1L, function(r) paste0(paste(r, collapse = " & "), " \\\\"))
  lines <- c(
    "\\documentclass{article}",
    "\\usepackage{booktabs}",
    "\\usepackage[margin=1in]{geometry}",
    "\\begin{document}",
    "",
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{%s}", colspec),
    "\\toprule",
    paste0(paste(header, collapse = " & "), " \\\\"),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\end{document}")
  writeLines(lines, file)
}

# ============================================================
# Table 1: reconstruction fidelity under constant follow-up
# ============================================================
t1 <- rs[rs$dropout == 0, ]
t1 <- t1[order(t1$phase, t1$nu), ]
t1_cells <- cbind(
  t1$phase,
  formatC(t1$nu, format = "g"),
  as.character(t1$n_C),
  fmt_sci(t1$med_abs_err_log_RR),
  fmt_sci(t1$max_abs_err_log_RR),
  fmt_sci(t1$med_abs_err_SE),
  fmt_sci(t1$med_rel_err))
write_latex_table(
  cells   = t1_cells,
  header  = c("Phase", "$\\nu$", "$n$/group",
              "med $|\\Delta \\log RR|$", "max $|\\Delta \\log RR|$",
              "med $|\\Delta SE|$", "med rel.\\ err."),
  colspec = "lrrrrrr",
  caption = paste("Reconstruction fidelity under constant within-group",
                  "follow-up (no dropout). The summary-based estimator",
                  "reproduces the negative binomial maximum likelihood fit",
                  "to numerical precision."),
  label   = "tab:recon-fidelity",
  file    = file.path(results_dir, "table1_reconstruction_fidelity.tex"))

# ============================================================
# Table 2: operating characteristics summary
# ============================================================
t2 <- os_alt[order(os_alt$phase, os_alt$nu, os_alt$alloc,
                   match(os_alt$method, METHOD_LEVELS)), ]
t2_cells <- cbind(
  t2$phase,
  formatC(t2$nu, format = "g"),
  ifelse(t2$alloc == 1, "1:1", "2:1"),
  METHOD_LABELS[match(t2$method, METHOD_LEVELS)],
  fmt_f(t2$bias, 4),
  fmt_f(t2$rmse, 4),
  ifelse(t2$method == "NB", "---", fmt_sci(t2$mad_vs_NB)),
  fmt_f(t2$coverage, 3),
  fmt_f(t2$time_per_fit, 5))
write_latex_table(
  cells   = t2_cells,
  header  = c("Phase", "$\\nu$", "Alloc.", "Method", "Bias", "RMSE",
              "med.\\ vs NB", "Coverage", "Time (s)"),
  colspec = "lrllrrrrr",
  caption = paste("Operating characteristics by scenario and method.",
                  "Bias and RMSE are relative to the true log rate ratio;",
                  "med.\\ vs NB is the median absolute deviation from the",
                  "negative binomial maximum likelihood estimate; Time is",
                  "the elapsed time per fit."),
  label   = "tab:oc-summary",
  file    = file.path(results_dir, "table2_operating_characteristics.tex"))

cat("Figures and tables written to:\n  ", results_dir, "\n")
