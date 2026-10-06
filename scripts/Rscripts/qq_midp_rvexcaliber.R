#!/usr/bin/env Rscript
#=======================================================================================================================
# qq_midp_rvexcaliber.R
#
# QQ plots for an RV-EXCALIBER gene-level summary file, in three versions:
#
#   1. standard : RV-EXCALIBER's own p-values against the uniform diagonal (as in the pipeline's plot)
#   2. midp     : Poisson mid-p values against the uniform diagonal
#                 mid-p = P(X > observed) + 0.5 * P(X = observed),  X ~ Poisson(expected)
#   3. discrete : RV-EXCALIBER's own p-values against the null that is actually achievable for discrete counts,
#                 simulated by drawing observed ~ Poisson(expected) for every gene (n_sim times)
#
# Usage:
#   Rscript qq_midp_rvexcaliber.R <SummaryAssociations_..._allele_filter_...txt> [out_prefix] [type] [n_sim]
#   type: standard | midp | discrete | all   (default all)
#
# In this script, some of the RV-EXCALIBER columns have been renamed:
#      Original column                           Renamed to
#      ===============                           ==========
#      Gene                                      Gene
#      test_allele_count                         observed
#      test_gnomAD_allele_count_iCFgCFadjust     expected
#      test_P_rvexcaliber_base_iCFgCFadjust      p_rvx
#=======================================================================================================================

Sys.setenv(TZ = "Australia/Sydney")
options(error = function() { traceback(2); quit(status = 1, save = "no") })
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

args       <- commandArgs(TRUE)
sum_file   <- args[1]
out_prefix <- if (length(args) >= 2) args[2] else sub("\\.txt$", "", sum_file)
plot_type  <- if (length(args) >= 3) args[3] else "all"
n_sim      <- if (length(args) >= 4) as.integer(args[4]) else 1000L
set.seed(2026)

d <- data.table::fread(sum_file)
data.table::setnames(d, 1:4, c("Gene", "observed", "expected", "p_rvx"))
n <- nrow(d)

d[, p_std := stats::ppois(observed - 1, expected, lower.tail = FALSE)]                                 # P(X >= o)
d[, p_mid := stats::ppois(observed, expected, lower.tail = FALSE) + 0.5 * stats::dpois(observed, expected)]

lambda_median <- function(p) stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)

#-----------------------------------------------------------------------------------------------------------------------
# Plot helpers
#-----------------------------------------------------------------------------------------------------------------------

qq_theme <- ggplot2::theme_bw(base_size = 16)

uniform_qq <- function(p, title, subtitle) {
  p   <- sort(p)
  exp <- -log10(stats::ppoints(length(p)))
  i   <- seq_along(p)
  df  <- data.table::data.table(expected = exp, observed = -log10(p),
                                lo = -log10(stats::qbeta(0.975, i, length(p) - i + 1)),
                                hi = -log10(stats::qbeta(0.025, i, length(p) - i + 1)))
  ggplot2::ggplot(df, ggplot2::aes(expected, observed)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::labs(x = expression(Expected ~ -log[10](p) ~ "(uniform)"), y = expression(Observed ~ -log[10](p)),
                  title = title, subtitle = subtitle) +
    qq_theme
}

save_plot <- function(g, suffix) {
  f <- paste0(out_prefix, "_QQ_", suffix, ".png")
  ggplot2::ggsave(f, g, width = 8, height = 8, dpi = 150)
  cat("Written:", f, "\n")
}

do_plot <- function(type) plot_type %in% c(type, "all")

#-----------------------------------------------------------------------------------------------------------------------
# 1. Standard (RV-EXCALIBER p-values vs uniform)
#-----------------------------------------------------------------------------------------------------------------------

if (do_plot("standard")) {
  save_plot(uniform_qq(d$p_rvx, "RV-EXCALIBER p-values",
                       sprintf("n genes = %d | median lambda = %.2f", n, lambda_median(d$p_rvx))), "standard")
}

#-----------------------------------------------------------------------------------------------------------------------
# 2. Mid-p (Poisson mid-p vs uniform)
#-----------------------------------------------------------------------------------------------------------------------

if (do_plot("midp")) {
  save_plot(uniform_qq(d$p_mid, "Poisson mid-p values",
                       sprintf("n genes = %d | median lambda = %.2f (standard Poisson: %.2f)",
                               n, lambda_median(d$p_mid), lambda_median(d$p_std))), "midp")
}

#-----------------------------------------------------------------------------------------------------------------------
# 3. Discrete null (RV-EXCALIBER p-values vs simulated achievable null)
#    Null order statistics come from Poisson-standard p-values of simulated counts; RV-EXCALIBER's p-values rank
#    genes like the Poisson test (Spearman ~0.999), so the Poisson null is used as the reference shape.
#-----------------------------------------------------------------------------------------------------------------------

if (do_plot("discrete")) {
  sim <- vapply(seq_len(n_sim), function(k) {
    o_sim <- stats::rpois(n, d$expected)
    sort(-log10(stats::ppois(o_sim - 1, d$expected, lower.tail = FALSE)), decreasing = TRUE)
  }, numeric(n))
  null_df <- data.table::data.table(expected = rowMeans(sim),
                                    lo = apply(sim, 1, stats::quantile, 0.025),
                                    hi = apply(sim, 1, stats::quantile, 0.975),
                                    observed_rvx = sort(-log10(d$p_rvx), decreasing = TRUE),
                                    observed_std = sort(-log10(d$p_std), decreasing = TRUE))
  lam_null <- stats::median(apply(sim, 2, function(s) lambda_median(10^(-s))))

  g <- ggplot2::ggplot(null_df, ggplot2::aes(expected)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
    ggplot2::geom_point(ggplot2::aes(y = observed_std), size = 1.5) +
    ggplot2::labs(x = expression(Expected ~ -log[10](p) ~ "(simulated discrete null)"),
                  y = expression(Observed ~ -log[10](p) ~ "(Poisson)"),
                  title = "Poisson p-values vs simulated discrete null",
                  subtitle = sprintf("n genes = %d | %d simulations | observed lambda = %.2f, null lambda = %.2f",
                                     n, n_sim, lambda_median(d$p_std), lam_null)) +
    qq_theme
  save_plot(g, "discrete_null")
}

#-----------------------------------------------------------------------------------------------------------------------
# Table of top genes with all three p-values
#-----------------------------------------------------------------------------------------------------------------------

top <- d[order(p_mid)][1:min(20, .N), .(Gene, observed, expected = round(expected, 2),
                                         ratio = round(observed / expected, 2),
                                         p_rvx = signif(p_rvx, 3), p_poisson = signif(p_std, 3), p_mid = signif(p_mid, 3))]
cat("\nTop genes by mid-p (Bonferroni threshold:", signif(0.05 / n, 3), ")\n")
print(top)
data.table::fwrite(d[order(p_mid)], paste0(out_prefix, "_with_midp.tsv"), sep = "\t")
cat("Written:", paste0(out_prefix, "_with_midp.tsv"), "\n")
