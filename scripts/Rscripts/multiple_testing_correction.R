#!/usr/bin/env Rscript
#=======================================================================================================================
# multiple_testing_correction.R
#
# Adds multiple-testing corrections to an RV-EXCALIBER gene-level summary file:
#   p_bonferroni : Bonferroni-adjusted p-value (p x number of genes tested, capped at 1)
#                  A gene is significant if it's below 0.05,
#                  which is the same as comparing the raw p-value with 0.05/1,131 = 4.4 × 10^-5.
#   q_fdr_bh     : Benjamini-Hochberg false discovery rate q-value
#                  It's less strict than Bonferroni, and usual for exploratory gene lists:
#                  q < 0.10 means you expect fewer than 10% of genes below that cutoff to be false positives.
#
# Correction is across the genes actually tested, i.e. the rows of the summary file supplied (normally the
# ..._SummaryAssociations_allele_filter_..._iCFgCFadjust_rvexcaliber_base.txt file).
#
# Usage:
#   Rscript multiple_testing_correction.R <summary_file> [out_prefix] [alpha] [n_show]
#     out_prefix : output path prefix (default: <summary_dir>/SummaryAssociations_multiple_testing)
#     alpha      : significance level for Bonferroni and FDR (default 0.05)
#     n_show     : number of top genes to print (default 20)
#
# Output: <out_prefix>.tsv  (all genes, sorted by p-value, with p_bonferroni and q_fdr_bh added)
#
# This script outputs the following to the screen - a table of the top genes (20 by default), with these columns:
#      Column           Meaning
#      ======           =======
#      Gene             gene symbol
#      observed         qualifying alleles in your cohort
#      expected         expected alleles from gnomAD (iCF/gCF-adjusted)
#      ratio            observed / expected
#      p_rvx            RV-EXCALIBER p-value
#      p_bonferroni     Bonferroni-adjusted p-value (significant if < 0.05)
#      q_fdr_bh         Benjamini–Hochberg FDR q-value
#=======================================================================================================================

Sys.setenv(TZ = "Australia/Sydney")
options(width = 150, error = function() { traceback(2); quit(status = 1, save = "no") })
suppressMessages(library(data.table))

args       <- commandArgs(TRUE)
if (length(args) < 1) stop("Usage: Rscript multiple_testing_correction.R <summary_file> [out_prefix] [alpha] [n_show]")
sum_file   <- args[1]
out_prefix <- if (length(args) >= 2) args[2] else file.path(dirname(sum_file), "SummaryAssociations_multiple_testing")
alpha      <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
n_show     <- if (length(args) >= 4) as.integer(args[4]) else 20L

#-----------------------------------------------------------------------------------------------------------------------
# Read and correct
#-----------------------------------------------------------------------------------------------------------------------

res <- data.table::fread(sum_file)
data.table::setnames(res, 1:4, c("Gene", "observed", "expected", "p_rvx"))

n_missing_p <- sum(is.na(res$p_rvx))
res <- res[!is.na(p_rvx)]
n   <- nrow(res)

res[, `:=`(ratio        = observed / expected,
           p_bonferroni = stats::p.adjust(p_rvx, method = "bonferroni"),
           q_fdr_bh     = stats::p.adjust(p_rvx, method = "BH"))]
res <- res[order(p_rvx)]

#-----------------------------------------------------------------------------------------------------------------------
# Report
#-----------------------------------------------------------------------------------------------------------------------

cat(sprintf("Summary file: %s\n", basename(sum_file)))
cat(sprintf("Genes tested: %d%s\n", n, if (n_missing_p > 0) sprintf(" (%d rows with missing p excluded)", n_missing_p) else ""))
cat(sprintf("Bonferroni threshold (alpha = %.2f): %.3g\n", alpha, alpha / n))
cat(sprintf("Bonferroni-significant genes: %d\n", sum(res$p_bonferroni < alpha)))
cat(sprintf("FDR < %.2f: %d | FDR < 0.10: %d\n", alpha, sum(res$q_fdr_bh < alpha), sum(res$q_fdr_bh < 0.10)))
cat(sprintf("Smallest p-value: %.3g (%s)\n\n", res$p_rvx[1], res$Gene[1]))

print(res[seq_len(min(n_show, .N)),
          .(Gene, observed, expected = round(expected, 2), ratio = round(ratio, 2),
            p_rvx = signif(p_rvx, 3), p_bonferroni = signif(p_bonferroni, 3), q_fdr_bh = signif(q_fdr_bh, 3))])

#-----------------------------------------------------------------------------------------------------------------------
# Write
#-----------------------------------------------------------------------------------------------------------------------

out_file <- paste0(out_prefix, ".tsv")
data.table::fwrite(res, out_file, sep = "\t")
cat("\nWritten:", out_file, "\n")

