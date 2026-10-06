#!/usr/bin/env Rscript
#=======================================================================================================================
# panel_burden_test.R
#
# Pooled (gene-set) burden test of a pre-specified cancer predisposition gene panel, built directly from the
# RV-EXCALIBER internal and gnomAD burden matrices, so every panel gene is included regardless of RV-EXCALIBER's
# per-gene filters (keep_genes, MIGen ranking, allele-count threshold).
#
#   observed  = qualifying alleles in the cohort (row sum of the internal burden matrix; 0 if the gene is absent)
#   expected  = raw gnomAD expected alleles (row sum of the gnomAD burden matrix) x adjustment factor, where the
#               factor is the median ratio of RV-EXCALIBER's iCF/gCF-adjusted expected to raw expected
#               across all genes in the summary file
#   test      = exact Poisson test of pooled observed vs pooled expected (one-sided, greater), plus mid-p
#   baseline  = rate ratio of panel genes vs all other genes, which absorbs any genome-wide sensitivity shift
#
# Usage: Rscript panel_burden_test.R [run_dir] [internal_dataset] [filter_string] [panel_file]
#        panel_file: optional, one gene symbol per line (e.g. Huang et al. 2018 Cell Table S1); default built-in list
#=======================================================================================================================
 
Sys.setenv(TZ = "Australia/Sydney")
options(width = 150)
suppressMessages(library(data.table))
 
#-----------------------------------------------------------------------------------------------------------------------
# Settings
#-----------------------------------------------------------------------------------------------------------------------
 
args    <- commandArgs(TRUE)
run_dir <- if (length(args) >= 1) args[1] else "/srv/scratch/z3531501/multiple_primaries/data/vcf_hg19_to_plink_then_amrevel_rvexcaliber_run21.DPge10.rcc"
ds      <- if (length(args) >= 2) args[2] else "ccvANDmpvardb.sqlite_to_vcf.sort_normalise.annovar.PASS.liftover_to_hg19.leftaligned.ided.filter_DPge10.5_samples_removed"
filt    <- if (length(args) >= 3) args[3] else "0.001_0.5_nfe90_amr7_eas3_rcc"
use_panel_file <- length(args) >= 4
 
builtin_panel <- c("ATM", "BRCA1", "BRCA2", "PALB2", "CHEK2", "BRIP1", "RAD51C", "RAD51D", "BARD1",
                   "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "TP53", "PTEN", "CDH1", "STK11",
                   "APC", "MUTYH", "CDKN2A", "BAP1", "POLE", "POLD1", "NBN", "SMAD4", "BMPR1A")
panel_genes <- if (use_panel_file) unique(trimws(data.table::fread(args[4], header = FALSE)[[1]])) else builtin_panel
 
out_dir <- file.path(run_dir, "panel_burden_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
#-----------------------------------------------------------------------------------------------------------------------
# Observed and expected per gene from the burden matrices
#-----------------------------------------------------------------------------------------------------------------------
 
int_mat <- data.table::fread(file.path(run_dir, paste0(ds, "_RVBurdenMatrix_", filt, ".txt.gz")))
gno_mat <- data.table::fread(file.path(run_dir, paste0("gnomAD_RVBurdenMatrix_", ds, "_", filt, ".txt.gz")))
 
counts <- merge(
  data.table::data.table(Gene = int_mat$Gene, observed = rowSums(as.matrix(int_mat[, -1]))),
  data.table::data.table(Gene = gno_mat$Gene, exp_raw  = rowSums(as.matrix(gno_mat[, -1]))),
  by = "Gene", all = TRUE)
counts[is.na(observed), observed := 0]
counts[is.na(exp_raw),  exp_raw  := 0]
 
#-----------------------------------------------------------------------------------------------------------------------
# Adjustment factor: RV-EXCALIBER adjusted expected / raw expected (median over summary genes)
#-----------------------------------------------------------------------------------------------------------------------
 
sum_file <- file.path(run_dir, paste0("rvexcaliber_testing_", ds, "_SummaryAssociations_", filt,
                                      "_iCFgCFadjust_rvexcaliber_base.txt"))
sum_all <- data.table::fread(sum_file)
data.table::setnames(sum_all, 1:4, c("Gene", "observed_rvx", "expected_rvx", "p_rvx"))
 
adj_tbl <- merge(counts, sum_all[, .(Gene, expected_rvx)], by = "Gene")[exp_raw > 0]
adj     <- stats::median(adj_tbl$expected_rvx / adj_tbl$exp_raw)
counts[, expected := exp_raw * adj]
cat(sprintf("Adjustment factor (median adjusted/raw expected over %d genes): %.3f\n\n", nrow(adj_tbl), adj))
 
#-----------------------------------------------------------------------------------------------------------------------
# Panel table
#-----------------------------------------------------------------------------------------------------------------------
 
counts[, panel := Gene %in% panel_genes]
panel_tbl <- counts[panel == TRUE][order(-(observed - expected))]
panel_tbl <- merge(panel_tbl, sum_all[, .(Gene, p_rvx)], by = "Gene", all.x = TRUE)[order(-(observed - expected))]
 
missing_genes <- setdiff(panel_genes, counts$Gene)
cat("Panel genes:", length(panel_genes), "| in matrices:", nrow(panel_tbl),
    "| absent from both matrices:", if (length(missing_genes)) paste(missing_genes, collapse = ", ") else "none", "\n\n")
print(panel_tbl[, .(Gene, observed, exp_raw = round(exp_raw, 2), expected = round(expected, 2),
                    ratio = round(observed / expected, 2), p_rvx = signif(p_rvx, 3))])
 
#-----------------------------------------------------------------------------------------------------------------------
# Pooled tests
#-----------------------------------------------------------------------------------------------------------------------
 
pool_test <- function(tbl, label) {
  o  <- sum(tbl$observed); e <- sum(tbl$expected)
  pt <- stats::poisson.test(o, T = e, alternative = "greater")
  ci <- stats::poisson.test(o, T = e)$conf.int
  p_mid <- stats::ppois(o, e, lower.tail = FALSE) + 0.5 * stats::dpois(o, e)
  data.table::data.table(test = label, n_genes = nrow(tbl), observed = o, expected = round(e, 2),
                         ratio = round(o / e, 2), ci_low = round(ci[1], 2), ci_high = round(ci[2], 2),
                         p_one_sided = signif(pt$p.value, 3), p_mid = signif(p_mid, 3))
}
 
tests <- data.table::rbindlist(list(
  pool_test(panel_tbl,                                              "all panel genes"),
  pool_test(panel_tbl[exp_raw >= 1],                                "PRIMARY: raw expected >= 1"),
  pool_test(panel_tbl[exp_raw >= 0.5],                              "raw expected >= 0.5"),
  pool_test(panel_tbl[exp_raw >= 1 & Gene != "BRCA1"],              "expected >= 1, without BRCA1"),
  pool_test(panel_tbl[exp_raw >= 1 & Gene != "BRCA2"],              "expected >= 1, without BRCA2"),
  pool_test(panel_tbl[exp_raw >= 1 & Gene != "ATM"],                "expected >= 1, without ATM"),
  pool_test(panel_tbl[exp_raw >= 1 & !(Gene %in% c("BRCA1", "BRCA2", "ATM"))], "expected >= 1, without BRCA1/BRCA2/ATM")
))
cat("\nPooled panel tests (exact Poisson, one-sided for excess; CI two-sided):\n")
print(tests)
 
#-----------------------------------------------------------------------------------------------------------------------
# Panel vs rest-of-genome rate ratio (absorbs genome-wide sensitivity differences)
#-----------------------------------------------------------------------------------------------------------------------
 
rr_tbl <- counts[exp_raw >= 1, .(observed = sum(observed), expected = sum(expected)), by = panel]
rr <- stats::poisson.test(x = c(rr_tbl[panel == TRUE, observed], rr_tbl[panel == FALSE, observed]),
                          T = c(rr_tbl[panel == TRUE, expected], rr_tbl[panel == FALSE, expected]),
                          alternative = "greater")
rr2 <- stats::poisson.test(x = c(rr_tbl[panel == TRUE, observed], rr_tbl[panel == FALSE, observed]),
                           T = c(rr_tbl[panel == TRUE, expected], rr_tbl[panel == FALSE, expected]))
cat(sprintf("\nPanel vs other genes (expected >= 1): panel O/E = %.2f, other O/E = %.2f, rate ratio = %.2f (95%% CI %.2f-%.2f), one-sided p = %.3g\n",
            rr_tbl[panel == TRUE, observed / expected], rr_tbl[panel == FALSE, observed / expected],
            rr$estimate, rr2$conf.int[1], rr2$conf.int[2], rr$p.value))
 
#-----------------------------------------------------------------------------------------------------------------------
# Write
#-----------------------------------------------------------------------------------------------------------------------
 
data.table::fwrite(panel_tbl, file.path(out_dir, "panel_genes.tsv"), sep = "\t")
data.table::fwrite(tests,     file.path(out_dir, "panel_tests.tsv"), sep = "\t")
cat("\nWritten to:", out_dir, "\n")
