#!/usr/bin/env Rscript
#=======================================================================================================================
# report_candidate_genes.R
#
# Builds a candidate-gene table from an RV-EXCALIBER run, with QC columns so artefacts are excluded by
# stated rules rather than case by case. Genes passing QC are labelled "significant" (Bonferroni-adjusted
# p < 0.05) or "suggestive" (nominal p < 0.05); genes failing QC are labelled "excluded_by_QC" regardless
# of their p-value.
#
# Outputs (in out_dir):
#   candidate_genes_all.tsv        every tested gene with estimates and QC columns
#   candidate_genes.tsv genes meeting the suggestive criteria and passing QC
#   candidate_genes_panel.tsv      pre-specified candidate genes, reported whatever their p-value
#
# In this script, artefact_regex flags gene families known to produce mapping artefacts.
#
#=======================================================================================================================

Sys.setenv(TZ = "Australia/Sydney")
options(width=150)

#-----------------------------------------------------------------------------------------------------------------------
# Settings: edit these
#-----------------------------------------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
run_dir = as.character(args[1])
ds = as.character(args[2])
filt = as.character(args[3])

#run_dir  <- "/srv/scratch/z3531501/multiple_primaries/data/vcf_hg19_to_plink_then_rvexcaliber_run14"
#ds       <- "ccvANDmpvardb.sqlite_to_vcf.sort_normalise.annovar.PASS.liftover_to_hg19.leftaligned.ided.filter_DPge10.5_samples_removed"
#filt     <- "0.001_0.025_nfe90_amr7_eas3_rcc"

out_dir  <- file.path(run_dir, "candidate_gene_report")

summary_file <- file.path(run_dir, paste0("rvexcaliber_testing_", ds, "_SummaryAssociations_allele_filter_", filt,
                                          "_iCFgCFadjust_rvexcaliber_base.txt"))
matrix_file  <- file.path(run_dir, paste0(ds, "_RVBurdenMatrix_", filt, ".txt.gz"))

# suggestive criteria (fix these before looking at results)
p_suggestive     <- 0.01   # nominal p-value threshold
top_n_fallback   <- 20     # always report at least the top N genes by p-value
min_observed     <- 3      # at least this many internal alleles

# QC thresholds
max_frac_hom        <- 0.5   # flag if more than half of carriers have >= 2 alleles in the gene
max_cohort_ratio    <- 4     # flag if per-sample allele rate differs > 4-fold between cohorts (either direction)
ccv_sample_regex    <- "^CCV"

# artefact-prone gene families (segmental duplications, VCNTRs, pseudogenes, highly polymorphic loci)
# artefact_regex <- "^(MUC[0-9]|USP17L|NOMO[0-9]|KRT[0-9]|PRSS[0-9]|CR1$|DRD4$|ACAN$|HLA-|OR[0-9]+[A-Z]|ANKRD36|FRG[0-9]|NBPF|TBC1D3|GOLGA[0-9]|POTE|CTAGE)"
# artefact_regex <- "^(MUC[0-9]|USP17L|NOMO[0-9]|KRT[0-9]|PRSS[0-9]|CR1$|DRD4$|ACAN$|HLA-|OR[0-9]+[A-Z]|ANKRD36|FRG[0-9]|NBPF|TBC1D3|GOLGA[0-9]|POTE|CTAGE|NPIP|ZNF717$)"
artefact_regex <- "^(MUC[0-9]|USP17L|NOMO[0-9]|KRT[0-9]|PRSS[0-9]|CR1$|DRD4$|ACAN$|CEL$|ENOSF1$|HLA-|OR[0-9]+[A-Z]|ANKRD36|FRG[0-9]|NBPF|TBC1D3|GOLGA[0-9]|POTE|CTAGE|NPIP|ZNF717$)"

# pre-specified candidate genes (edit to your panel; define before looking at results)
panel_genes <- c("ATM", "BRCA1", "BRCA2", "PALB2", "CHEK2", "BRIP1", "RAD51C", "RAD51D", "BARD1",
                 "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "TP53", "PTEN", "CDH1", "STK11",
                 "APC", "MUTYH", "CDKN2A", "BAP1", "POLE", "POLD1", "NBN", "SMAD4", "BMPR1A")

#-----------------------------------------------------------------------------------------------------------------------
# Read inputs
#-----------------------------------------------------------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sum_df <- readr::read_tsv(summary_file, show_col_types = FALSE)
colnames(sum_df)[1:4] <- c("Gene", "observed", "expected", "p")

burden <- data.table::fread(matrix_file)
sample_cols <- setdiff(colnames(burden), "Gene")
is_ccv      <- grepl(ccv_sample_regex, sample_cols)
bm          <- as.matrix(burden[, ..sample_cols])

carrier_qc <- tibble::tibble(
  Gene             = burden$Gene,
  n_carriers       = rowSums(bm > 0),
  n_carriers_ge2   = rowSums(bm >= 2),
  alleles_ccv      = rowSums(bm[, is_ccv, drop = FALSE]),
  alleles_mpvardb  = rowSums(bm[, !is_ccv, drop = FALSE]),
  rate_ccv         = rowSums(bm[, is_ccv, drop = FALSE]) / sum(is_ccv),
  rate_mpvardb     = rowSums(bm[, !is_ccv, drop = FALSE]) / sum(!is_ccv)
)

#-----------------------------------------------------------------------------------------------------------------------
# Estimates, FDR and QC flags
#-----------------------------------------------------------------------------------------------------------------------

ratio_ci <- function(o, e) {
  ci <- stats::poisson.test(o, T = e)$conf.int
  c(ci[1], ci[2])
}

all_df <- sum_df |>
  dplyr::mutate(
    ratio_obs_exp = observed / expected,
    q_BH          = stats::p.adjust(p, method = "BH"),
    rank_p        = rank(p, ties.method = "min")
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    ratio_ci_low  = ratio_ci(observed, expected)[1],
    ratio_ci_high = ratio_ci(observed, expected)[2]
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(carrier_qc, by = "Gene") |>
  dplyr::mutate(
    frac_hom_carriers  = dplyr::if_else(n_carriers > 0, n_carriers_ge2 / n_carriers, NA_real_),
    cohort_rate_ratio  = (rate_ccv + 1e-6) / (rate_mpvardb + 1e-6),
    flag_artefact_gene = grepl(artefact_regex, Gene),
    flag_many_hom      = !is.na(frac_hom_carriers) & frac_hom_carriers > max_frac_hom,
    flag_cohort_imbal  = cohort_rate_ratio > max_cohort_ratio | cohort_rate_ratio < 1 / max_cohort_ratio,
    qc_pass            = !flag_artefact_gene & !flag_many_hom & !flag_cohort_imbal
  ) |>
  dplyr::arrange(p)

n_rows        <- nrow(all_df)
n_with_p      <- sum(!is.na(all_df$p))
n_p_exactly_1 <- sum(all_df$p == 1, na.rm = TRUE)

cat("rows in all_df:", n_rows,
    "| rows with a p-value:", n_with_p,
    "| rows with p == 1:", n_p_exactly_1, "\n")

# if the script has a tested/status column, cross-tabulate it against p being present
if ("tested" %in% names(all_df)) print(table(tested = all_df$tested, has_p = !is.na(all_df$p)))

#-----------------------------------------------------------------------------------------------------------------------
# Multiple-testing correction, candidate list (significant / suggestive / excluded_by_QC) and pre-specified panel
#-----------------------------------------------------------------------------------------------------------------------

alpha_mtc <- 0.05                                     # family-wise alpha for Bonferroni

all_df <- all_df |>
  dplyr::mutate(
    p_bonferroni = stats::p.adjust(p, method = "bonferroni"),
    q_fdr_bh     = stats::p.adjust(p, method = "BH")
  )

candidate_df <- all_df |>
  dplyr::filter(observed > expected,
                p_bonferroni < alpha_mtc |
                  (observed >= min_observed & (p < p_suggestive | rank_p <= top_n_fallback))) |>
  dplyr::mutate(report_status = dplyr::case_when(
    !qc_pass                    ~ "excluded_by_QC",
    p_bonferroni < alpha_mtc    ~ "significant",
    TRUE                        ~ "suggestive"
  )) |>
  dplyr::arrange(p)

panel_df <- all_df |>
  dplyr::filter(Gene %in% panel_genes)

missing_panel <- setdiff(panel_genes, all_df$Gene)

#-----------------------------------------------------------------------------------------------------------------------
# Write and summarise
#-----------------------------------------------------------------------------------------------------------------------

readr::write_tsv(all_df,        file.path(out_dir, "candidate_genes_all.tsv"))
readr::write_tsv(candidate_df, file.path(out_dir, "candidate_genes.tsv"))
readr::write_tsv(panel_df,      file.path(out_dir, "candidate_genes_panel.tsv"))

cat("Genes tested:", nrow(all_df), "| Bonferroni threshold:", signif(0.05 / nrow(all_df), 3),
    "| min q_BH:", signif(min(all_df$q_BH), 3), "\n")
cat("Suggestive (passing QC):", sum(candidate_df$report_status == "suggestive"),
    "| excluded by QC:", sum(candidate_df$report_status == "excluded_by_QC"), "\n")
cat("Panel genes tested:", nrow(panel_df), "| not tested (below allele filter or no qualifying variants):",
    paste(missing_panel, collapse = ", "), "\n\n")

print(as.data.frame(candidate_df |>
  dplyr::select(Gene, observed, expected, ratio_obs_exp, ratio_ci_low, ratio_ci_high, p, q_BH,
                n_carriers, frac_hom_carriers, cohort_rate_ratio, report_status)), digits = 3)

