#!/usr/bin/env Rscript
# Per-gene variant counts and frequencies in cases and gnomAD controls, for genes without suspicious indel pairs.
# Usage: Rscript gene_variant_frequencies.R run_dir dataset filter_string suspicious_genes_file [n_top]

options(width=250)
suppressMessages(library(data.table))
args <- commandArgs(TRUE)
run_dir   <- args[1]; ds <- args[2]; filt <- args[3]; suspicious_file <- args[4]
n_top     <- if (length(args) >= 5) as.integer(args[5]) else 50
id_col    <- 1   # variant ID column in the variant lists
gene_col  <- 4   # gene column in the variant lists

in_run <- function(f) file.path(run_dir, f)

#summary_df <- data.table::fread(in_run(paste0("rvexcaliber_testing_", ds, "_SummaryAssociations_", filt, "_iCFgCFadjust_rvexcaliber_base.txt")))
summary_df <- data.table::fread(in_run(paste0("rvexcaliber_testing_", ds, "_SummaryAssociations_allele_filter_", filt, "_iCFgCFadjust_rvexcaliber_base.txt")))
data.table::setnames(summary_df, 1:4, c("Gene", "observed", "expected", "p_rvx"))

int_mat <- data.table::fread(in_run(paste0(ds, "_RVBurdenMatrix_", filt, ".txt.gz")))
gno_mat <- data.table::fread(in_run(paste0("gnomAD_RVBurdenMatrix_", ds, "_", filt, ".txt.gz")))
n_cases <- ncol(int_mat) - 1

int_counts <- data.table::data.table(Gene          = int_mat$Gene,
                                     alleles_cases  = rowSums(as.matrix(int_mat[, -1])),
                                     carriers_cases = rowSums(as.matrix(int_mat[, -1]) > 0))
gno_counts <- data.table::data.table(Gene = gno_mat$Gene, expected_raw = rowSums(as.matrix(gno_mat[, -1])))

count_variants <- function(file, label) {
  vl <- data.table::fread(in_run(file), header = FALSE)
  vl[, .(n = data.table::uniqueN(get(paste0("V", id_col)))), by = c(Gene = paste0("V", gene_col))] |>
    data.table::setnames(c("Gene", label))
}
int_vars <- count_variants(paste0(ds, "_variant_list_extracted_", filt, ".txt.gz"), "n_variants_cases")
gno_vars <- count_variants(paste0("gnomAD_variant_list_extracted_", filt, ".txt.gz"), "n_variants_controls")

suspicious <- unique(data.table::fread(suspicious_file, header = FALSE)[[1]])

res <- Reduce(function(a, b) merge(a, b, by = "Gene", all.x = TRUE),
              list(summary_df[, .(Gene, observed, expected, p_rvx)], int_counts, gno_counts, int_vars, gno_vars))
res[is.na(n_variants_cases),    n_variants_cases := 0]
res[is.na(n_variants_controls), n_variants_controls := 0]
res[, freq_cases    := alleles_cases / (2 * n_cases)]
res[, freq_controls := expected_raw  / (2 * n_cases)]

res <- res[!Gene %in% suspicious][order(p_rvx)]
out <- res[seq_len(min(n_top, .N)),
           .(Gene, p_rvx = signif(p_rvx, 3), observed, expected = round(expected, 2),
             n_variants_cases, carriers_cases, freq_cases = signif(freq_cases, 3),
             n_variants_controls, freq_controls = signif(freq_controls, 3))]

cat(sprintf("Cases: %d | genes after removing %d with suspicious indel pairs: %d\n\n",
            n_cases, length(suspicious), nrow(res)))
print(out)
out_file <- in_run(paste0("gene_variant_frequencies_", ds, "_", filt, ".tsv"))
data.table::fwrite(res, out_file, sep = "\t")
cat("\nFull table written to:", out_file, "\n")
