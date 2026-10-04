args <- commandArgs(TRUE)
bim_file  <- args[1]
frq_file  <- args[2]
raw_file  <- args[3]
out_file  <- args[4]

options(error = function() {
  traceback(2)
  quit(status = 1, save = "no")
})

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

cat("Loading bim...\n"); flush.console()
bim <- fread(bim_file, header = FALSE, col.names = c("chrom", "id", "gdist", "pos", "A1", "A2"))
bim[, `:=`(REF = A2, ALT = A1)]
cat("bim loaded:", nrow(bim), "rows\n"); flush.console()

cat("Loading frq...\n"); flush.console()
frq <- fread(frq_file, header = TRUE)
cat("frq loaded:", nrow(frq), "rows\n"); flush.console()

is_indel <- (abs(nchar(bim$REF) - nchar(bim$ALT)) %% 3L) != 0L   # frameshift-length indels only; in-frame indels get own AF like SNVs
indel_bim <- bim[is_indel]
snp_ids   <- bim$id[!is_indel]

if (nrow(indel_bim) == 0) {

  cat("No indels found - using frq as-is\n"); flush.console()
  merged_af <- frq[, .(id = SNP, AF_int = MAF)]

} else {

  cat("Number of indels:", nrow(indel_bim), "\n"); flush.console()

  indel_bim[, `:=`(
    start = pos,
    end = ifelse(nchar(REF) > nchar(ALT), pos + nchar(REF) - 1L, pos + 1L)
  )]
  setorder(indel_bim, chrom, start)

  n <- nrow(indel_bim)
  parent <- seq_len(n)
  find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  union <- function(x, y) { rx <- find(x); ry <- find(y); if (rx != ry) parent[rx] <<- ry }

  for (chr in unique(indel_bim$chrom)) {
    idx <- which(indel_bim$chrom == chr)
    ord <- idx[order(indel_bim$start[idx])]
    max_end <- indel_bim$end[ord[1]]
    max_end_idx <- ord[1]
    if (length(ord) > 1) {
      for (k in 2:length(ord)) {
        i <- ord[k]
        if (indel_bim$start[i] <= max_end) {
          union(i, max_end_idx)
          max_end <- max(max_end, indel_bim$end[i])
        } else {
          max_end <- indel_bim$end[i]
        }
        max_end_idx <- i
      }
    }
  }

  indel_bim[, group_id := sapply(seq_len(n), find)]
  n_groups <- length(unique(indel_bim$group_id))
  cat("Number of indel groups:", n_groups, "\n"); flush.console()

  cat("Loading raw genotype file (this may take a while for large files)...\n"); flush.console()
  raw <- fread(raw_file, header = TRUE)
  cat("raw loaded:", nrow(raw), "samples x", ncol(raw), "columns\n"); flush.console()

  geno_cols   <- colnames(raw)[-(1:6)]
  geno_id_map <- sub("_[^_]+$", "", geno_cols)

  # keep only genotype columns that correspond to an indel (aligned to groups via match() below)
  keep <- geno_id_map %in% indel_bim$id
  col_ids <- geno_id_map[keep]

  # Confirm every indel in the bim has a genotype column in the raw file
  missing_ids <- setdiff(indel_bim$id, col_ids)
  cat("Indels in bim:", nrow(indel_bim),
      "| matched in raw:", length(unique(col_ids)),
      "| missing:", length(missing_ids), "\n"); flush.console()

  if (length(missing_ids) > 0) {
    cat("First missing IDs:", head(missing_ids, 5), "\n"); flush.console()
    stop("Indels in bim not found in raw file - check indel extraction step")
  }

  # geno_mat <- as.matrix(raw[, geno_cols[keep], with = FALSE])	# this can produce duplicate IDs, so instead do:
  geno_idx <- 6L + which(keep)
  geno_mat <- as.matrix(raw[, ..geno_idx])

  storage.mode(geno_mat) <- "double"

  # map each genotype column to its group_id, in the SAME order as geno_mat's columns
  col_group <- indel_bim$group_id[match(col_ids, indel_bim$id)]

  cat("Collapsing genotypes by group (vectorized)...\n"); flush.console()

  # sum genotype dosage across all columns sharing a group_id, per sample - vectorized via matrix %*% indicator
  groups_sorted <- sort(unique(col_group))

  # sparse indicator matrix: columns (variants) x groups, one 1 per row
  ind_mat <- Matrix::sparseMatrix(
    i    = seq_along(col_group),
    j    = match(col_group, groups_sorted),
    x    = 1,
    dims = c(ncol(geno_mat), length(groups_sorted))
  )

  # replace NA with 0 for summing dosage, but track missingness separately
  geno_mat_0 <- geno_mat
  geno_mat_0[is.na(geno_mat_0)] <- 0
  not_na_mat <- !is.na(geno_mat)

  dosage_sum   <- as.matrix(geno_mat_0 %*% ind_mat)   # samples x groups: summed dosage
  n_nonmissing <- as.matrix(not_na_mat %*% ind_mat)   # samples x groups: count of non-missing member variants per sample

  collapsed <- pmin(dosage_sum, 2)                # cap at 2 (adjust here if you prefer a different compound-het rule)
  collapsed[n_nonmissing == 0] <- NA               # sample missing at every member variant -> NA for the group

  # per-group merged AF_int
  group_AF <- colSums(collapsed, na.rm = TRUE) / (2 * colSums(!is.na(collapsed)))
  names(group_AF) <- groups_sorted

  indel_af_lookup <- indel_bim[, .(id, AF_int = group_AF[as.character(group_id)])]

  snp_af <- frq[SNP %in% snp_ids, .(id = SNP, AF_int = MAF)]

  merged_af <- rbindlist(list(snp_af, indel_af_lookup))
}

cat("Writing output...\n"); flush.console()

out <- merged_af %>%
  separate(id, into = c("Chr", "Pos", "Ref", "Alt"), sep = ":") %>%
  select(Chr, Pos, Ref, Alt, AF_int)

write.table(out, out_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Done.\n")
