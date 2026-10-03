#!/usr/bin/env Rscript
#=======================================================================================================================
# format_gnomad211_amrevel_pruned.R
#
# Converts one chromosome of ANNOVAR output (refGeneWithVer, gnomad211_exome, dbnsfp47a; -vcfinput) into the
# RV-EXCALIBER gnomAD control layout:
#   Chr Pos Ref Alt Alt2 Func.refGene Gene.refGene ExonicFunc.refGene AF_nfe AF_afr AF_sas AF_eas AF_amr pathogenic_missense
#
# pathogenic_missense = 1 if AlphaMissense_pred == "P" or REVEL_score >= 0.644, otherwise 0.
# Run RV-EXCALIBER with an MCAP threshold between 0 and 1 (e.g. 0.5).
#
# Usage: Rscript format_gnomad211_amrevel_pruned.R <chrN.hg19_multianno.txt> <chr> <out_dir> <hcc_bed>
#=======================================================================================================================

Sys.setenv(TZ = "Australia/Sydney")
options(error = function() { traceback(2); quit(status = 1, save = "no") })
suppressMessages(library(data.table))

args      <- commandArgs(TRUE)
in_file   <- args[1]; chr <- args[2]; out_dir <- args[3]; hcc_bed <- args[4]
revel_min <- 0.644
qual_exonic_func <- c("nonsynonymous_SNV", "stopgain", "stoploss", "startloss",
                      "frameshift_deletion", "frameshift_insertion")
splice_func <- c("splicing", "ncRNA_exonic;splicing")
pops      <- c("nfe", "afr", "sas", "eas", "amr")

d  <- data.table::fread(in_file, sep = "\t", colClasses = "character", quote = "", na.strings = NULL)
nm <- colnames(d)

pick <- function(pattern, what) {
  hit <- grep(pattern, nm, value = TRUE)
  if (length(hit) == 0) stop("No column for ", what, " (pattern '", pattern, "'). Columns: ", paste(nm, collapse = ", "))
  hit[1]
}

func_col  <- pick("^Func\\.refGeneWithVer$",       "Func")
gene_col  <- pick("^Gene\\.refGeneWithVer$",       "Gene")
exfun_col <- pick("^ExonicFunc\\.refGeneWithVer$", "ExonicFunc")
d[[exfun_col]] <- gsub(" ", "_", d[[exfun_col]], fixed = TRUE)
af_cols   <- vapply(pops, function(p) pick(paste0("^AF_", p, "$"), paste("AF", p)), "")
revel_col <- pick("^REVEL_score$",        "REVEL_score")
am_col    <- pick("^AlphaMissense_pred$", "AlphaMissense_pred")

# VCF-style coordinates and alleles from the ID written into the VCF (chr:pos:ref:alt), carried in Otherinfo
id_col <- NA_character_
for (cc in grep("^Otherinfo", nm, value = TRUE)) {
  if (is.na(id_col) && mean(grepl("^[0-9]+:[0-9]+:[ACGTN]+:[ACGTN]+$", head(d[[cc]], 1000))) > 0.9) id_col <- cc
}
if (is.na(id_col)) stop("Could not find the chr:pos:ref:alt ID among Otherinfo columns")
idp <- data.table::tstrsplit(d[[id_col]], ":", fixed = TRUE)

revel <- suppressWarnings(as.numeric(d[[revel_col]]))
flag  <- as.integer((!is.na(revel) & revel >= revel_min) | d[[am_col]] == "P")

out <- data.table::data.table(
  Chr                 = idp[[1]],
  Pos                 = idp[[2]],
  Ref                 = idp[[3]],
  Alt                 = idp[[4]],
  Alt2                = idp[[4]],
  Func.refGene        = d[[func_col]],
  Gene.refGene        = d[[gene_col]],
  ExonicFunc.refGene  = d[[exfun_col]],
  AF_nfe              = d[[af_cols["nfe"]]],
  AF_afr              = d[[af_cols["afr"]]],
  AF_sas              = d[[af_cols["sas"]]],
  AF_eas              = d[[af_cols["eas"]]],
  AF_amr              = d[[af_cols["amr"]]],
  pathogenic_missense = flag,
  annovar_start       = as.integer(d[["Start"]])
)
out <- out[(Func.refGene %in% c("exonic", "exonic;splicing") & ExonicFunc.refGene %in% qual_exonic_func) | Func.refGene %in% splice_func]

# hcc = sites inside high-coverage coding regions (ANNOVAR start, so deletions use the first deleted base)
bed <- data.table::fread(hcc_bed, header = FALSE, select = 1:3, col.names = c("Chr", "s0", "e"))
bed <- bed[as.character(Chr) == chr, .(Chr = as.character(Chr), start = as.integer(s0) + 1L, end = as.integer(e))]
data.table::setkey(bed, Chr, start, end)
pts <- out[, .(Chr, start = annovar_start, end = annovar_start, row = .I)]
in_hcc <- unique(data.table::foverlaps(pts, bed, nomatch = NULL)$row)

cols_out <- setdiff(colnames(out), "annovar_start")
f_rcc <- file.path(out_dir, paste0(chr, "_gnomAD_pruned_annotation_for_R_input_rcc.txt.gz"))
f_hcc <- file.path(out_dir, paste0(chr, "_gnomAD_pruned_annotation_for_R_input_hcc.txt.gz"))
data.table::fwrite(out[, ..cols_out],       f_rcc, sep = "\t", quote = FALSE)
data.table::fwrite(out[in_hcc, ..cols_out], f_hcc, sep = "\t", quote = FALSE)

cat(sprintf("chr%s: %d annotated | %d kept (rcc) | %d hcc | pathogenic_missense=1: %d (%.1f%% of nonsynonymous)\n",
            chr, nrow(d), nrow(out), length(in_hcc), sum(out$pathogenic_missense),
            100 * mean(out[ExonicFunc.refGene == "nonsynonymous_SNV", pathogenic_missense])))
cat("Columns used:", func_col, gene_col, exfun_col, paste(af_cols, collapse = ","), revel_col, am_col, id_col, "\n")
