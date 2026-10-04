#!/usr/bin/env Rscript
#=======================================================================================================================
# pool_indel_tract_AF.R
#
# Pools gnomAD allele frequencies across all FRAMESHIFT indels in the same repeat tract, so the gnomAD MAF filter
# treats every frameshift in a tract as one locus (the caller may place the same event anywhere in the tract).
# In-frame indels (length change a multiple of 3) are neither pooled nor used to build tracts: they keep their own
# AF, exactly like SNVs.
#
# Tract  = chain of overlapping frameshift-indel intervals
# Pooled = per population, sum of member AFs, capped at 1 (NA only if every member is NA)
#
# Usage:
#   Rscript pool_indel_tract_AF.R build_annovar <gnomad211_exome_indels.txt.gz> <tracts_out.txt.gz>   (recommended)
#   Rscript pool_indel_tract_AF.R build    <gnomAD_dir> <coverage> <tracts_out.txt.gz>   (pruned files: incomplete)
#   Rscript pool_indel_tract_AF.R annotate <tracts.txt.gz> <in_R_input.txt[.gz]> <out_R_input.txt.gz>
#
# 'annotate' appends AF_nfe_tract ... AF_amr_tract:
#   - SNVs, and indels not overlapping any gnomAD tract: own AF (unchanged)
#   - indels overlapping a gnomAD tract:                max(own AF, pooled tract AF)
# The original AF_* columns are left untouched; get_Varlist uses them for the gnomAD CMAC.
#=======================================================================================================================

options(error = function() { traceback(2); quit(status = 1, save = "no") })

suppressMessages(library(data.table))

args    <- commandArgs(TRUE)
mode    <- args[1]
af_cols <- c("AF_nfe", "AF_afr", "AF_sas", "AF_eas", "AF_amr")

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pool_af <- function(x) if (all(is.na(x))) NA_real_ else min(1, sum(x, na.rm = TRUE))
max_af  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

# Indel intervals in anchored (VCF-style) coordinates. Accepts VCF-style alleles (CCAG>C)
# or ANNOVAR-style alleles (CAG>- / ->CAG), which are shifted back onto the anchor base.
get_intervals <- function(chr, pos, ref, alt) {

  pos <- as.integer(pos); ref <- as.character(ref); alt <- as.character(alt)

  # RV-EXCALIBER renamed indels: Ref = "I"/"D", Alt = inserted/deleted length (optionally "_var_N"),
  # Pos = VCF anchor base
  id_ins  <- ref == "I"
  id_del  <- ref == "D"
  id_len  <- suppressWarnings(as.integer(sub("_var_.*$", "", alt)))
  id_code <- (id_ins | id_del) & !is.na(id_len)

  av_ins  <- !id_code & ref == "-"
  av_del  <- !id_code & alt == "-"
  vcf     <- !id_code & !av_ins & !av_del & nchar(ref) != nchar(alt)
  vcf_del <- vcf & nchar(ref) > nchar(alt)
  vcf_ins <- vcf & !vcf_del

  start <- rep(NA_integer_, length(pos)); end <- start

  start[vcf_del] <- pos[vcf_del];     end[vcf_del] <- pos[vcf_del] + nchar(ref[vcf_del]) - 1L
  start[vcf_ins] <- pos[vcf_ins];     end[vcf_ins] <- pos[vcf_ins] + 1L
  start[av_del]  <- pos[av_del] - 1L; end[av_del]  <- pos[av_del] + nchar(ref[av_del]) - 1L
  start[av_ins]  <- pos[av_ins];      end[av_ins]  <- pos[av_ins] + 1L

  d <- id_code & id_del; i <- id_code & id_ins
  start[d] <- pos[d]; end[d] <- pos[d] + id_len[d]
  start[i] <- pos[i]; end[i] <- pos[i] + 1L

  # only frameshift-length indels (length change not a multiple of 3) take part in tract pooling;
  # in-frame indels are treated individually, like SNVs
  indel_len <- ifelse(id_code, id_len,
               ifelse(av_ins, nchar(alt), ifelse(av_del, nchar(ref), abs(nchar(ref) - nchar(alt)))))
  is_fs     <- (id_code | av_ins | av_del | vcf) & (indel_len %% 3L != 0L)

  data.table::data.table(Chr = sub("^chr", "", chr), start = start, end = end,
                         is_indel = is_fs)
}

valid_mode <- mode %in% c("build", "build_annovar", "annotate")
if (!valid_mode) stop("First argument must be 'build', 'build_annovar' or 'annotate'")

# Collapse a table of indels (Chr, Pos, Ref, Alt, AF_*) into tracts with pooled AFs
make_tracts <- function(dt) {

  iv <- cbind(get_intervals(dt$Chr, dt$Pos, dt$Ref, dt$Alt), dt[, ..af_cols])
  iv <- iv[is_indel == TRUE]
  for (a in af_cols) data.table::set(iv, j = a, value = to_num(iv[[a]]))

  # a new tract starts whenever an indel begins after every earlier indel on the chromosome has ended
  data.table::setorder(iv, Chr, start, end)
  iv[, prev_max_end := data.table::shift(cummax(end), fill = -1L), by = Chr]
  iv[, tract := cumsum(start > prev_max_end), by = Chr]

  tr <- iv[, c(list(start = min(start), end = max(end), n_indels = .N), lapply(.SD, pool_af)),
           by = .(Chr, tract), .SDcols = af_cols]

  cat("  indels:", nrow(iv), "| tracts:", nrow(tr),
      "| tracts with >1 indel:", sum(tr$n_indels > 1L), "\n"); flush.console()
  cat("  tract size distribution (5 = 5+):\n"); print(table(pmin(tr$n_indels, 5L))); flush.console()

  tr[, tract := NULL]
}


#-----------------------------------------------------------------------------------------------------------------------
# build_annovar: tracts from the full ANNOVAR gnomAD 2.1.1 exome database (indel rows)
#   Rscript pool_indel_tract_AF.R build_annovar <gnomad211_exome_indels.txt.gz> <tracts_out.txt.gz>
#-----------------------------------------------------------------------------------------------------------------------

if (mode == "build_annovar") {

  in_file <- args[2]; out_file <- args[3]

  if (is.na(out_file) || !nzchar(out_file)) stop("No output file given (3rd argument after 'build_annovar')")

  cat("Reading", in_file, "\n"); flush.console()
  dt <- data.table::fread(in_file, sep = "\t", colClasses = "character")
  data.table::setnames(dt, c("#Chr", "Start"), c("Chr", "Pos"), skip_absent = TRUE)

  missing_cols <- setdiff(c("Chr", "Pos", "Ref", "Alt", af_cols), colnames(dt))
  if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  tracts <- make_tracts(dt)
  data.table::fwrite(tracts, out_file, sep = "\t", na = "NA")
  cat("Tracts written:", nrow(tracts), "->", out_file, "\nDone.\n")
}


#-----------------------------------------------------------------------------------------------------------------------
# build: one row per gnomAD indel tract, with pooled AFs
#-----------------------------------------------------------------------------------------------------------------------

if (mode == "build") {

  gnomAD_dir <- args[2]; coverage <- args[3]; out_file <- args[4]
  tr_list    <- vector("list", 22)

  for (chr in 1:22) {

    f <- file.path(gnomAD_dir, paste0(chr, "_gnomAD_pruned_annotation_for_R_input_", coverage, ".txt.gz"))
    cat("Reading", f, "\n"); flush.console()

    dt <- data.table::fread(f, sep = "\t", select = c("Chr", "Pos", "Ref", "Alt", af_cols),
                            colClasses = "character")

    tr_list[[chr]] <- make_tracts(dt)
  }

  tracts <- data.table::rbindlist(tr_list)
  data.table::fwrite(tracts, out_file, sep = "\t", na = "NA")
  cat("Tracts written:", nrow(tracts), "->", out_file, "\nDone.\n")
}


#-----------------------------------------------------------------------------------------------------------------------
# annotate: add *_tract AF columns to an R-input file (gnomAD or internal)
#-----------------------------------------------------------------------------------------------------------------------

if (mode == "annotate") {

  tracts_file <- args[2]; in_file <- args[3]; out_file <- args[4]

  tracts <- data.table::fread(tracts_file, sep = "\t", colClasses = list(character = "Chr"))
  data.table::setkey(tracts, Chr, start, end)

  # read everything as character so untouched columns are written back exactly as they were
  dt <- data.table::fread(in_file, sep = "\t", colClasses = "character", na.strings = "")

  iv <- get_intervals(dt$Chr, dt$Pos, dt$Ref, dt$Alt)
  iv[, row_id := .I]
  q  <- iv[is_indel == TRUE, .(Chr, start, end, row_id)]

  ov  <- data.table::foverlaps(q, tracts, by.x = c("Chr", "start", "end"), type = "any", nomatch = NULL)
  hit <- ov[, lapply(.SD, max_af), by = row_id, .SDcols = af_cols]

  idx <- match(seq_len(nrow(dt)), hit$row_id)
  has <- !is.na(idx)

  for (a in af_cols) {
    own         <- to_num(dt[[a]])
    pooled      <- own
    pooled[has] <- pmax(own[has], hit[[a]][idx[has]], na.rm = TRUE)
    data.table::set(dt, j = paste0(a, "_tract"),
                    value = ifelse(is.na(pooled), ".", as.character(pooled)))
  }

  cat("Rows:", nrow(dt),
      "| indels:", nrow(q),
      "| indels overlapping a gnomAD tract:", sum(has), "\n"); flush.console()

  data.table::fwrite(dt, out_file, sep = "\t", quote = FALSE, na = "")
  cat("Written:", out_file, "\nDone.\n")
}
