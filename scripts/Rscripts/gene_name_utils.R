# Split an ANNOVAR Gene.refGene string into its genes.
# "A2M;A2M" -> "A2M"; "GENE1;GENE2" -> c("GENE1","GENE2"); "TEN1;TEN1,TEN1-CDK3" -> "TEN1"
split_genes <- function(gene_str, drop_readthrough = TRUE) {
  genes <- unique(unlist(strsplit(gene_str, "[;,]")))
  genes <- genes[!is.na(genes) & genes != "" & genes != "."]
  if (drop_readthrough && length(genes) > 1) {
    # readthrough X-Y is dropped when X or Y is also listed (keeps HLA-A, NKX2-1, etc.)
    is_readthrough <- vapply(genes, function(g) {
      parts <- strsplit(g, "-", fixed = TRUE)[[1]]
      length(parts) == 2 && any(parts %in% setdiff(genes, g))
    }, logical(1))
    genes <- genes[!is_readthrough]
  }
  if (length(genes) == 0) genes <- gene_str
  genes
}

# One row per (variant, gene): a variant annotated to two genes counts toward both
explode_genes <- function(ann, gene_col, drop_readthrough = TRUE) {
  gene_strs  <- unique(ann[[gene_col]])
  gene_lists <- lapply(gene_strs, split_genes, drop_readthrough = drop_readthrough)
  lookup     <- data.table::data.table(.gene_str = rep(gene_strs, lengths(gene_lists)),
                                       .gene     = unlist(gene_lists))
  ann <- data.table::as.data.table(ann)
  ann[, .gene_str := get(gene_col)]
  ann <- merge(ann, lookup, by = ".gene_str", allow.cartesian = TRUE)
  ann[, (gene_col) := .gene][, c(".gene_str", ".gene") := NULL]
  ann[]
}
