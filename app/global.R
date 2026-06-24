# global.R — loaded once at app startup
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(data.table)
  library(plotly)
})

source("gene_model.R")  # Ensembl exon mapping + gene-model/sashimi track

DATA_DIR <- file.path("data")

exons       <- readRDS(file.path(DATA_DIR, "exons.rds"))      # keyed data.table on gene_name
sample_meta <- fread(file.path(DATA_DIR, "sample_meta.csv"))
gene_index  <- fread(file.path(DATA_DIR, "gene_index.csv"))
setkey(exons, gene_name)

SAMPLE_COLS <- sample_meta$SRR
GROUP_LEVELS <- c("Young", "Old", "AD")
GROUP_COLORS <- c(Young = "#2c7fb8", Old = "#fec44f", AD = "#d7301f")
sample_meta[, group := factor(group, levels = GROUP_LEVELS)]
setkey(sample_meta, SRR)

# Genes offered in the dropdown: those with any expression (keeps the list useful).
GENE_CHOICES <- gene_index[total_count > 0][order(-total_count), gene_name]

NORM_CHOICES <- c(
  "Relative exon usage (splicing)" = "rel",
  "CPM (expression-normalized)"    = "cpm",
  "Raw counts"                     = "raw"
)
NORM_YLAB <- c(rel = "Fraction of gene's reads (PSI proxy)",
               cpm = "Counts per million",
               raw = "Raw read count")

# ---- per-gene long table: one row per (exon bin x sample), 5'->3' ordered ----
get_gene <- function(gene) {
  g <- exons[gene]
  if (nrow(g) == 0L || is.na(g$gene_name[1])) return(NULL)
  strand <- g$strand[1]
  ord <- if (identical(strand, "-")) order(-g$feature_start) else order(g$feature_start)
  g <- g[ord]
  g[, exon_idx := seq_len(.N)]
  g[, exon_label := sprintf("E%02d", exon_idx)]

  long <- melt(g,
               id.vars = c("exon_idx","exon_label","feature_id","chr","strand",
                           "feature_start","feature_stop","length"),
               measure.vars = SAMPLE_COLS, variable.name = "SRR", value.name = "count")
  long[, SRR := as.character(SRR)]
  long <- sample_meta[long, on = "SRR"]           # adds group, lib_size, GSM

  # gene total per sample (denominator for relative usage)
  long[, gene_total := sum(count), by = SRR]
  long[, rel := fifelse(gene_total > 0, count / gene_total, 0)]
  long[, cpm := count / lib_size * 1e6]
  long[, raw := as.numeric(count)]
  long[, exon_label := factor(exon_label, levels = unique(exon_label[order(exon_idx)]))]
  long[, group := factor(group, levels = GROUP_LEVELS)]
  long[]
}

# Per-exon Kruskal-Wallis across groups on a given measure (exploratory DEU hint)
exon_kw <- function(long, measure = "rel") {
  long <- long[!is.na(group)]
  res <- long[, {
    v <- get(measure); ok <- is.finite(v)
    p <- tryCatch(
      if (length(unique(group[ok])) >= 2 && sum(ok) >= 3)
        kruskal.test(v[ok], group[ok])$p.value else NA_real_,
      error = function(e) NA_real_)
    .(p = p)
  }, by = .(exon_idx, exon_label)]
  res[, fdr := p.adjust(p, "BH")]
  res[]
}
