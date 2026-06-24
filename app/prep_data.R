#!/usr/bin/env Rscript
# prep_data.R — build compact app data from samples_exons_counts.tab
# Run once on HPC:  module load R/4.5.2 && Rscript splice/app/prep_data.R
#
# Inputs : ../samples_exons_counts.tab  (exon-bin count matrix, 30 SRA samples)
# Outputs: data/exons.rds        keyed data.table (gene_name) of exon bins + 30 counts
#          data/sample_meta.csv  SRR, GSM, group (Young/Old/AD), lib_size
#          data/gene_index.csv   gene_name, gene_id, n_exons, total_count (for dropdown)

suppressPackageStartupMessages({
  library(data.table)
})

app_dir  <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(app_dir) || !nzchar(app_dir)) app_dir <- "."
in_tab   <- normalizePath(file.path(app_dir, "..", "samples_exons_counts.tab"), mustWork = TRUE)
data_dir <- file.path(app_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading ", in_tab, " ...")
dt <- fread(in_tab, sep = "\t", header = TRUE)

meta_cols   <- c("gene_id","gene_name","chr","strand","feature_start","feature_stop","length","feature_id")
sample_cols <- setdiff(names(dt), meta_cols)
stopifnot(all(meta_cols %in% names(dt)), length(sample_cols) == 30L)

# ---- sample -> group map (decoded from SRA SRP287843 / GEO titles) ----
old_ids   <- paste0("SRR128508", 30:39)
ad_ids    <- paste0("SRR128508", c(40,41,43,44,45,46,47,48,49,50,51,52))
young_ids <- paste0("SRR128508", c(42,53,54,55,56,57,58,59))
group_of <- function(x) fifelse(x %in% old_ids, "Old",
                         fifelse(x %in% ad_ids, "AD",
                         fifelse(x %in% young_ids, "Young", NA_character_)))
stopifnot(setequal(sample_cols, c(old_ids, ad_ids, young_ids)))

# GSM accessions are sequential with the SRR order (SRR12850830 -> GSM4837814)
srr_sorted <- sort(sample_cols)
gsm_map    <- setNames(paste0("GSM", 4837814 + seq_along(srr_sorted) - 1L), srr_sorted)

lib_size <- vapply(sample_cols, function(s) sum(dt[[s]], na.rm = TRUE), numeric(1))
sample_meta <- data.table(
  SRR      = sample_cols,
  GSM      = gsm_map[sample_cols],
  group    = group_of(sample_cols),
  lib_size = as.numeric(lib_size[sample_cols])
)
setorder(sample_meta, group, SRR)
stopifnot(!anyNA(sample_meta$group), all(sample_meta$lib_size > 0))
fwrite(sample_meta, file.path(data_dir, "sample_meta.csv"))
message("Groups: ", paste(sample_meta[, .N, by = group][, paste0(group, "=", N)], collapse = ", "))

# ---- gene index (for the searchable dropdown) ----
dt[, total_count := rowSums(.SD), .SDcols = sample_cols]
gene_index <- dt[, .(gene_id = gene_id[1], n_exons = .N, total_count = sum(total_count)),
                 by = gene_name]
setorder(gene_index, -total_count)
fwrite(gene_index, file.path(data_dir, "gene_index.csv"))
dt[, total_count := NULL]

# ---- compact keyed exon table ----
keep <- c("gene_name","gene_id","chr","strand","feature_start","feature_stop","length","feature_id", sample_cols)
exons <- dt[, ..keep]
for (s in sample_cols) set(exons, j = s, value = as.integer(exons[[s]]))
setkey(exons, gene_name)
saveRDS(exons, file.path(data_dir, "exons.rds"), compress = "xz")

# ---- sanity asserts ----
n_grin2a <- exons["GRIN2A", .N]
n_grin1  <- exons["GRIN1",  .N]
message(sprintf("GRIN2A bins=%d  GRIN1 bins=%d  genes=%d  rows=%d",
                n_grin2a, n_grin1, nrow(gene_index), nrow(exons)))
stopifnot(n_grin2a == 52L, n_grin1 == 44L)

sz <- file.info(file.path(data_dir, "exons.rds"))$size
message(sprintf("exons.rds = %.1f MB", sz/1e6))
message("prep_data.R done.")
