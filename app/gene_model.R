# gene_model.R — map exon-count bins onto the real Ensembl (GRCh38) gene model,
# draw a genome-browser-style track with per-group usage coverage and schematic
# (annotation-based) sashimi arcs. Live REST, cached per gene.
# NB: do NOT library(jsonlite) — its validate() would mask shiny::validate().
# plotly/data.table are attached by global.R; we use jsonlite:: with ::.

.model_cache <- new.env(parent = emptyenv())

# Fetch transcript/exon structure for one gene symbol from Ensembl REST (GRCh38).
# Returns list(chr, strand(+/-), gene_start, gene_end, canonical=data.table(exon_num,start,end),
#              transcripts=list of data.table(start,end) ) or NULL on failure.
get_gene_model <- function(gene) {
  if (exists(gene, envir = .model_cache)) return(get(gene, envir = .model_cache))
  url <- sprintf("https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s?expand=1;content-type=application/json",
                 utils::URLencode(gene))
  d <- tryCatch(jsonlite::fromJSON(url, simplifyDataFrame = FALSE), error = function(e) NULL)
  if (is.null(d) || is.null(d$Transcript)) { assign(gene, NULL, envir = .model_cache); return(NULL) }
  strand <- if (identical(d$strand, -1L) || identical(d$strand, -1)) "-" else "+"
  ex_dt <- function(tx) {
    e <- data.table::rbindlist(lapply(tx$Exon, function(x) list(start = x$start, end = x$end)))
    if (identical(strand, "-")) e[order(-start)] else e[order(start)]
  }
  tx <- d$Transcript
  can_i <- which(vapply(tx, function(t) isTRUE(t$is_canonical == 1), logical(1)))
  if (!length(can_i)) can_i <- 1L
  can <- ex_dt(tx[[can_i[1]]])
  can[, exon_num := seq_len(.N)]
  res <- list(
    gene = gene, chr = as.character(d$seq_region_name), strand = strand,
    gene_start = d$start, gene_end = d$end,
    canonical_tx = tx[[can_i[1]]]$id,
    canonical = can,
    transcripts = lapply(tx, function(t) {
      e <- ex_dt(t); list(id = t$id, biotype = t$biotype, is_canonical = isTRUE(t$is_canonical==1), exons = e)
    })
  )
  assign(gene, res, envir = .model_cache)
  res
}

# Tag each count bin (one row per bin) with the canonical exon it best overlaps.
# `bins` must have feature_start, feature_stop. Returns same rows + exon_num, exon_map.
map_bins_to_exons <- function(bins, model) {
  bins <- data.table::copy(bins)
  if (is.null(model)) { bins[, c("exon_num","exon_map") := list(NA_integer_, "no annotation")]; return(bins) }
  can <- model$canonical
  ov <- function(s, e) {
    o <- pmax(0L, pmin(e, can$end) - pmax(s, can$start) + 1L)
    if (all(o <= 0)) return(c(NA_integer_, NA_real_))
    i <- which.max(o); c(can$exon_num[i], o[i] / (e - s + 1L))
  }
  m <- t(mapply(ov, bins$feature_start, bins$feature_stop))
  bins[, exon_num := as.integer(m[,1])]
  bins[, frac := m[,2]]
  # union of all transcripts' exons, to tell "alternative exon" from true "intronic"
  alt <- data.table::rbindlist(lapply(model$transcripts, `[[`, "exons"))
  hits_alt <- function(s, e) any(pmin(e, alt$end) - pmax(s, alt$start) >= 0L)
  bins[, in_any := mapply(hits_alt, feature_start, feature_stop)]
  bins[, exon_map := data.table::fifelse(!is.na(exon_num),
         data.table::fifelse(frac >= 0.95, sprintf("exon %d", exon_num),
                                            sprintf("part of exon %d", exon_num)),
         data.table::fifelse(in_any, "alternative-isoform exon", "intronic"))]
  bins[, c("frac","in_any") := NULL]
  bins[]
}

# Build the genome-browser figure: canonical exons + arcs on top, per-group usage
# coverage tracks below, all on a shared genomic x-axis.
plot_gene_model <- function(gene, long, measure, groups, model, ylab) {
  # one value per bin per group (mean of measure), keep genomic extent
  agg <- long[group %in% groups,
              .(start = feature_start[1], stop = feature_stop[1], feature_id = feature_id[1],
                exon_label = exon_label[1], val = mean(get(measure))),
              by = .(exon_idx, group)]
  bins1 <- unique(long[, .(exon_idx, feature_start, feature_stop, feature_id)])
  if (!is.null(model)) bins1 <- map_bins_to_exons(bins1, model) else bins1[, exon_map := ""]
  agg <- bins1[, .(exon_idx, exon_map)][agg, on = "exon_idx"]

  vmax <- max(agg$val, na.rm = TRUE); if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  band_h <- 0.9
  floors <- c(AD = 2.2, Old = 1.1, Young = 0.0)
  exon_floor <- 3.5; exon_h <- 0.55; intron_y <- exon_floor + exon_h/2

  p <- plot_ly()

  # ---- per-group usage coverage (hoverable bars at genomic position) ----
  for (g in intersect(c("AD","Old","Young"), groups)) {
    a <- agg[group == g]
    if (!nrow(a)) next
    p <- add_trace(p, type = "bar",
      x = (a$start + a$stop)/2, y = a$val/vmax*band_h, width = pmax(a$stop - a$start, 30),
      base = floors[[g]], name = g, legendgroup = g,
      marker = list(color = GROUP_COLORS[[g]], line = list(width = 0)), opacity = 0.75,
      customdata = a$exon_map,
      hovertemplate = paste0(g, " — %{customdata}<br>", ylab,
                             ": ", sprintf("%.3g", a$val), "<extra></extra>"))
  }

  # ---- canonical exon boxes + numbers ----
  if (!is.null(model)) {
    can <- model$canonical
    p <- add_trace(p, type = "bar", x = (can$start + can$end)/2, y = rep(exon_h, nrow(can)),
      width = pmax(can$end - can$start, 40), base = exon_floor, name = "exon",
      marker = list(color = "#444"), showlegend = FALSE,
      customdata = can$exon_num,
      hovertemplate = paste0("exon %{customdata}<extra></extra>"))
    # always-visible exon numbers above each canonical exon
    p <- add_trace(p, type = "scatter", mode = "text",
      x = (can$start + can$end)/2, y = rep(exon_floor + exon_h + 0.12, nrow(can)),
      text = can$exon_num, textfont = list(size = 9, color = "#333"),
      hoverinfo = "skip", showlegend = FALSE)
    # intron line spanning the gene
    p <- add_trace(p, type = "scatter", mode = "lines",
      x = c(model$gene_start, model$gene_end), y = c(intron_y, intron_y),
      line = list(color = "#888", width = 1), hoverinfo = "skip", showlegend = FALSE)
    # ---- schematic sashimi arcs over introns (annotation, not read-supported) ----
    ce <- can[order(start)]
    if (nrow(ce) >= 2) {
      ax <- ay <- numeric(0)
      for (i in 1:(nrow(ce)-1)) {
        x0 <- ce$end[i]; x1 <- ce$start[i+1]
        if (x1 <= x0) next
        t <- seq(0, 1, length.out = 24)
        ax <- c(ax, x0 + (x1-x0)*t, NA)
        ay <- c(ay, (exon_floor+exon_h) + 0.7*sin(pi*t), NA)
      }
      p <- add_trace(p, type = "scatter", mode = "lines", x = ax, y = ay,
        line = list(color = "#7a7a7a", width = 1), hoverinfo = "skip",
        name = "intron (annotation)", showlegend = TRUE)
    }
  }

  # band labels
  anns <- list()
  lab <- function(y, txt, col) list(x = model$gene_start %||% min(agg$start), y = y, text = txt,
    xref = "x", yref = "y", xanchor = "right", showarrow = FALSE,
    font = list(size = 11, color = col))
  for (g in intersect(c("AD","Old","Young"), groups))
    anns <- c(anns, list(lab(floors[[g]] + band_h/2, g, GROUP_COLORS[[g]])))
  if (!is.null(model)) anns <- c(anns, list(lab(intron_y, "gene model", "#444")))

  ttl <- if (!is.null(model))
    sprintf("%s — chr%s (%s), canonical %s, %d exons", gene, model$chr, model$strand,
            model$canonical_tx, nrow(model$canonical))
  else sprintf("%s — (Ensembl model unavailable; showing usage by genomic position)", gene)

  layout(p, barmode = "overlay",
    title = list(text = ttl, x = 0.02, font = list(size = 13)),
    xaxis = list(title = "genomic position (bp, GRCh38)", zeroline = FALSE),
    yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE,
                 range = c(-0.2, exon_floor + exon_h + 0.9)),
    annotations = anns, legend = list(orientation = "h"),
    margin = list(l = 70))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
