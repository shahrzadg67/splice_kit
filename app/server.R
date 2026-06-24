# server.R
function(input, output, session) {

  updateSelectizeInput(session, "gene", choices = GENE_CHOICES,
                       selected = if ("GRIN2A" %in% GENE_CHOICES) "GRIN2A" else GENE_CHOICES[1],
                       server = TRUE)

  geneData <- reactive({
    req(input$gene)
    d <- get_gene(input$gene)
    shiny::validate(shiny::need(!is.null(d), "Gene not found."))
    d
  })

  # Ensembl gene model (cached, live REST) — used by every tab for exon labels
  geneModel <- reactive({
    req(input$gene)
    withProgress(message = "Fetching Ensembl gene model...", value = 0.5,
                 get_gene_model(input$gene))
  })

  # bin -> exon mapping for the current gene (exon_num, exon_map text, short display tag)
  mapInfo <- reactive({
    d <- geneData(); m <- geneModel()
    bins <- unique(d[, .(exon_idx, exon_label, feature_start, feature_stop)])
    mp <- if (!is.null(m)) map_bins_to_exons(bins, m)
          else { b <- copy(bins); b[, c("exon_num","exon_map") := list(NA_integer_, "")]; b }
    mp[, tag := fifelse(!is.na(exon_num), paste0("ex", exon_num),
                 fifelse(exon_map == "alternative-isoform exon", "alt",
                 fifelse(exon_map == "intronic", "int", "")))]
    mp[, exon_disp := fifelse(nzchar(tag), paste0(as.character(exon_label), " ·", tag),
                                            as.character(exon_label))]
    setorder(mp, exon_idx)
    mp[, .(exon_idx, exon_num, exon_map, exon_disp)]
  })

  # geneData enriched with mapping columns + ordered exon_disp factor
  mappedData <- reactive({
    mi <- mapInfo()
    d <- mi[geneData(), on = "exon_idx"]
    lev <- mi$exon_disp[order(mi$exon_idx)]
    d[, exon_disp := factor(exon_disp, levels = lev)]
    d[]
  })

  selData <- reactive(mappedData()[group %in% input$groups])

  measure <- reactive(input$norm)
  ylab    <- reactive(NORM_YLAB[[measure()]])

  sigTab <- reactive({
    if (!isTRUE(input$show_sig)) return(NULL)
    exon_kw(geneData(), measure())
  })

  output$gene_info <- renderUI({
    d <- geneData(); m <- geneModel()
    gi <- gene_index[gene_name == input$gene]
    tagList(
      strong(input$gene), br(),
      sprintf("%s | chr%s (%s)", gi$gene_id[1], d$chr[1], d$strand[1]), br(),
      sprintf("%d count bins", uniqueN(d$exon_idx)),
      if (!is.null(m)) tagList(br(), sprintf("canonical %s: %d exons", m$canonical_tx, nrow(m$canonical)))
    )
  })

  # ---- Tab 1: exon-usage track (x-axis & hover carry exon numbers) ----
  output$track <- renderPlotly({
    d <- selData(); m <- measure()
    req(nrow(d) > 0)
    agg <- d[, .(mean = mean(get(m)), sd = sd(get(m)),
                 exon_map = exon_map[1], exon_disp = exon_disp[1]),
             by = .(exon_idx, group)]
    setorder(agg, exon_idx)
    xlevels <- levels(d$exon_disp)

    p <- plot_ly()
    for (g in intersect(GROUP_LEVELS, input$groups)) {
      a <- agg[group == g]
      p <- add_trace(p, data = a, x = ~exon_disp, y = ~mean, type = "scatter",
                     mode = "lines+markers", name = g,
                     line = list(color = GROUP_COLORS[[g]]),
                     marker = list(color = GROUP_COLORS[[g]]),
                     error_y = list(array = a$sd, color = GROUP_COLORS[[g]], thickness = 1),
                     customdata = a$exon_map,
                     hovertemplate = paste0(g, " — %{x}<br>maps to: %{customdata}<br>",
                                            ylab(), ": %{y:.3g}<extra></extra>"))
    }
    st <- sigTab()
    if (!is.null(st)) {
      sig <- merge(st[is.finite(fdr) & fdr < 0.05, .(exon_idx)],
                   unique(d[, .(exon_idx, exon_disp)]), by = "exon_idx")
      if (nrow(sig) > 0) {
        ymax <- max(agg$mean + ifelse(is.na(agg$sd), 0, agg$sd), na.rm = TRUE)
        p <- add_trace(p, x = sig$exon_disp, y = rep(ymax * 1.06, nrow(sig)),
                       type = "scatter", mode = "text", text = "*",
                       textfont = list(size = 18, color = "black"),
                       name = "FDR<0.05", hoverinfo = "skip", showlegend = TRUE)
      }
    }
    layout(p,
           title = list(text = paste0(input$gene, " — exon usage (5'→3'); ·exN = Ensembl exon"), x = 0.02),
           xaxis = list(title = "count bin · exon", categoryorder = "array",
                        categoryarray = xlevels, tickangle = -45),
           yaxis = list(title = ylab()),
           legend = list(orientation = "h"))
  })

  # ---- Tab: exon map / gene model ----
  output$model <- renderPlotly({
    d <- selData(); req(nrow(d) > 0)
    plot_gene_model(input$gene, d, measure(), input$groups, geneModel(), ylab())
  })

  output$map_tbl <- renderTable({
    d <- geneData(); m <- geneModel()
    bins <- unique(d[, .(exon_idx, exon_label, chr, feature_start, feature_stop, length, feature_id)])
    bins <- if (!is.null(m)) map_bins_to_exons(bins, m) else { bins[, exon_map := "n/a"]; bins }
    setorder(bins, exon_idx)
    bins[, .(`bin` = exon_label, `chr` = chr, `start` = feature_start, `stop` = feature_stop,
             `bp` = length, `maps to` = exon_map)]
  }, striped = TRUE, spacing = "xs", width = "100%")

  # ---- Tab 2: heatmap (rows labeled with exon number) ----
  output$heatmap <- renderPlotly({
    d <- selData(); m <- measure()
    req(nrow(d) > 0)
    smeta <- sample_meta[group %in% input$groups]
    setorder(smeta, group, SRR)
    col_order <- smeta$SRR
    d[, exon_disp := factor(exon_disp, levels = rev(levels(d$exon_disp)))]
    mat <- dcast(d, exon_disp ~ SRR, value.var = m)
    rn  <- as.character(mat$exon_disp); mat[, exon_disp := NULL]
    M <- as.matrix(mat)[, col_order, drop = FALSE]; rownames(M) <- rn
    col_grp <- as.character(smeta$group)
    plot_ly(x = paste0(col_order, " (", col_grp, ")"), y = rn, z = M,
            type = "heatmap", colors = "YlOrRd",
            colorbar = list(title = ylab()),
            hovertemplate = "%{y}<br>%{x}<br>%{z:.3g}<extra></extra>") |>
      layout(title = list(text = paste0(input$gene, " — bin (·exon) x sample"), x = 0.02),
             xaxis = list(title = "sample (group)", tickangle = -60),
             yaxis = list(title = "count bin · exon"))
  })

  # ---- Tab 3: single-exon detail (dropdown shows exon number) ----
  observeEvent(mappedData(), {
    mi <- mapInfo(); setorder(mi, exon_idx)
    d <- geneData()
    # named choices: label shows "E03 ·ex2", value is the bin label "E03"
    ch <- setNames(as.character(unique(d[order(exon_idx)]$exon_label)), mi$exon_disp)
    updateSelectInput(session, "exon", choices = ch, selected = ch[1])
  })

  output$exon_coords <- renderText({
    req(input$exon)
    d <- mappedData()[exon_label == input$exon]
    req(nrow(d) > 0)
    sprintf("%s  →  %s  |  chr%s:%s-%s (%s)  |  %s bp  |  feature_id %s",
            input$exon, d$exon_map[1], d$chr[1], d$feature_start[1], d$feature_stop[1],
            d$strand[1], d$length[1], d$feature_id[1])
  })

  output$exon_box <- renderPlotly({
    req(input$exon)
    d <- selData()[exon_label == input$exon]; m <- measure()
    req(nrow(d) > 0)
    ttl <- sprintf("%s — %s (%s)", input$gene, input$exon, d$exon_map[1])
    plot_ly(d, x = ~group, y = ~get(m), color = ~group, colors = GROUP_COLORS,
            type = "box", boxpoints = "all", jitter = 0.4, pointpos = 0,
            text = ~SRR, hovertemplate = "%{text}<br>%{y:.3g}<extra></extra>") |>
      layout(title = list(text = ttl, x = 0.02),
             xaxis = list(title = ""), yaxis = list(title = ylab()), showlegend = FALSE)
  })

  # ---- download (includes exon mapping) ----
  output$dl_csv <- downloadHandler(
    filename = function() paste0(input$gene, "_exon_usage.csv"),
    content = function(file) {
      d <- selData()
      out <- d[, .(gene = input$gene, bin = exon_label, maps_to = exon_map, feature_id, chr, strand,
                   feature_start, feature_stop, length, SRR, group, count, rel, cpm)]
      fwrite(out, file)
    }
  )
}
