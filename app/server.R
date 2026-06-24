# server.R
function(input, output, session) {

  updateSelectizeInput(session, "gene", choices = GENE_CHOICES,
                       selected = if ("GRIN2A" %in% GENE_CHOICES) "GRIN2A" else GENE_CHOICES[1],
                       server = TRUE)

  geneData <- reactive({
    req(input$gene)
    d <- get_gene(input$gene)
    validate(need(!is.null(d), "Gene not found."))
    d
  })

  # groups currently selected
  selData <- reactive({
    d <- geneData()
    d[group %in% input$groups]
  })

  measure <- reactive(input$norm)
  ylab    <- reactive(NORM_YLAB[[measure()]])

  # per-exon significance across ALL groups (independent of the display filter)
  sigTab <- reactive({
    if (!isTRUE(input$show_sig)) return(NULL)
    exon_kw(geneData(), measure())
  })

  output$gene_info <- renderUI({
    d <- geneData()
    gi <- gene_index[gene_name == input$gene]
    tagList(
      strong(input$gene), br(),
      sprintf("%s | chr%s (%s)", gi$gene_id[1], d$chr[1], d$strand[1]), br(),
      sprintf("%d exon bins", uniqueN(d$exon_idx))
    )
  })

  # ---- Tab 1: exon-usage track ----
  output$track <- renderPlotly({
    d <- selData(); m <- measure()
    req(nrow(d) > 0)
    agg <- d[, .(mean = mean(get(m)), sd = sd(get(m)), n = .N),
             by = .(exon_idx, exon_label, group)]
    setorder(agg, exon_idx)
    xlevels <- levels(d$exon_label)

    p <- plot_ly()
    for (g in intersect(GROUP_LEVELS, input$groups)) {
      a <- agg[group == g]
      p <- add_trace(p, data = a, x = ~exon_label, y = ~mean, type = "scatter",
                     mode = "lines+markers", name = g,
                     line = list(color = GROUP_COLORS[[g]]),
                     marker = list(color = GROUP_COLORS[[g]]),
                     error_y = list(array = a$sd, color = GROUP_COLORS[[g]], thickness = 1),
                     hovertemplate = paste0(g, " — %{x}<br>", ylab(), ": %{y:.3g}<extra></extra>"))
    }
    # significance flags
    st <- sigTab()
    if (!is.null(st)) {
      sig <- st[is.finite(fdr) & fdr < 0.05]
      if (nrow(sig) > 0) {
        ymax <- max(agg$mean + ifelse(is.na(agg$sd), 0, agg$sd), na.rm = TRUE)
        p <- add_trace(p, x = sig$exon_label, y = rep(ymax * 1.06, nrow(sig)),
                       type = "scatter", mode = "text", text = "*",
                       textfont = list(size = 18, color = "black"),
                       name = "FDR<0.05", hoverinfo = "skip", showlegend = TRUE)
      }
    }
    layout(p,
           title = list(text = paste0(input$gene, " — exon usage (5'→3')"), x = 0.02),
           xaxis = list(title = "exon bin", categoryorder = "array", categoryarray = xlevels,
                        tickangle = -45),
           yaxis = list(title = ylab()),
           legend = list(orientation = "h"))
  })

  # ---- Tab 2: heatmap (exon x sample) ----
  output$heatmap <- renderPlotly({
    d <- selData(); m <- measure()
    req(nrow(d) > 0)
    smeta <- sample_meta[group %in% input$groups]
    setorder(smeta, group, SRR)
    col_order <- smeta$SRR
    d[, exon_label := factor(exon_label, levels = rev(levels(d$exon_label)))]
    mat <- dcast(d, exon_label ~ SRR, value.var = m)
    rn  <- as.character(mat$exon_label); mat[, exon_label := NULL]
    M <- as.matrix(mat)[, col_order, drop = FALSE]; rownames(M) <- rn
    col_grp <- as.character(smeta$group)
    plot_ly(x = paste0(col_order, " (", col_grp, ")"), y = rn, z = M,
            type = "heatmap", colors = "YlOrRd",
            colorbar = list(title = ylab()),
            hovertemplate = "exon %{y}<br>%{x}<br>%{z:.3g}<extra></extra>") |>
      layout(title = list(text = paste0(input$gene, " — exon x sample"), x = 0.02),
             xaxis = list(title = "sample (group)", tickangle = -60),
             yaxis = list(title = "exon bin"))
  })

  # ---- Tab 3: single-exon detail ----
  observeEvent(geneData(), {
    d <- geneData()
    lv <- levels(d$exon_label)
    updateSelectInput(session, "exon", choices = lv, selected = lv[1])
  })

  output$exon_coords <- renderText({
    req(input$exon)
    d <- geneData()[exon_label == input$exon]
    req(nrow(d) > 0)
    sprintf("%s  |  chr%s:%s-%s (%s)  |  %s bp  |  feature_id %s",
            input$exon, d$chr[1], d$feature_start[1], d$feature_stop[1], d$strand[1],
            d$length[1], d$feature_id[1])
  })

  output$exon_box <- renderPlotly({
    req(input$exon)
    d <- selData()[exon_label == input$exon]; m <- measure()
    req(nrow(d) > 0)
    plot_ly(d, x = ~group, y = ~get(m), color = ~group, colors = GROUP_COLORS,
            type = "box", boxpoints = "all", jitter = 0.4, pointpos = 0,
            text = ~SRR, hovertemplate = "%{text}<br>%{y:.3g}<extra></extra>") |>
      layout(title = list(text = paste0(input$gene, " — ", input$exon), x = 0.02),
             xaxis = list(title = ""), yaxis = list(title = ylab()), showlegend = FALSE)
  })

  # ---- download ----
  output$dl_csv <- downloadHandler(
    filename = function() paste0(input$gene, "_exon_usage.csv"),
    content = function(file) {
      d <- selData()[, .(gene = input$gene, exon = exon_label, feature_id, chr, strand,
                         feature_start, feature_stop, length, SRR, group, count, rel, cpm)]
      fwrite(d, file)
    }
  )
}
