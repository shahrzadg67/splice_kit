# ui.R
fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  titlePanel("Human brain exon-usage explorer — aging & Alzheimer's"),
  div(style = "color:#555; margin-bottom:10px;",
      "Exon-level read usage across ", strong("Young"), ", ", strong("Old (aged control)"),
      " and ", strong("AD (Alzheimer's)"), " lateral temporal lobe (SRA SRP287843, n=30). ",
      "Pick a gene to see which exons are included/skipped per group. ",
      em("Note: based on exon-bin counts — shows exon usage, not junction (sashimi) arcs.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectizeInput("gene", "Gene", choices = NULL,
                     options = list(placeholder = "type a gene, e.g. GRIN2A",
                                    maxOptions = 200)),
      radioButtons("norm", "Y-axis / normalization", choices = NORM_CHOICES, selected = "rel"),
      checkboxGroupInput("groups", "Groups", choices = GROUP_LEVELS, selected = GROUP_LEVELS),
      checkboxInput("show_sig", "Flag exons differing across groups (Kruskal–Wallis FDR<0.05)", TRUE),
      hr(),
      uiOutput("gene_info"),
      hr(),
      downloadButton("dl_csv", "Download data (CSV)", class = "btn-sm"),
      br(), br(),
      div(style = "font-size:11px;color:#888;",
          "Relative usage = each exon's fraction of the gene's reads in a sample ",
          "(removes overall expression differences, isolating splicing).")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("Exon-usage track",
                 br(),
                 plotlyOutput("track", height = "460px")),
        tabPanel("Exon map (gene model)",
                 br(),
                 div(style = "color:#666;font-size:12px;margin-bottom:6px;",
                     "Count bins mapped onto the real Ensembl GRCh38 gene model. Top: canonical ",
                     "transcript exons with intron arcs (annotation skeleton, ",
                     em("not"), " read-supported junctions). Below: per-group usage as coverage ",
                     "at true genomic position. Hover a bar for its exon."),
                 plotlyOutput("model", height = "520px"),
                 br(),
                 h5("Bin → exon mapping"),
                 tableOutput("map_tbl")),
        tabPanel("Sample heatmap",
                 br(),
                 plotlyOutput("heatmap", height = "560px")),
        tabPanel("Single-exon detail",
                 br(),
                 fluidRow(
                   column(4, selectInput("exon", "Exon bin", choices = NULL)),
                   column(8, div(style="margin-top:28px;color:#666;", textOutput("exon_coords")))),
                 plotlyOutput("exon_box", height = "420px"))
      )
    )
  )
)
