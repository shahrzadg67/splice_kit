#!/usr/bin/env Rscript
# deploy.R — publish the app to shinyapps.io
#
# One-time setup (get these from https://www.shinyapps.io/admin/#/tokens):
#   export SHINYAPPS_NAME=<account>
#   export SHINYAPPS_TOKEN=<token>
#   export SHINYAPPS_SECRET=<secret>
# Then:  module load R/4.5.2 && Rscript splice/app/deploy.R
#
# rsconnect is not in the R module; install once into your personal library if needed.

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  message("Installing rsconnect into personal library ...")
  install.packages("rsconnect", repos = "https://cloud.r-project.org")
}
library(rsconnect)

acc <- Sys.getenv("SHINYAPPS_NAME"); tok <- Sys.getenv("SHINYAPPS_TOKEN"); sec <- Sys.getenv("SHINYAPPS_SECRET")
if (!nzchar(acc) || !nzchar(tok) || !nzchar(sec))
  stop("Set SHINYAPPS_NAME / SHINYAPPS_TOKEN / SHINYAPPS_SECRET env vars first.")

setAccountInfo(name = acc, token = tok, secret = sec)

app_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(app_dir) || !nzchar(app_dir)) app_dir <- "."

deployApp(
  appDir   = app_dir,
  appName  = "brain-exon-usage",
  appFiles = c("global.R", "ui.R", "server.R", "gene_model.R",
               "data/exons.rds", "data/sample_meta.csv", "data/gene_index.csv"),
  forceUpdate = TRUE
)
