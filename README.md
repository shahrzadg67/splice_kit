# splice_kit — human brain exon-usage explorer

Interactive R/Shiny app for exploring **exon-level usage / splicing** of any gene across
**Young**, **Old (aged control)** and **AD (Alzheimer's)** human brain, built from a public
exon-bin count matrix. Built with the GRIN1/GRIN2A glutamate-receptor splicing question in mind,
but works for any gene.

## Data

- Source matrix: `samples_exons_counts.tab` (exon-bin counts × 30 samples). **Not committed** —
  it is 165 MB; regenerate the app data from it with `app/prep_data.R`.
- Samples: SRA study **SRP287843 / PRJNA670209** ("integrated multi-omics ... Alzheimer's
  disease"), lateral temporal lobe postmortem brain, single-end 75 bp.
- Groups (decoded from GEO titles): **Young** n=8, **Old** n=10, **AD** n=12.

## What it shows (and doesn't)

This is exon-bin **count** data — no junction reads / BAMs. So the app shows **exon usage**
(inclusion vs skipping), the correct view for "which exons change," as a sashimi-style coverage
track. It does **not** draw junction arcs (a true sashimi plot needs the raw alignments).

- **Relative exon usage** = each exon's fraction of the gene's reads in a sample → isolates
  splicing from overall expression (the default view).
- **CPM** / **Raw** toggles for expression-level views.
- Exploratory per-exon **Kruskal–Wallis** across groups flags differing exons (FDR<0.05) — a hint,
  not a substitute for a formal DEXSeq/limma-diffSplice model.

## Layout

```
app/
  prep_data.R   build app data from ../samples_exons_counts.tab (run once)
  global.R      load data, get_gene(), normalization + KW helpers
  ui.R          gene picker, group/normalization controls, 3 tabs
  server.R      exon-usage track, sample heatmap, single-exon detail, CSV download
  deploy.R      publish to shinyapps.io
  data/         exons.rds, sample_meta.csv, gene_index.csv (built by prep_data.R)
  README_for_M.md  one-page guide for non-technical users
```

## Run locally

```bash
module load R/4.5.2
Rscript app/prep_data.R          # build app/data/ from the .tab (~3 min)
Rscript -e 'shiny::runApp("app", launch.browser=FALSE, port=8731)'
```

## Deploy

See `app/deploy.R` — set `SHINYAPPS_NAME/TOKEN/SECRET` then run it.
