# Human brain exon-usage explorer — quick start

**Live app:** https://shahrzadghs.shinyapps.io/brain-exon-usage/

1. Open the link in any web browser (no install needed).
2. In the **Gene** box, type a gene name — it defaults to **GRIN2A**. Try **GRIN1** too.
3. Look at the three tabs:
   - **Exon-usage track** — each exon along the gene (5'→3'); the three colored lines are
     **Young**, **Old (aged)** and **AD (Alzheimer's)** brain. A high point = that exon is used a
     lot; a dip in one group = that exon is relatively skipped in that group. A `*` marks exons
     that differ across groups.
   - **Exon map (gene model)** — your data lined up against the *real* Ensembl gene:
     the canonical transcript's exons with intron arcs on top, and Young/Old/AD usage as
     coverage below, at true genomic position. The table tells you which exon each bin is
     (e.g. "exon 5", "part of exon 2", or "alternative-isoform exon"). The arcs show the
     splice structure from annotation — they are not read-counted junctions.
   - **Sample heatmap** — every one of the 30 brain samples, exon by exon.
   - **Single-exon detail** — pick one exon, compare the three groups directly.
4. Top-left **Y-axis** toggle:
   - *Relative exon usage* (default) = the splicing view (each exon as a share of the gene).
   - *CPM* = overall expression level. *Raw* = raw counts.
5. **Download data (CSV)** saves the numbers behind the current gene.

**Data:** 30 human lateral-temporal-lobe brain samples (public study SRP287843): 8 Young,
10 Old, 12 Alzheimer's. This shows **exon usage** (inclusion/skipping), which is what these data
support — it is not a junction "sashimi" plot (that needs the raw alignments).
