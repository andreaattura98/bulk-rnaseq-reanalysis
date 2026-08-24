# Bulk RNA-seq analysis — ASD vs control neural progenitors

Re-analysis of the bulk RNA-seq data from **Mariani et al. (2015), *Cell*** —
*"FOXG1-Dependent Dysregulation of GABA/Glutamate Neuron Differentiation in Autism
Spectrum Disorders"* — using the [recount3](https://rna.recount.bio/) re-quantification of
the SRA study **`SRP047194`**.

The project starts from iPSC-derived neural progenitors of **ASD probands** and their
**unaffected relatives (controls)**, profiled at three terminal-differentiation time points
(day 0, day 11, day 31), and reproduces the differential-expression analysis of the paper.

## Goal

Recover the transcriptional signature reported by the authors — most notably the
up-regulation of **FOXG1** and of the GABAergic-fate genes (DLX family, GAD1) — starting from
the publicly available recount3 counts, with a clean and reproducible edgeR pipeline.

## Analysis pipeline (`asd_rnaseq_reanalysis.Rmd`)

1. **Data** — gene-level counts of `SRP047194` loaded from recount3 (48 samples, 63 856 genes).
2. **Sample annotation** — group (ASD vs control) and differentiation day parsed from the
   sample metadata.
3. **Exploratory analysis** — clinical table of the patients (`mmc2.xlsx`), per-sample
   log-count distributions, RLE plot, PCA (coloured by group and day).
4. **Differential expression** — edgeR on the raw filtered counts:
   `filterByExpr` → TMM library-size normalization → `glmQLFit` on a `~0 + group.day`
   design, testing **ASD vs control at TD11 and TD31**.
5. **Downstream** — result tables with gene symbols, explicit check of the paper's genes,
   MD / volcano / p-value plots, heatmap + silhouette of the top genes, and KEGG
   over-representation analysis (ORA).

## Two key methodological points

The analysis hinges on two decisions that are essential to reproduce the paper:

- **The comparison must be ASD vs control at matched time points**, *not* a time course
  (day 0 vs day 11/31) within the ASD group. FOXG1 and the GABAergic genes are
  *ASD-vs-control* differences; a within-group temporal contrast cannot surface them. Using
  the correct contrast is what brings the paper's genes back into the results.
- **No extra normalization is applied to the counts fed to edgeR.** edgeR performs its own
  (TMM) library-size normalization internally; the raw filtered counts are given to it
  directly. An additional GC-content / between-sample normalization is unnecessary here (and
  incorrect if pre-normalized counts are passed to edgeR).

## Results

Running the pipeline (verified end-to-end in R 4.6.1 / edgeR 4.10):

- **FOXG1 is recovered and up-regulated in ASD**, logFC ≈ **+1.3** at both TD11 and TD31 —
  the **same direction** reported by the authors (who measured an 8.5- and 13-fold increase).
- The GABAergic-fate genes of the paper (DLX1/2/5/6, GAD1, EOMES, NRXN1) also appear
  up-regulated in ASD, with **DLX6** the strongest signal (logFC ≈ +3.6 at TD31,
  FDR ≈ 0.03–0.06).
- FOXG1 itself does **not** cross the FDR < 0.05 threshold on this re-quantification
  (FDR ≈ 0.5–0.65), and the number of FDR-significant DEGs is smaller than the
  1062 (TD11) / 2203 (TD31) reported in the paper.

These differences are expected and are stated openly in the report: recount3 re-aligns and
re-annotates the reads with a different pipeline (Gencode G026) than the original study
(Tophat + GencodeV7), the cohort is small with high inter-individual variability, and the
authors used a more elaborate family/network model. What is reproduced here is the
**biological signal and its direction** (FOXG1 and the GABAergic programme up in ASD), rather
than the exact DEG count.

## How to reproduce

Open `asd_rnaseq_reanalysis.Rmd` in RStudio and click **Knit**, or from R:

```r
rmarkdown::render("asd_rnaseq_reanalysis.Rmd")
```

Requirements: R (≥ 4.4) with Bioconductor packages `recount3`, `edgeR`, `limma`, `EDASeq`,
`SummarizedExperiment`, `org.Hs.eg.db`, `clusterProfiler`, `enrichplot`, plus CRAN
`ggplot2`, `ggfortify`, `dplyr`, `corrplot`, `readxl`, `pheatmap`, `cluster`. The first-run
install chunks in the notebook are set to `eval=FALSE`. An internet connection is needed
(recount3 downloads the counts; KEGG ORA queries the KEGG API).

## Files

| File | Description |
|------|-------------|
| `asd_rnaseq_reanalysis.Rmd` | Analysis notebook (source) |
| `asd_rnaseq_reanalysis.html` | Rendered report (self-contained HTML) |
| `data/mmc2/mmc2.xlsx` | Supplementary clinical table from the paper |
| `DEG.RData` | Saved edgeR results — **not committed** (large; `.gitignore`d, regenerated on knit) |

## Credits

This analysis originated as a group project (Group 6) at the University of Padova.
Repository prepared and cleaned for publication by **Andrea Attura**.
Original contributors: _add names here_.

## Reference

Mariani J. *et al.* (2015). *FOXG1-Dependent Dysregulation of GABA/Glutamate Neuron
Differentiation in Autism Spectrum Disorders.* **Cell** 162(2):375-390.
[doi:10.1016/j.cell.2015.06.034](https://doi.org/10.1016/j.cell.2015.06.034) ·
[PMC4519016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4519016/) ·
Data: recount3 study `SRP047194`.
