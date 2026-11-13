# suchilactone-ctsb-ad-as-code

Analysis code for the study **"Hawthorn-derived suchilactone engages cathepsin B to
attenuate NF-κB/NLRP3-driven cerebrovascular endothelial injury"**.

This repository contains the 41 R scripts behind every computational result in the
manuscript: a transcriptomic meta-analysis of eight public GEO cohorts, network
pharmacology of hawthorn (*Crataegus pinnatifida*), nine-algorithm diagnostic
modelling with SHAP interpretation, single-cell transcriptomics of Alzheimer's
disease and atherosclerotic tissue, in silico CTSB perturbation, and the
molecular-docking preparation.

All scripts sit in one place, named with a two-digit run-order prefix. Run them
top to bottom; each states its input and output in the header comment.

---

## Quick start

```r
# 1. install the R packages (once; R >= 4.4.3)
source("install_packages.R")

# 2. download the public data from GEO
source("00_download_data.R")

# 3. run the scripts in numerical order
source("01_bulk_AD_GSE5281_normalize.R")
...
```

`00_config.R` holds the shared settings (paths, seeds, DEG thresholds) and **must
stay in the repository root**: every script begins with a bootstrap block that
walks up from the working directory until it finds that file, then sets
`PROJECT_ROOT` and creates `work/` for intermediate files.

Notation: `AD` = Alzheimer's disease, `AS` = atherosclerosis, `sc` = single-cell,
`knk` = scTenifoldKnk, `ML` = machine learning.

---

## Phase 1 — bulk transcriptomic meta-analysis (manuscript 2.1–2.2)

One cohort per script; the eight cohorts are normalised independently, then
concatenated and batch-corrected.

| # | File | Step | Input → output |
| --- | --- | --- | --- |
| 01 | `01_bulk_AD_GSE5281_normalize.R` | 1.1 | raw series matrix → `GSE5281.normalize.txt` |
| 02 | `02_bulk_AD_GSE138260_normalize.R` | 1.1 | → `GSE138260.normalize.txt` |
| 03 | `03_bulk_AD_GSE29378_normalize.R` | 1.1 | → `GSE29378.normalize.txt` |
| 04 | `04_bulk_AD_GSE36980_normalize.R` | 1.1 | → `GSE36980.normalize.txt` |
| 05 | `05_bulk_AD_GSE37263_normalize.R` | 1.1 | → `GSE37263.normalize.txt` |
| 06 | `06_bulk_AS_GSE100927_normalize.R` | 1.1 | → `GSE100927.normalize.txt` |
| 07 | `07_bulk_AS_GSE100927_merge_platforms.R` | 1.1b | sub-series, probe → gene symbol |
| 08 | `08_bulk_AS_GSE28829_normalize.R` | 1.1 | → `GSE28829.normalize.txt` |
| 09 | `09_bulk_AS_GSE43292_normalize.R` | 1.1 | → `GSE43292.normalize.txt` |
| 10 | `10_bulk_AD_combat_batch_correction.R` | 1.2 | five AD cohorts → `merge.preNorm.txt`, `merge.normalize.txt` |
| 11 | `11_bulk_AS_combat_batch_correction.R` | 1.2 | three AS cohorts → same two files |
| 12 | `12_bulk_AD_PCA_before_after.R` | 1.3 | `merge.preNorm/normalize.txt` → PCA panels |
| 13 | `13_bulk_AS_PCA_before_after.R` | 1.3 | same for AS |
| 14 | `14_bulk_AD_limma_DEG.R` | 1.4 | `merge.normalize.txt` → `all.txt`, `diff.xls`, volcano/heatmap |
| 15 | `15_bulk_AS_limma_DEG.R` | 1.4 | same for AS |
| 16 | `16_bulk_venn_shared_DEGs.R` | 1.5 | both `diff.xls` → `interGenes.List.txt`, Venn diagrams |

## Phase 2 — network pharmacology (manuscript 2.3–2.4)

| # | File | Step | Input → output |
| --- | --- | --- | --- |
| 17 | `17_netpharm_swissADME_filter.R` | 2.1 | 26 SwissADME exports → `swissadme.filter.csv` |
| 18 | `18_netpharm_merge_compound_targets.R` | 2.2 | TCMSP + SwissTargetPrediction + SEA → merged target table |
| 19 | `19_netpharm_drug_disease_intersection.R` | 2.3 | herb targets ∩ `interGenes.txt` → shared targets |
| 20 | `20_netpharm_deduplicate_drug_targets.R` | 2.4 | removes duplicated compound–target edges |

The PPI network itself is built in the Cytoscape GUI (STRING, *Homo sapiens*,
combined score > 0.4); the exported session and edge list are part of the full
reproducibility package.

## Phase 3 — functional enrichment and GSEA (manuscript 2.4, 2.8)

| # | File | Step | Input → output |
| --- | --- | --- | --- |
| 21 | `21_bulk_GO_enrichment.R` | 1.6 | `interGenes.txt` → GO bubble/bar plots |
| 22 | `22_bulk_KEGG_enrichment.R` | 1.6 | `interGenes.txt` → KEGG bubble/bar plots |
| 23 | `23_bulk_AD_GSEA_core_genes.R` | 1.7 | `merge.normalize.txt` + KEGG gmt → GSEA for CTSB / P2RX7 |
| 24 | `24_bulk_AS_GSEA_core_genes.R` | 1.7 | same for AS |

## Phase 4 — diagnostic modelling (manuscript 2.5)

| # | File | Step | Input → output |
| --- | --- | --- | --- |
| 25 | `25_bulk_ML_prepare_training_data.R` | 1.8 | expression matrices + `interGenes.txt` → `data.train.txt`, `data.test.txt` |
| 26 | `26_bulk_ML_train_nine_algorithms.R` | 1.9 | nine classifiers (`caret`) → `importanceGene.<algorithm>.txt`, ROC/box panels |
| 27 | `27_bulk_ML_SHAP_importance.R` | 1.10 | best model → SHAP beeswarm and ranking |
| 28 | `28_bulk_ML_ROC_curves.R` | 1.11 | `merge.normalize.txt` → diagnostic ROC curves |
| 29 | `29_bulk_ML_nomogram_decision_curve.R` | 1.12 | nomogram, calibration, decision-curve analysis |

## Phase 5 — single-cell transcriptomics (manuscript 2.6)

| # | File | Step | Input → output |
| --- | --- | --- | --- |
| 30 | `30_sc_QC_clustering_annotation.R` | 3.1 | GSE157827 / GSE159677 → QC, clustering, SingleR labels |
| 31 | `31_sc_manual_cell_annotation.R` | 3.2 | manual refinement of cell types |
| 32 | `32_sc_ctsb_p2rx7_expression.R` | 3.3 | feature / dot / violin plots for CTSB and P2RX7 |
| 33 | `33_sc_myeloid_subclustering.R` | 3.4 | macrophage (AS) and microglia (AD) subclusters |
| 34 | `34_sc_ctsb_high_vs_low_differential.R` | 3.5 | CTSB-high vs CTSB-low DEGs and enrichment |
| 35 | `35_sc_pathway_activity_correlation.R` | 3.6 | module scores and Spearman correlation with CTSB |
| 36 | `36_sc_cellchat_communication.R` | 3.7 | CellChat ligand–receptor networks |
| 37 | `37_sc_gene_expression_visualisation.R` | 3.8 | publication-quality expression panels |
| 38 | `38_sc_full_pipeline_reference.R` | 3.9 | consolidated reference copy of the whole single-cell run |

## Phase 6 — in silico perturbation (manuscript 2.8)

| # | File | Step | Input → output |
| --- | --- | --- | --- |
| 39 | `39_knk_sctenifoldknk_ctsb.R` | 4.1 | CTSB-centred network → `sigDiff` tables |
| 40 | `40_knk_enrichment_after_perturbation.R` | 4.2 | GO / KEGG / Reactome / DisGeNET / GSEA of the perturbed network |

## Step 0 — data download

| # | File | Purpose |
| --- | --- | --- |
| 00 | `00_download_data.R` | fetches the eight bulk cohorts and the two single-cell series from GEO |

---

## Data

Nothing large is stored in this repository. The expression data are public:

| Data | Accession |
| --- | --- |
| AD bulk cohorts | GSE138260, GSE29378, GSE36980, GSE37263, GSE5281 (356 samples) |
| AS bulk cohorts | GSE100927, GSE28829, GSE43292 (197 samples) |
| AD single-cell | GSE157827 (21 samples) |
| AS single-cell | GSE159677 (6 samples) |
| Compound–target and SwissADME tables, DEG tables, STRING edge list | available from the corresponding author on request |

`00_download_data.R` retrieves the series matrices; the sample grouping of every
cohort is in the `metadata/` folder of the full reproducibility package.

---

## Requirements

* **R ≥ 4.4.3**
* Bioconductor and CRAN packages: `limma`, `affy`, `sva`, `clusterProfiler`,
  `org.Hs.eg.db`, `enrichplot`, `fgsea`, `msigdbr`, `Seurat`, `SingleR`,
  `CellChat`, `scTenifoldKnk`, `caret`, `randomForest`, `e1071`, `glmnet`,
  `xgboost`, `gbm`, `rpart`, `kernlab`, `nnet`, `shapviz`, `pROC`, `rms`,
  `ggDCA`, `ggplot2`, `pheatmap`, `VennDiagram`. `install_packages.R` installs the
  whole set and writes the versions actually resolved on your machine to
  `package_versions.txt`.
* External programs for the steps that are not scripted: Cytoscape 3.10.3 with the
  CytoNCA plugin, AutoDock Vina 1.1.2 / AutoDockTools 1.5.7, PyMOL 2.5,
  Discovery Studio Visualizer 21.1.

## Licence and citation

Released under the MIT licence — see [`LICENSE`](LICENSE).
If you use this code, please cite the article.

## Contact

Corresponding author: De-Hui Qian — doctorqiandehui@163.com
