# Reviewer Comments

This directory contains scripts used for single-cell RNA-seq quality control, preprocessing, integration, and downstream analysis.

## Directory structure

reviewers-comments/
├── README.md
├── seurat_SC-vs_human_260717_gh_integrated.v6.R
├── seurat_SC-vs_human_260624_gh_integrated.v7.R
└── utils/
    ├── malat1_function.R
    ├── seurat_doubletfinder_20-2-2026_gh.r
    └── seurat_soupx_22-2-2026_gh.r

## Scripts

### seurat_SC-vs_human_260717_gh_integrated.v6.R

Main comparative Seurat analysis script for 3 primary islet datasets and 5 stem cell-derived islet datasets, including Fiancette et al.

The script performs:

- Loading and merging Seurat objects.
- Metadata assignment.
- Quality-control plotting.
- SCTransform normalization.
- PCA, clustering, and UMAP.
- scVI and FastMNN integration.
- Cell-type annotation.
- Marker and differential-expression analysis.
- Beta-cell and enteroendocrine-cell analyses.
- Summary plots and pseudobulk-style comparisons.

#### Input

The script reads processed Seurat objects from:

- ../../processed/seurat_qc_260722_gh_cli/
  
It also requires supporting files from the project, including:

- ../../raw_data/Lickert_contaminants.csv
- ../../raw_data/*.grp

#### Output

Results are written to:

- ../../processed/seurat_SC-vs-Human_260717_gh_integrated.v6/
  
The output directory contains plots, tables, heatmaps, PCA results, differential-expression results, and the final integrated Seurat object:

- integrated.v6.RDS

### seurat_SC-vs_human_260824_gh_integrated.v7.R

The QC and analysis used for `seurat_SC-vs_human_260717_gh_integrated.v6.R` was reapplied to the cohort used for the original Fiancette et al analysis, which does not include SC islets from other labs.

### utils/malat1_function.R

Defines a function for estimating a normalized MALAT1 expression threshold. The function uses density estimation and diagnostic plots to identify a threshold for filtering low-quality cells.

The same utility is also available at:

https://github.com/BaderLab/malatin/tree/v1.0.0

### utils/seurat_doubletfinder_20-2-2026_gh.r

Defines a function for detecting and removing predicted doublets from a Seurat object using DoubletFinder.

The function:

- Identifies the optimal DoubletFinder parameters.
- Estimates the expected doublet rate.
- Removes cells classified as doublets.
- Returns the filtered Seurat object and selected parameters.

The same utility is also available at:

https://github.com/chris-mcginnis-ucsf/DoubletFinder

### utils/seurat_soupx_22-2-2026_gh.r

Defines a function for correcting ambient RNA contamination using SoupX.

The function:

- Creates a SoupX object from raw and filtered counts.
- Estimates or sets the contamination fraction.
- Adjusts the expression matrix.
- Creates and returns a new Seurat object with corrected counts.
  
The same utility is also available at:

https://github.com/constantAmateur/SoupX

## Notes

The scripts contain project-specific paths, sample names, metadata mappings, and cluster labels. These may need to be updated when applying the workflow to a different dataset or directory structure.
