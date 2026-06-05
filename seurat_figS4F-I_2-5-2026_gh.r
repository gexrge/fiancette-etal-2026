library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(tibble)
library(data.table)
library(DESeq2)
options(future.globals.maxSize = 100 * 1024^3) # 1024^3 = 1Gb

# ---- settings ----
FindClusters.res <- 0.4
FindNeighbors.dims <- 1:30

# ---- read in data ----
path <- here::here()
indir <- normalizePath(file.path(path, "../../processed/seurat_droplet-qc_26-5-2026_gh_cli"))
outdir <- normalizePath(file.path(path, "../../processed/seurat_figS4F-I_2-5-2026_gh"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data.dirs <- list.dirs(indir, recursive = FALSE)
rds.dirs <- list.files(data.dirs, pattern = ".rds", full.names = TRUE)
data.list <- sapply(rds.dirs, readRDS)

name.dirs <- sapply(data.dirs, basename)
names(data.list) <- sapply(strsplit(name.dirs, "_"), `[`, 1)

# ---- Create combined seurat object ----
# Reset to raw counts for merging
data.list <- lapply(data.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x[["SCT"]] <- NULL
  return(x)
})

# Merge data for comparison 
merged <- merge(
  x = data.list[[1]],
  y = data.list[-1],
  add.cell.ids = names(data.list),
  project = "SC-vs-Human"
)

# extracting experiment from cell id
merged@meta.data$experiment <- sapply(strsplit(rownames(filtered@meta.data), "_"), "[", 1)

# extracting condition
merged@meta.data$condition <- ifelse(
  merged@meta.data$experiment == "akerman",
  sapply(strsplit(rownames(merged@meta.data), "_"), "[", 2),
  merged@meta.data$experiment
)

# ---- Start analysis ----
pdf(file.path(outdir, "seurat_SC-vs-Human_plots.pdf"), width = 10, height = 6)

merged <- SCTransform(merged, vars.to.regress = "percent.mt", verbose = FALSE)
merged <- RunPCA(merged, verbose = FALSE)
print(ElbowPlot(merged, ndims = 30))
merged <- FindNeighbors(merged, dims = FindNeighbors.dims, verbose = FALSE)
merged <- FindClusters(merged, resolution = FindClusters.res, cluster.name = "unintegrated_clusters", verbose = FALSE)
merged <- RunUMAP(merged, dims = FindNeighbors.dims, reduction.name = "umap.unintegrated", verbose = FALSE)

# ---- Plotting ----
print(DimPlot(merged, reduction = "umap.unintegrated", group.by = "unintegrated_clusters", raster = FALSE))
print(DimPlot(merged, reduction = "umap.unintegrated", group.by = "experiment", raster = FALSE))
print(DimPlot(merged, reduction = "umap.unintegrated", group.by = "condition", raster = FALSE))
print(DimPlot(merged, reduction = "umap.unintegrated", label = TRUE, raster = FALSE) + NoLegend())

gene_markers <- c(
  "ISL1", "INS", "IAPP", "GCG", "ARX", "SST", "GHRL", "PPY","TPH1",
  "KRT19", # ductal
  "PRSS1", # acinar
  "PECAM1", # vasc
  "PTPRC", # imm
  "COL1A1" # stromal
)
qc_markers <- c("MALAT1", "nFeature_RNA", "nCount_RNA", "percent.mt")
markers <- c(gene_markers, qc_markers)

for (m in markers) {
  p1 <- FeaturePlot(merged, reduction = "umap.unintegrated", features = m, raster = FALSE)
  p2 <- VlnPlot(merged, features = m) + NoLegend()
  print(p1 | p2)
}

print(DotPlot(merged, features = gene_markers))
print(RidgePlot(merged, features = gene_markers, ncol = length(gene_markers)))
print(RidgePlot(merged, features = qc_markers, ncol = length(qc_markers)))

# ---- Removing empty droplets ----
empty_clusters <- c(2, 4, 12, 21, 22)

## remove empty droplets and try UMAPing from fresh (avoid crash)
filtered <- subset(merged, !unintegrated_clusters %in% empty_clusters)
DefaultAssay(filtered) <- "RNA"
filtered[["SCT"]] <- NULL
filtered <- JoinLayers(filtered)

# re run the pipline
filtered <- SCTransform(filtered, vars.to.regress = "percent.mt", verbose = FALSE)
filtered <- RunPCA(filtered, verbose = FALSE)
print(ElbowPlot(filtered, ndims = 30))
filtered <- FindNeighbors(filtered, dims = FindNeighbors.dims, verbose = FALSE)
filtered <- FindClusters(filtered, resolution = FindClusters.res, cluster.name = "filtered_clusters", verbose = FALSE)
filtered <- RunUMAP(filtered, dims = FindNeighbors.dims, reduction.name = "umap.filtered", verbose = FALSE)

# plot new data
DimPlot(filtered, reduction = "umap.filtered", group.by = "filtered_clusters", raster = FALSE)
DimPlot(filtered, reduction = "umap.filtered", group.by = "experiment", raster = FALSE)
DimPlot(filtered, reduction = "umap.filtered", group.by = "condition", raster = FALSE)
DimPlot(filtered, reduction = "umap.filtered", label = TRUE, raster = FALSE) + NoLegend()

for (m in markers) {
  p1 <- FeaturePlot(filtered, reduction = "umap.filtered", features = m, raster = FALSE)
  ggsave(file.path(outdir, "FeaturePlots", paste0(m,".png")), plot = p1, width = 11, height = 10)
  p2 <- VlnPlot(filtered, features = m) + NoLegend()
  print(p1 | p2)
  ggsave(file.path(outdir, "FeatureVlnPlots", paste0(m, ".png")))
}

# ---- labelling clusters ----

# annotate cell types
filtered@meta.data <- filtered@meta.data %>% 
  mutate(cell_type = case_when(
    filtered$filtered_clusters %in% c(1, 9, 5, 18, 6, 23) ~ "beta",
    filtered$filtered_clusters %in% c(0, 2, 16, 17, 11, 22) ~ "alpha",
    filtered$filtered_clusters %in% c(8,19,20) ~ "delta",
    filtered$filtered_clusters %in% c(3) ~ "ductal",
    filtered$filtered_clusters %in% c(10) ~ "ec",
    filtered$filtered_clusters %in% c(7, 14) ~ "acinar",
    filtered$filtered_clusters %in% c(12) ~ "vasc",
    filtered$filtered_clusters %in% c(15) ~ "immune",
    filtered$filtered_clusters %in% c(4) ~ "stromal",
    filtered$filtered_clusters %in% c(13) ~ "ppy",
    TRUE ~ "other"
  ))

DimPlot(filtered, reduction = "umap.filtered", group.by = "cell_type", label = TRUE, raster = FALSE)

# annotate cell type origin
primary <- c("bandesh", "fasolino", "kang")
filtered$beta_cell <- case_when(
  filtered$cell_type == "beta" & filtered$experiment == "akerman" ~ "sc-beta",
  filtered$cell_type == "beta" & filtered$experiment %in% primary ~ "primary-beta",
  TRUE ~ "other"
)
filtered$alpha_cell <- case_when(
  filtered$cell_type == "alpha" & filtered$experiment == "akerman" ~ "sc-alpha",
  filtered$cell_type == "alpha" & filtered$experiment %in% primary ~ "primary-alpha",
  TRUE ~ "other"
)
filtered$delta_cell <- case_when(
  filtered$cell_type == "delta" & filtered$experiment == "akerman" ~ "sc-delta",
  filtered$cell_type == "delta" & filtered$experiment %in% primary ~ "primary-delta",
  TRUE ~ "other"
)

print(DimPlot(
  filtered, 
  reduction = "umap.filtered",
  group.by = "beta_cell",
  cols = c(
    "sc-beta" = "tomato",
    "primary-beta" = "darkblue",
    "other" = "grey"),
  raster = FALSE
))

## UMAP, all experiments but combine H1s and H3s
filtered@meta.data <- filtered@meta.data %>% 
  mutate(remi = case_when(
    condition %in% c("D45-H1", "D39-ZKSCAN1-WT37") ~ "H1",
    condition %in% c("D39-SIM1-WT8", "D39-SIM1-WT18") ~ "H3",
    TRUE ~ condition
  ))
filtered$remi <- factor(filtered$remi, levels = c("bandesh", "fasolino", "kang", "H1", "H3"))
DimPlot(filtered, reduction = "umap.filtered", group.by = "remi", raster = FALSE)

# group H1 and H3 beta cell sub populations for cleaner plotting
filtered$beta_cell_condition <- ifelse(filtered$cell_type == "beta", as.character(filtered$remi), NA)

# factorise beta_cell_type column for cleaner plotting
filtered$beta_cell_condition <- factor(
  filtered$beta_cell_condition,
  levels = c("bandesh", "fasolino", "kang", "H1", "H3")
)

print(DimPlot(
  filtered,
  reduction = "umap.filtered",
  group.by = "beta_cell_condition",
  raster = FALSE
))

# subset beta cells first to remove NA
beta <- subset(filtered, cell_type == "beta")

# recalculate SCT layer after subsetting
DefaultAssay(beta) <- "RNA"
beta[["SCT"]] <- NULL
beta <- SCTransform(beta, vars.to.regress = "percent.mt", verbose = FALSE)

# generate log normalised for violin plots for easier comparison
beta <- NormalizeData(beta, assay = "RNA", verbose = FALSE)

print(VlnPlot(
    beta, 
    assay = "RNA", 
    layer = "data", 
    group.by = "beta_cell_condition", 
    features = "INS",
    pt.size = 0,
    cols = c(
      "bandesh" = "#990000",
      "fasolino" = "#CC0000",
      "kang" = "#FF6666",
      "H1" = "#38BDF8",
      "H3" = "#003e86"
    )
  ) +
  labs(title = t) +
  NoLegend()
)

# create human vs h1 vs h3 for pairwise analysis and lickert contamination
beta$beta_cell_origin <- ifelse(
  beta$beta_cell_condition %in% primary, 
  "primary", as.character(beta$beta_cell_condition)
)

beta$beta_cell_origin <- factor(beta$beta_cell_origin)

# calculate pairwise significance values inbetween each group
pairs <- combn(levels(beta$beta_cell_origin), 2, simplify = FALSE)
done <- lapply(pairs, function(p) {
  FindMarkers(
    beta,
    assay = "RNA",
    layer = "data",
    features = "INS",
    verbose = FALSE,
    group.by = "beta_cell_origin",
    ident.1 = p[1],
    ident.2 = p[2],
    logfc.threshold = 0
  )
})
names(done) <- vapply(pairs, function(p) paste(p, collapse = "-vs-"), "")

# Actual differences between INS expr across conditions
avg.ins.expr <- AverageExpression(
  beta, 
  group.by = "beta_cell_condition", 
  assays = "RNA", 
  layer = "data",
  features = "INS", 
  verbose = FALSE
)$RNA

# AverageExpression un-logs values before averaging, re-log
avg.ins.expr_df <- avg.ins.expr %>%
  as.matrix() %>%              
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "condition", values_to = "avg.INS") %>% 
  select(-(gene)) %>% 
  mutate(avg.INS_ln = log1p(avg.INS))

# plot genes of interest
beta_markers <- c(
  "PDX1", "NKX6-1", "MAFA", "MNX1", "GLIS3", "HNF1A",  # beta cell
  "NEUROD1", "INSM1", "ISL1", "NKX2-2", "PAX6", "RFX6", "MYT1", "MYT1L", # beta cell/ neuronal
  "MAFB", "PAX4", "SIX3", "FOXO1", "FOXA2", "SOX9" # alpha cell / other 
)

# reverse order for plotting
beta$beta_cell_condition_rev <- factor(
  beta$beta_cell_condition,
  levels = rev(levels(beta$beta_cell_condition))
)

DotPlot(beta, assay = "RNA", features = beta_markers, group.by = "beta_cell_condition_rev") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(
  file.path(outdir, "DotPlots", "DotPlot_beta-markers_combined.rev.png"),
  width = 14, height = 4
)

VlnPlot(beta, assay = "RNA", layer = "data", group.by = "beta_cell_condition_rev", features = beta_markers, stack = TRUE) +
  NoLegend() 
ggsave(
  file.path(outdir, "VlnPlots", "VlnPlot_beta-markers_combined.rev.png"),
  width = 14, height = 4
)

# Load in lickert contaminants for 
lickert_contaminants <- fread(
  file.path(indir, "../..", "raw_data", "Lickert_contaminants.csv"),
  select = c("gene", "category"),
  verbose = FALSE
)

# for SBCs and PBCs, calculate pct expressed for each contaminant gene
for (bc in unique(beta$beta_cell_origin)) {
  
  beta.bc <- subset(beta, beta_cell_origin == bc)
  genes <- intersect(lickert_contaminants$gene, rownames(beta.bc))
  lc_expr <- GetAssayData(beta.bc[genes, ], assay = "RNA", layer = "data")
  
  # set pct thresholds
  for (pct in c(0, 0.5, 1, 2)) {
    
    lc_expr.pct <- rowMeans(lc_expr > pct) * 100
    col_name <- paste("pct", pct, bc, sep = "_")
    lickert_contaminants[[col_name]] <- lc_expr.pct[lickert_contaminants$gene]
    
  }
  
  remove(beta.bc, genes, lc_expr, lc_expr.pct, col_name)
}

# save and plot results
write.table(
  lickert_contaminants,
  file.path(outdir, "pheatmap", "lickert-contaminants_SBC-vs-PBC.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

lickert_contaminants <- lickert_contaminants[order(lickert_contaminants$category),]

mat <- lickert_contaminants |>
  column_to_rownames("gene") |>
  select(starts_with("pct")) |>
  as.matrix()

annotation_row <- data.frame(Category = lickert_contaminants$category)
rownames(annotation_row) <- rownames(mat)

annotation_col <- data.frame(beta_cell = sub(".*\\_", "", colnames(mat)))
rownames(annotation_col) <- colnames(mat)
  
ann_colours <- list(
  beta_cell = c(
    primary = "#CC0000",
    H1 = "#38BDF8",
    H3 = "#003e86" 
  )
)

pheatmap::pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = ann_colours,
  display_numbers = TRUE,
  #fontsize = 15,
  filename = file.path(outdir, "pheatmap", "pheatmap_pct.expr.lickert_contaminants.png"),
  width = 7,
  height = 12
)

# create new subdivisions for pcc analysis
filtered$pcc_clusters <- case_when(
  filtered$cell_type == "beta" ~ filtered$beta_cell_condition,
  filtered$cell_type == "alpha" ~ filtered$alpha_cell,
  filtered$cell_type == "delta" ~ filtered$delta_cell,
  TRUE ~ filtered$cell_type
)

filtered$pcc_clusters.reps <- case_when(
  filtered$cell_type == "beta" ~ filtered$condition,
  filtered$cell_type == "alpha" ~ filtered$alpha_cell,
  filtered$cell_type == "delta" ~ filtered$delta_cell,
  TRUE ~ filtered$cell_type
)

# factorise pcc_cluster column for cleaner plotting
filtered$pcc_clusters <- factor(
  filtered$pcc_clusters,
  levels = c(
    "bandesh", "fasolino", "kang", 
    "H1", "H3", 
    "primary-alpha", "sc-alpha", 
    "primary-delta", "sc-delta",
    "ductal"
  )
)

filtered$pcc_clusters.reps <- factor(
  filtered$pcc_clusters.reps,
  levels = c(
    "bandesh", "fasolino", "kang", 
    "D39-ZKSCAN1-WT37", "D45-H1", 
    "D39-SIM1-WT8", "D39-SIM1-WT18", 
    "primary-alpha", "sc-alpha", 
    "primary-delta", "sc-delta",
    "ductal"
  )
)

DimPlot(filtered, reduction = "umap.filtered", group.by = "pcc_clusters", raster = FALSE)
DimPlot(filtered, reduction = "umap.filtered", group.by = "pcc_clusters.reps", raster = FALSE)

# ---- Pearson correlation coefficients (pcc) ----
pcc <- subset(filtered, cell_type %in% c("beta", "alpha", "delta", "ductal"))
pcc$cell_type <- factor(pcc$cell_type)

# recalculate variable features after subsetting
DefaultAssay(pcc) <- "RNA"
pcc[["SCT"]] <- NULL
pcc <- SCTransform(pcc, assay = "RNA", verbose = FALSE)

# calculate normalised counts RNA layer for avg expr
pcc <- NormalizeData(pcc, assay = "RNA", verbose = FALSE)

# build a function for data generation
make_pcc_step <- function(step, pcc, pcc_clust, beta) {
  
  if (step == "top3000-beta-vfs") {
    
    vfs <- VariableFeatures(beta)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else if (step == "top3000-vfs") {
    
    vfs <- VariableFeatures(pcc)
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else if (step == "vfs-per-celltype") {
    
    hvgs <- lapply(levels(pcc$cell_type), function(ct) {
      obj <- subset(pcc, cell_type == ct)
      obj <- SCTransform(obj, assay = "RNA", verbose = FALSE)
      VariableFeatures(obj)
    })
    vfs <- unique(unlist(hvgs))
    
    agg <- AggregateExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else if (step == "vfs-after-downsampling") {
    
    smallest_sample <- min(table(pcc$cell_type))
    
    meta <- pcc@meta.data
    meta$cell <- rownames(meta)
    
    downsample <- meta %>%
      dplyr::group_by(cell_type) %>%
      dplyr::slice_sample(n = smallest_sample) %>%
      dplyr::pull(cell)
    
    pcc.ds <- subset(pcc, cells = downsample)
    pcc.ds <- SCTransform(pcc.ds, assay = "RNA", verbose = FALSE)
    vfs <- VariableFeatures(pcc.ds)
    
    agg <- AggregateExpression(pcc.ds, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    avg <- AverageExpression(pcc.ds, features = vfs, group.by = pcc_clust, assays = "RNA", verbose = FALSE)
    
  } else {
    stop("Unknown step")
  }
  
  ## colours
  groups <- factor(colnames(agg$RNA), levels = levels(pcc@meta.data[[pcc_clust]]))
  annotation_col <- data.frame(Group = groups)
  rownames(annotation_col) <- colnames(agg$RNA)
  
  ## correlations
  agg.corr2 <- cor(as.matrix(agg$RNA), method = "pearson")^2
  avg.corr2 <- cor(as.matrix(avg$RNA), method = "pearson")^2
  
  list(
    title = step,
    features = vfs,
    agg = agg,
    avg = avg,
    agg.corr2 = agg.corr2,
    avg.corr2 = avg.corr2,
    annotation_col = annotation_col
  )
}

# build a function for colour matching
make_pcc_colors <- function(pcc, pcc_clust) {
  
  lvls <- levels(pcc@meta.data[[pcc_clust]])
  gg_cols <- scales::hue_pal()(length(lvls))
  
  list(
    levels = lvls,
    annotation_colors = list(
      Group = setNames(gg_cols, lvls)
    )
  )
}

# build data for heatmaps
pcc_results <- list()

steps <- c(
  "top3000-beta-vfs",
  "top3000-vfs",
  "vfs-per-celltype",
  "vfs-after-downsampling"
)

for (pcc_clust in c("pcc_clusters", "pcc_clusters.reps")) {
  
  outdir_pheatmap <- file.path(outdir, "pheatmap", pcc_clust)
  dir.create(outdir_pheatmap, recursive = TRUE, showWarnings = FALSE)
  
  colors <- make_pcc_colors(pcc, pcc_clust)
  
  pcc_results[[pcc_clust]] <- list(
    outdir = outdir_pheatmap,
    colors = colors,
    steps  = lapply(
      steps, make_pcc_step,
      pcc = pcc,
      pcc_clust = pcc_clust,
      beta = beta
    )
  )
  
  names(pcc_results[[pcc_clust]]$steps) <- steps
}

# plot heatmaps
for (pcc_clust in names(pcc_results)) {
  
  if (pcc_clust == "pcc_clusters") {fon <- 15; wid <- 9; hei <- 7} 
  else {fon <- 15; wid <- 11; hei <- 8}
  
  res <- pcc_results[[pcc_clust]]
  
  for (step in names(res$steps)) {
    
    obj <- res$steps[[step]]
    
    obj$agg.corr2[lower.tri(obj$agg.corr2)] <- NA 
    pheatmap::pheatmap(
      obj$agg.corr2,
      annotation_col = obj$annotation_col,
      annotation_colors = res$colors$annotation_colors,
      main = sprintf("pearson r2 agg expr (%s)", obj$title),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      fontsize = fon,
      filename = file.path(res$outdir, sprintf("pheatmap_agg_%s.png", obj$title)),
      width = wid, height = hei,
      na_col = "white"
    )
    
    obj$avg.corr2[lower.tri(obj$avg.corr2)] <- NA
    pheatmap::pheatmap(
      obj$avg.corr2,
      annotation_col = obj$annotation_col,
      annotation_colors = res$colors$annotation_colors,
      main = sprintf("pearson r2 avg expr (%s)", obj$title),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      fontsize = fon,
      filename = file.path(res$outdir, sprintf("pheatmap_avg_%s.png", obj$title)),
      width = wid, height = hei,
      na_col = "white"
    )
  }
}

# calculate expression sets for each pcc cluster
pcc$expr.sets <- case_when(
  pcc$cell_type == "beta" ~ pcc$beta_cell,
  pcc$cell_type == "alpha" ~ pcc$alpha_cell,
  pcc$cell_type == "delta" ~ pcc$delta_cell,
  TRUE ~ NA
)

avg.expr.sets <- AverageExpression(pcc, group.by = "expr.sets", assays = "RNA", layer = "data", verbose = FALSE)$RNA
colnames(avg.expr.sets) <- sub("-", "_", colnames(avg.expr.sets))

for (coln in colnames(avg.expr.sets)) {
  avg.expr.sets %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    select(gene, coln) %>% 
    filter(.data[[coln]] > 0) %>% 
    arrange(desc(.data[[coln]])) %>% 
    write.table(
      file.path(outdir, "pcc-clusters_avg-expr-sets", paste0("avg-expr-sets_linear_", coln, ".rnk")),
      quote = FALSE,
      col.names = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
}

# ---- Perform DEA ----
if (length(filtered[["SCT"]]@SCTModel.list) > 1) {
  filtered <- PrepSCTFindMarkers(filtered, verbose = FALSE)
}

filtered.markers <- FindAllMarkers(filtered, group.by = "cell_type", verbose = FALSE)
write.table(
  filtered.markers,
  file.path(outdir, "FindAllMarkers_cell-type_output.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# find DEGs for stem-cell beta cells vs primary beta cells
filtered$beta_cell <- factor(filtered$beta_cell)
sbc.vs.pbc <- FindMarkers(
  filtered, 
  group.by = "beta_cell",
  ident.1 = "sc-beta", 
  ident.2 = "primary-beta"
)

# rename columns for clarity and create gene column
sbc.vs.pbc <- dplyr::rename(sbc.vs.pbc, pct.sbc = pct.1, pct.pbc = pct.2)
sbc.vs.pbc$gene <- row.names(sbc.vs.pbc)

write.table(
  sbc.vs.pbc,
  file.path(outdir, "FindMarkers_SBC-vs-PBC_output.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# create ranked gene set for gsea
sbc.vs.pbc %>% 
  select(gene, avg_log2FC) %>% 
  arrange(desc(avg_log2FC)) %>% 
  write.table(
    file.path(outdir, "FindMarkers_SBC-vs-PBC_output.rnk"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

# create filtered marker table
sbc.vs.pbc_sig2fc <- sbc.vs.pbc %>% 
  filter(avg_log2FC > 2 | avg_log2FC < -2, p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC))
  
write.table(
  sbc.vs.pbc_sig2fc,
  file.path(outdir, "FindMarkers_SBC-vs-PBC_output.sig2fc.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# create filtered ranked gene set
sbc.vs.pbc_sig2fc %>% 
  select(gene, avg_log2FC) %>% 
  write.table(
    file.path(outdir, "FindMarkers_SBC-vs-PBC_output.sig2fc.rnk"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

# find markers for sc and primary beta cells
filtered.markers.beta_cell <- FindAllMarkers(filtered, group.by = "beta_cell", verbose = FALSE)
filtered.markers.beta_cell <- dplyr::rename(filtered.markers.beta_cell, beta_cell = cluster)

write.table(
  filtered.markers.beta_cell,
  file.path(outdir, "FindAllMarkers_beta-cell_output.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# create preranked gene set lists for markers of SBCs and PBCs
for (bc in unique(filtered.markers.beta_cell$beta_cell)) {
  if (bc == "other") next
  filtered.markers.beta_cell %>% 
    filter(beta_cell == bc) %>% 
    select(gene, avg_log2FC) %>% 
    arrange(desc(avg_log2FC)) %>% 
    write.table(
      file.path(outdir, sprintf("FindAllMarkers_beta-cell_output.%s.rnk", bc)),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
}

# create pseudo bulk count matrix for deseq2
condition.expr <- AggregateExpression(
  beta,
  group.by = "condition",
  assay = "RNA",
  layers = "counts",
  normalization.method = NULL,
  scale.factor = NULL,
  margin = NULL,
  verbose = FALSE
)$RNA

deseq2.info <- data.frame(
  sample = colnames(condition.expr),
  condition = c("primary-beta", replicate(4, "sc-beta"), replicate(2, "primary-beta")),
  replicate = c("pri1", "sc1", "sc2", "sc3", "sc4", "pri2", "pri3")
)
rownames(deseq2.info) <- deseq2.info$sample

dds <- DESeqDataSetFromMatrix(
  countData = round(condition.expr), 
  colData = deseq2.info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "sc-beta", "primary-beta")
)

vsd <- vst(dds)
plotPCA(vsd, intgroup = "condition")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df %>%
  select(gene, log2FoldChange, padj) %>%
  filter(!is.na(padj)) %>%
  mutate(
    neglog10p = -log10(pmax(padj, 1e-300)),
    neglog10p.sign = neglog10p * ifelse(log2FoldChange >= 0, 1, -1)
  ) %>%
  select(gene, neglog10p.sign) %>%
  arrange(desc(neglog10p.sign)) %>%
  write.table(
    file.path(outdir, "DESeq2_SBC-vs-PBC.neg-log10-sign.rnk"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

# create ranked linear expression of all genes in SBCs
## aggregate expression
sbc.agg.expr <- AggregateExpression(beta, group.by = "beta_cell", assay = "RNA", verbose = FALSE)$RNA
sbc.agg.expr_df <- sbc.agg.expr %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rename("sc-beta" = "sc_beta") %>% 
  rownames_to_column(var = "gene") %>%
  select(gene, sc_beta) %>% 
  arrange(desc(sc_beta))

write.table(
  sbc.agg.expr_df,
  file.path(outdir, "AggregateExpression_SBCs.rnk"),
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

sbc.agg.expr_df %>% 
  mutate(sc_beta_int = round(sc_beta)) %>% 
  select(gene, sc_beta_int) %>% 
  write.table(
    file.path(outdir, "AggregateExpression_SBCs.int.rnk"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

## average expression
sbc.avg.expr <- AverageExpression(beta, group.by = "beta_cell", assay = "RNA", verbose = FALSE)$RNA
sbc.avg.expr_df <- sbc.avg.expr %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rename("sc-beta" = "sc_beta") %>% 
  rownames_to_column(var = "gene") %>%
  select(gene, sc_beta) %>% 
  arrange(desc(sc_beta)) %>% 
  mutate(sc_beta_log = log1p(sc_beta))

sbc.avg.expr_df %>% 
  select(gene, sc_beta) %>% 
  write.table(
    file.path(outdir, "AverageExpression_SBCs.rnk"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

sbc.avg.expr_df %>% 
  select(gene, sc_beta_log) %>% 
  write.table(
    file.path(outdir, "AverageExpression_SBCs.log1p.rnk"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
