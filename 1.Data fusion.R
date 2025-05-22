# Data Fusion and Integration Evaluation Using Harmony, LIGER, and LISI

# Clear environment
rm(list = ls())

# Load required packages
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(harmony)
library(rliger)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(readr)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(future)
library(clustree)
library(reshape2)
library(lisi)
library(ggthemes)

# Set global options
options(timeout = 10000)
plan("multisession", workers = 2)
options(future.globals.maxSize = 10 * 1024^3)  # Set memory usage

# ---- 1. HARMONY Integration ----

# Load dataset (assumed predefined)
ifnb.data <- scRNA

# Preprocessing: normalization, HVGs, scaling, PCA
ifnb.data <- ifnb.data %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

# Run Harmony integration on 'orig.ident'
ifnb.data <- RunHarmony(ifnb.data, "orig.ident", plot_convergence = TRUE)

# Dimensionality reduction and clustering
n.pcs <- 20
ifnb.data <- ifnb.data %>%
  RunUMAP(reduction = "harmony", dims = 1:n.pcs, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:n.pcs, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:n.pcs) %>%
  FindClusters(resolution = 0.5, algorithm = 1)

# UMAP visualization by stage group and cluster
p1 <- DimPlot(ifnb.data, reduction = "umap", group.by = "stage_group")
p2 <- DimPlot(ifnb.data, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
P.total <- p1 + p2
ggsave(P.total, filename = "integrated_snn_res.harmony.pdf", width = 15, height = 6)

# Save integrated Seurat object
saveRDS(ifnb.data, file = "harmony.rds")

# ---- 2. LIGER Integration ----

# Reload raw data
ifnb.data <- scRNA

# Preprocessing: normalization, HVGs, scaling, PCA
ifnb.data <- ifnb.data %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE, do.center = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

# Run LIGER integration using iNMF
ifnb.data <- RunOptimizeALS(ifnb.data, k = 30, lambda = 5, split.by = "HCC")
ifnb.data <- RunQuantileNorm(ifnb.data, split.by = "HCC")

# Dimensionality reduction and clustering
n.pcs <- ncol(ifnb.data[["iNMF"]])
ifnb.data <- ifnb.data %>%
  RunUMAP(reduction = "iNMF", dims = 1:n.pcs, verbose = FALSE) %>%
  RunTSNE(reduction = "iNMF", dims = 1:n.pcs, verbose = FALSE) %>%
  FindNeighbors(reduction = "iNMF", k.param = 10, dims = 1:30) %>%
  FindClusters(resolution = 0.5, algorithm = 1)

# Visualization
p1 <- DimPlot(ifnb.data, reduction = "umap", group.by = "stage_group", raster = FALSE)
p2 <- DimPlot(ifnb.data, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = FALSE)
P.total <- p1 + p2
ggsave(P.total, filename = "integrated_snn_res.LIGER.pdf", width = 15, height = 6)

# Save integrated object
saveRDS(ifnb.data, file = "Output/integrated.LIGER.rds")

# ---- 3. LISI Evaluation ----

# Load integrated datasets
ifnb.data.liger <- read_rds("./Output/integrated.LIGER.rds")
ifnb.data.harmony <- read_rds("./Output/integrated.harmony.rds")
ifnb.data.RPCA <- read_rds("./Output/integrated.RPCA.rds")
ifnb.data.CCA <- read_rds("./Output/integrated.CCA.rds")

# LISI score calculations
pca.score <- compute_lisi(Embeddings(ifnb.data.harmony, "pca"),
                          meta_data = ifnb.data.harmony@meta.data,
                          label_colnames = "stim") %>%
  mutate(type = "raw")

harmony.score <- compute_lisi(Embeddings(ifnb.data.harmony, "harmony"),
                              meta_data = ifnb.data.harmony@meta.data,
                              label_colnames = "stage_group") %>%
  mutate(type = "harmony")

CCA.score <- compute_lisi(Embeddings(ifnb.data.CCA, "pca"),
                          meta_data = ifnb.data.CCA@meta.data,
                          label_colnames = "stim") %>%
  mutate(type = "CCA")

RPCA.score <- compute_lisi(Embeddings(ifnb.data.RPCA, "pca"),
                           meta_data = ifnb.data.RPCA@meta.data,
                           label_colnames = "stim") %>%
  mutate(type = "RPCA")

LIGER.score <- compute_lisi(Embeddings(ifnb.data.liger, "iNMF"),
                            meta_data = ifnb.data.liger@meta.data,
                            label_colnames = "stage_group") %>%
  mutate(type = "LIGER")

FastMNN.score <- compute_lisi(Embeddings(scRNA, "mnn"),
                              meta_data = scRNA@meta.data,
                              label_colnames = "stage_group") %>%
  mutate(type = "FastMNN")

BBKNN.score <- compute_lisi(Embeddings(ifnb.data, "pca"),
                            meta_data = ifnb.data@meta.data,
                            label_colnames = "stage_group") %>%
  mutate(type = "BBKNN")

# Combine results for visualization
lisi.res <- rbind(pca.score, CCA.score, RPCA.score, harmony.score, LIGER.score)
# Optional alternative:
# lisi.res <- rbind(BBKNN.score, harmony.score, FastMNN.score, LIGER.score)

# Reshape and plot
lisi.res <- tidyr::gather(lisi.res, key, val, "stage_group")
ggboxplot(lisi.res, x = "type", y = "val", color = "type", palette = "jco", short.panel.labs = FALSE) +
  labs(y = "LISI score", x = NULL, title = "Integration method comparison")
