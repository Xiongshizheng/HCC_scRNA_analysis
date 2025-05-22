# 5.ST Processing.R
# Script for processing spatial transcriptomics data:
# - Data reading
# - Quality control and visualization
# - Normalization and clustering
# - Marker gene identification
# - Spatial deconvolution using RCTD with single-cell reference data

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyverse)
library(kableExtra)
library(RCTD)  # For spatial deconvolution
library(SeuratWrappers)  # Custom functions (make sure this script is sourced correctly)

# Set working directories
data_dir <- "D:/HCC-SC/SC-SG/data/5.ST Data/"
result_dir <- "D:/HCC-SC/SC-SG/result/5.ST/"

# Create result directory if it does not exist
if(!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

# -----------------------------
# 1. Load spatial transcriptomics data
# -----------------------------

# Read ST data samples (assuming 4 samples)
sample_ids <- c("S1","S2","S3","S4")

# Load Seurat objects into a list
st_samples <- lapply(sample_ids, function(sample) {
  sample_path <- file.path(data_dir, paste0(sample, "_filtered_feature_bc_matrix.h5"))
  Load10X_Spatial(data.dir = dirname(sample_path), filename = basename(sample_path))
})
names(st_samples) <- sample_ids

# -----------------------------
# 2. Quality Control and Visualization
# -----------------------------

for(sample in sample_ids){
  seurat_obj <- st_samples[[sample]]
  
  # Calculate percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # QC violin plot for nFeature_RNA, nCount_RNA, percent.mt
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    ggtitle(paste0(sample, " QC Metrics"))
  
  # Save QC plot
  ggsave(filename = paste0(result_dir, sample, "_QC_violin.png"))
  
  # Spatial feature plot for percent.mt
  SpatialFeaturePlot(seurat_obj, features = "percent.mt") +
    ggtitle(paste0(sample, " Spatial Mitochondrial %")) +
    theme(legend.position = "right")
  ggsave(filename = paste0(result_dir, sample, "_Spatial_percent_mt.png"))
  
  # Save back updated object
  st_samples[[sample]] <- seurat_obj
}

# -----------------------------
# 3. Data Normalization and Clustering
# -----------------------------

for(sample in sample_ids){
  seurat_obj <- st_samples[[sample]]
  
  # SCTransform normalization with mitochondrial percentage as variable to regress out
  seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", vars.to.regress = "percent.mt", verbose = FALSE)
  
  # Run PCA for dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
  
  # Find neighbors and clusters
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
  
  # Save clustering UMAP plot
  p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste0(sample, " UMAP Clusters"))
  ggsave(plot = p_umap, filename = paste0(result_dir, sample, "_UMAP_clusters.png"))
  
  # Spatial cluster plot
  p_spatial <- SpatialDimPlot(seurat_obj, group.by = "seurat_clusters") +
    ggtitle(paste0(sample, " Spatial Clusters"))
  ggsave(plot = p_spatial, filename = paste0(result_dir, sample, "_Spatial_clusters.png"))
  
  # Save updated object
  st_samples[[sample]] <- seurat_obj
}

# -----------------------------
# 4. Marker Gene Identification per Cluster
# -----------------------------

marker_list <- list()

for(sample in sample_ids){
  seurat_obj <- st_samples[[sample]]
  
  # Find cluster markers using default Wilcoxon test
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Save markers to CSV
  write.csv(markers, file = paste0(result_dir, sample, "_cluster_markers.csv"), row.names = FALSE)
  
  # Store in list for summary if needed
  marker_list[[sample]] <- markers
}

# -----------------------------
# 5. Spatial Deconvolution (RCTD)
# -----------------------------

# Load single-cell reference data (assumed preprocessed Seurat objects)
sc_reference_paths <- list(
  HCC = "D:/HCC-SC/SC-SG/data/Reference/HCC_SC.rds",
  Ma2021 = "D:/HCC-SC/SC-SG/data/Reference/Ma2021_SC.rds",
  Guilliams2022 = "D:/HCC-SC/SC-SG/data/Reference/Guilliams2022_SC.rds"
)

sc_refs <- lapply(sc_reference_paths, readRDS)

# Prepare RCTD reference objects
rctd_refs <- list()
for(ref_name in names(sc_refs)){
  sc_obj <- sc_refs[[ref_name]]
  
  # Extract count matrix and metadata
  counts <- GetAssayData(sc_obj, assay = "RNA", slot = "counts")
  cell_types <- sc_obj$celltype
  
  # Create Reference object for RCTD
  rctd_refs[[ref_name]] <- Reference(counts = counts, cell_types = cell_types)
}

# Run RCTD for each spatial sample and each reference
for(sample in sample_ids){
  seurat_obj <- st_samples[[sample]]
  
  # Extract spatial count matrix and coordinates
  spatial_counts <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts")
  coords <- GetTissueCoordinates(seurat_obj)
  
  for(ref_name in names(rctd_refs)){
    cat(paste0("Running RCTD for sample ", sample, " using reference ", ref_name, "...\n"))
    
    # Create SpatialRNA object
    spatial_rna <- SpatialRNA(coords, spatial_counts)
    
    # Create RCTD object
    rctd_obj <- create.RCTD(spatial_rna, rctd_refs[[ref_name]], max_cores = 4)
    
    # Run RCTD pipeline
    rctd_obj <- run.RCTD(rctd_obj, doublet_mode = "full")
    
    # Save results
    saveRDS(rctd_obj, file = paste0(result_dir, sample, "_RCTD_", ref_name, ".rds"))
  }
}

# End of script
