##### Load required R libraries #####
require(Seurat)
require(scales)
require(ggplot2)
require(viridis)
require(dplyr)
require(GeneTrajectory)
require(Matrix)
require(plot3D)

# GeneTrajectory constructs gene trajectories (not cell trajectories),
# which allows extracting genes that contribute continuously to a biological process
# and organizing them into trajectories revealing continuous gene activity patterns.

# Load Seurat object containing scRNA-seq data
load("E:/rdk/scRNA.RData")

# Set working directory for analysis
setwd("E:/R/xsz/GeneTrajectory")

# Check sample groups
table(scRNA$stage_group)

# Subset data to include only HCC and Healthy groups
scRNA <- subset(scRNA, stage_group == 'HCC' | stage_group == 'Healthy')

# Further subset to hepatocytes only
HCC <- subset(scRNA, celltype == "Hepatocyte")

# Plot cell clusters before reclustering
pdf("before_Hepatocyte.pdf", width = 5, height = 3)
DimPlot(HCC, group.by = "celltype", shuffle = TRUE)
dev.off()

## Re-clustering hepatocytes ##
HCC <- NormalizeData(HCC)
HCC <- FindVariableFeatures(HCC)
HCC <- ScaleData(HCC)
HCC <- RunPCA(HCC, npcs = 5, verbose = FALSE)
HCC <- RunUMAP(HCC, dims = 1:5)
HCC <- RunTSNE(HCC, dims = 1:5)
HCC <- FindNeighbors(HCC, dims = 1:5)
HCC <- FindClusters(HCC, resolution = 0.3)
HCC$cluster <- HCC$RNA_snn_res.0.3

# Visualize clusters on UMAP
DimPlot(HCC, reduction = "umap", label = TRUE, label.size = 5, group.by = "cluster")

## Garnett cell type classification setup ##
library(Seurat)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db)
library(garnett)

# Extract raw counts data and metadata from Seurat object
data <- GetAssayData(HCC, assay = 'RNA', layer = 'counts')
cell_metadata <- HCC@meta.data

# Create gene annotation dataframe
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Create monocle3 cell_data_set object for Garnett classification
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Preprocess the cell_data_set object
cds <- preprocess_cds(cds, 
                      num_dim = 5,
                      norm_method = "none")

# Load marker gene file for classification
cat(readLines("test.txt"), sep = "\n")

# Check marker genes in the dataset
marker_check <- check_markers(cds, marker_file = "test.txt",
                              db = org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

pdf("plot_markers.pdf", width = 8, height = 6)
plot_markers(marker_check)
dev.off()

# Train Garnett classifier or load existing one
class_file <- "my_classifier.Rdata"
if(!file.exists(class_file)){
  my_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = "test.txt",
                                         db = org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
  save(my_classifier, file = class_file)
} else {
  load(class_file)
}

# Assign Garnett classifications to cell metadata
pData(cds)$garnett_cluster <- pData(cds)$seurat_clusters
set.seed(1234)
cds <- classify_cells(cds,
                      my_classifier, 
                      db = org.Hs.eg.db,
                      cluster_extend = TRUE,
                      cluster_extend_max_frac_incorrect = 1,
                      cds_gene_id_type = "SYMBOL")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(viridis)

# Extract relevant metadata from cell_data_set for integration
cds.meta <- subset(pData(cds), select = c("cell_type", "cluster_ext_type")) %>% as.data.frame()

# Add Garnett classifications back to Seurat object
HCC <- AddMetaData(HCC, metadata = cds.meta)

# Plot classification results on UMAP
p1 <- DimPlot(HCC, reduction = "umap", group.by = "cell_type") + theme_bw()
p2 <- DimPlot(HCC, reduction = "umap", group.by = "cluster_ext_type") + theme_bw()

pdf("DimPlot_rescale_umap.pdf", width = 12, height = 6)
p1 + p2
dev.off()


### Metabolism activity comparison ###
devtools::install_github("wu-yc/scMetabolism")
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(ggpubr)

# Set cell identities for metabolism analysis
Idents(HCC) <- HCC$cluster_ext_type

# Compute metabolism scores with KEGG pathways using AUCell method
countexp.Seurat <- sc.metabolism.Seurat(obj = HCC, method = "AUCell", imputation = FALSE, ncores = 2, metabolism.type = "KEGG")

# Uncomment to visualize specific pathway (example: Glycolysis / Gluconeogenesis)
# DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "tsne", dimention.reduction.run = FALSE, size = 1)

# Get all metabolism pathways for downstream plotting
input.pathway <- rownames(countexp.Seurat@assays$METABOLISM$score)

# Boxplot of metabolism pathways across cluster_ext_type phenotypes
BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "cluster_ext_type", ncol = 1)

# Save DotPlot of metabolism pathways to PDF
pdf("pingfeng1.pdf", width = 7, height = 20)
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "cluster_ext_type", norm = "y")
dev.off()

# Custom function for boxplot with statistical significance annotation
BoxPlot.xsz <- function (obj, pathway, phenotype, ncol = 1) {
  input.pathway <- pathway
  input.parameter <- phenotype
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  
  metadata <- countexp.Seurat@meta.data
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
  
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, ])
  
  gg_table <- c()
  for (i in 1:length(input.pathway)) {
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], input.pathway[i], metabolism.matrix_sub[, i]))
  }
  gg_table <- data.frame(gg_table)
  gg_table[, 3] <- as.numeric(as.character(gg_table[, 3]))
  
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  
  ggplot(data = gg_table, aes(x = gg_table[, 1], y = gg_table[, 3], fill = gg_table[, 1])) +
    geom_boxplot(outlier.shape = NA) +
    geom_signif(comparisons = list(c("EP+diff", "HC meta"), c("HC meta", "MT immu"))) + 
    ylab("Metabolic Pathway") + xlab(input.parameter) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) + 
    facet_wrap(~gg_table[, 2], ncol = ncol, scales = "free") + 
    labs(fill = input.parameter) + NULL
}

### Assign hepatocyte subtype labels to Seurat object ###
# "HCC$cell_type" has unknown and malignant cells
# "HCC$cluster_ext_type" defines 4 subtypes

x <- as.data.frame(HCC$cluster_ext_type == 'HC meta')
x$orig <- rownames(x)
x <- x[x$`HCC$cluster_ext_type == "HC meta"` == TRUE,]
rownames(x)

scRNA$celltype_inclued_sub_hepa <- scRNA$celltype
cell <- as.data.frame(scRNA$celltype_inclued_sub_hepa)
colnames(cell) <- "celltype"
cell$celltype <- as.character(cell$celltype)
cell[rownames(cell) %in% rownames(x), ] <- "HC meta"

# Repeat assignment for other subtypes
for (type in c('EP+diff', 'MT immu', 'Malignance Cell')) {
  x <- as.data.frame(HCC$cluster_ext_type == type)
  x$orig <- rownames(x)
  x <- x[x[, 1] == TRUE, ]
  cell[rownames(cell) %in% rownames(x), ] <- type
}

table(cell)

# Update Seurat metadata with new hepatocyte subtypes
scRNA$celltype_inclued_sub_hepa <- as.factor(cell$celltype)

# Plot tSNE by hepatocyte subtypes
DimPlot(scRNA, label = TRUE, label.size = 2, reduction = "tsne", group.by = "celltype_inclued_sub_hepa",
        repel = TRUE, raster = FALSE) +
  theme(plot.title = element_blank(),
        axis.line = element_line(size = 0.75),
        axis.title = element_text(size = 6, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank())

ggsave("sub_Hepatocyte.pdf", width = 12, height = 9, units = "cm")


### Gene trajectory inference ###
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")

# Plot UMAP by cluster_ext_type before trajectory inference
pdf("before_GeneTrajectory_umap.pdf", height = 5, width = 5)
DimPlot(HCC, reduction = 'umap', group.by = "cluster_ext_type", cols = cbPalette) & NoAxes()
dev.off()

Idents(HCC) <- "cluster_ext_type"

# Find markers for all clusters
cluster_markers <- FindAllMarkers(HCC, only.pos = TRUE, min.diff.pct = 0.1)
cluster_markers$pct.diff <- cluster_markers$pct.1 - cluster_markers$pct.2
cluster_markers <- cluster_markers[cluster_markers$p_val_adj <= 0.05, ]

# Select top 10 markers per cluster based on average log2 fold change
top5 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Heatmap of top markers
pdf("1.pdf", width = 7, height = 5)
DoHeatmap(HCC, features = top5$gene, group.colors = cbPalette[c(3,4,2,1)], slot = "scale.data", size = 2.5) +
  theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
dev.off()

# Subset HCC cells and identify variable genes for trajectory
temp_HCC <- HCC
HCC <- subset(HCC, stage_group == 'HCC')
Idents(HCC)
DimPlot(HCC, group.by = "cluster_ext_type", shuffle = TRUE)

HCC$celltype <- HCC$cluster_ext_type
assay <- "RNA"
DefaultAssay(HCC) <- assay

HCC <- FindVariableFeatures(HCC, nfeatures = 2000)  # Top 2000 variable genes

all_genes <- HCC@assays[[assay]]@var.features
expr_percent <- apply(as.matrix(HCC[[assay]]@data[all_genes, ]) > 0, 1, sum) / ncol(HCC)

# Select genes expressed in 1% to 20% of cells
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.2)]
length(genes) # e.g., 1175

# Run gene trajectory dimensionality reduction
HCC <- GeneTrajectory::RunDM(HCC, dims = 1:5)

# Compute cell graph distances for trajectory inference
cell.graph.dist <- GetGraphDistance(HCC, dims = 1:5, K = 30)  # k and dims can be tuned

# Coarse grain cells and genes for trajectory analysis
cg_output <- CoarseGrain(HCC, cell.graph.dist, dims = 1:5, genes, N = 1000)

# Set up Python virtual environment and import trajectory function
if(!reticulate::virtualenv_exists('gene_trajectory')) {
  reticulate::virtualenv_create('gene_trajectory', packages = c('gene_trajectory'))
}
reticulate::use_virtualenv('gene_trajectory')

cal_ot_mat_from_numpy <- reticulate::import('gene_trajectory.compute_gene_distance_cmd')$cal_ot_mat_from_numpy

# Calculate gene-gene Wasserstein distance matrix
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cg_output[["graph.dist"]],
                                      gene_expr = cg_output[["gene.expression"]],
                                      num_iter_max = 50000,
                                      show_progress_bar = TRUE)

rownames(gene.dist.mat) <- cg_output[["features"]]
colnames(gene.dist.mat) <- cg_output[["features"]]
dim(gene.dist.mat)

### Gene trajectory visualization ###
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 10)$diffu.emb
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(20,21,21), K = 5)
table(gene_trajectory$selected)

gene_labels <- paste("-------", rownames(gene_embedding))
names(gene_labels) <- rownames(gene_embedding)

# Select genes for labeling on plot (you can adjust this list)
genes_selected <- top5$gene
gene_labels[!(names(gene_labels) %in% genes_selected)] <- ""

# 3D scatter plot of gene embedding colored by trajectory clusters
pdf("genetrajectory.pdf", width = 6, height = 6)
par(mar = c(1.5, 1.5, 1.5, 1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2",
          colvar = as.integer(as.factor(gene_trajectory$selected)) - 1,
          main = "Trajectory",
          pch = 19,
          cex = 1,
          theta = 45,
          phi = 0,
          col = ramp.col(c(hue_pal()(3))))
text3D(gene_embedding[,1],
       gene_embedding[,2],
       gene_embedding[,3],
       labels = gene_labels,
       add = TRUE,
       colkey = FALSE,
       cex = 0.5)
dev.off()

# Combined plots for each pseudoorder of gene trajectory
pdf("combine.pdf", width = 12, height = 6)
p1 <- scatter3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3], alpha = 1,
                bty = "b2", colvar = gene_trajectory$Pseudoorder1,
                main = "", pch = 19, cex = 1, theta = 45, phi = 0,
                col = ramp.col(rev(viridis(12)[2:11])))
p2 <- scatter3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3], alpha = 1,
                bty = "b2", colvar = gene_trajectory$Pseudoorder2,
                main = "", pch = 19, cex = 1, theta = 45, phi = 0,
                col = ramp.col(rev(viridis(12)[2:11])))
p3 <- scatter3D(gene_embedding[,1], gene_embedding[,2], gene_embedding[,3], alpha = 1,
                bty = "b2", colvar = gene_trajectory$Pseudoorder3,
                main = "", pch = 19, cex = 1, theta = 45, phi = 0,
                col = ramp.col(rev(viridis(12)[2:11])))
print(gridExtra::grid.arrange(p1, p2, p3, nrow = 1))
dev.off()

### Visualize expression of specific genes along trajectory ###
FeaturePlot(HCC, features = c("SLC25A11", "NDUFA4L2", "GLS"), reduction = "umap")

# Example: violin plot for gene "APOE" by cluster_ext_type
VlnPlot(HCC, features = "APOE", group.by = "cluster_ext_type")
