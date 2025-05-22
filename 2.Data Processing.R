# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scDblFinder)
library(reshape2)
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(clustree)
library(cowplot)
library(stringr)
library(SeuratDisk)
library(SeuratWrappers)

#### Filtering & QC ####
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^(MT|mt)(-|.)")  # Calculate mitochondrial gene percentage
levels(scRNA)
Idents(scRNA) <- scRNA$stage_group

# Visualize basic QC metrics
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), raster=FALSE, ncol = 3, pt.size = 0)

# Apply QC thresholds based on violin plot
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nCount_RNA < 30000 & nFeature_RNA < 6000 & percent.mt < 25)

# QC plot after filtering
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), raster=FALSE, ncol = 3, pt.size = 0)
ggsave("vlnplot.pdf", width = 15, height = 7)

save(scRNA, file = "111.Rdata")

#### Normalization ####
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#### Variable Gene Selection ####
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# Plot top 10 HVGs
top10 <- head(VariableFeatures(scRNA), 10)
LabelPoints(plot = VariableFeaturePlot(scRNA), 
            points = top10, repel = TRUE, cex = 2) +
  theme(legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.75))

ggsave("VariableFeatures_dot.pdf", width = 6, height = 9, units = "cm")

# Feature correlation
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)
FeatureScatter(scRNA, feature1 = "RNA_snn_res.0.9", feature2 = "percent.mt", raster = FALSE)

#### Data Scaling ####
# scRNA <- ScaleData(scRNA, features = rownames(scRNA)) # optional

#### Linear Dimensionality Reduction: PCA ####
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA), npcs = 50)

#### Determine the number of dimensions ####
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:20)

# JackStraw Plot
JackStrawPlot(scRNA, dims = 1:20, xmax = 0.8, ymax = 0.8) +
  theme(legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.75))
ggsave("JackStrawPlot.pdf", width = 9, height = 6, units = "cm")

# Elbow Plot to determine number of PCs
ElbowPlot(scRNA, ndims = 30) +
  theme(text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.75))
ggsave("ElbowPlot.pdf", width = 9, height = 6, units = "cm")

pca.num <- 1:20  # Set PC dimensions according to elbow plot

#### Cell Clustering ####
scRNA <- FindNeighbors(scRNA, dims = pca.num)
# Use clustree to explore resolution
clustree(scRNA@meta.data, prefix = "stage_group")
scRNA <- FindClusters(scRNA, resolution = 0.8)  # default resolution = 0.8

#### Non-linear Dimensionality Reduction ####
scRNA <- RunUMAP(scRNA, dims = pca.num)
scRNA <- RunTSNE(scRNA, dims = pca.num, check_duplicates = FALSE)

# Visualization of clusters using tSNE
DimPlot(Heal, label = TRUE, label.size = 2,
        reduction = "tsne", group.by = "seurat_clusters",
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
ggsave("seurat_clusters_tsne.pdf", width = 7, height = 6, units = "cm")

#########
#### LIGER integration ####
library()  # ??? empty call — placeholder?

# Re-normalization & PCA for LIGER
scRNA <- scRNA %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE, do.center = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

# Run LIGER: use stage_group as batch
levels(scRNA)
Idents(scRNA) <- scRNA$stage_group
scRNA <- RunOptimizeALS(scRNA, k = 30, lambda = 5, split.by = "stage_group")
scRNA <- RunQuantileNorm(scRNA, split.by = "stage_group")

#### UMAP & clustering after LIGER ####
n.pcs <- ncol(scRNA[["iNMF"]])

# Elbow Plot for iNMF dims
ElbowPlot(scRNA, ndims = 30) +
  theme(text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.75))
ggsave("ElbowPlot.pdf", width = 9, height = 6, units = "cm")
pca.num <- 1:20

# Perform UMAP, TSNE, and clustering
scRNA <- scRNA %>%
  RunUMAP(reduction = "iNMF", dims = 1:n.pcs, verbose = FALSE) %>%
  RunTSNE(reduction = "iNMF", dims = 1:n.pcs, verbose = FALSE) %>%
  FindNeighbors(reduction = "iNMF", k.param = 10, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.5, algorithm = 1)

# Check clustering resolution using clustree
clustree(scRNA@meta.data, prefix = "RNA_snn_res.")

#### UMAP visualization after LIGER ####
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "clusters", raster = FALSE)
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "stage_group", raster = FALSE)
P.total <- p1 + p2

pdf("P_umap.pdf", width = 9, height = 6)
p1
dev.off()

ggsave(P.total, filename = "integrated_LIGER_umap.pdf", width = 15, height = 6)
# saveRDS(ifnb.data,file = "Output/integrated.LIGER.rds") ### (optional save)
save(scRNA, file = "LIGER.Rdata")

#### SingleR-based automatic cell type annotation ####
library(SingleR)
library(celldex)
library(pheatmap)

# Load reference
ref <- HumanPrimaryCellAtlasData()

# Run SingleR using multiple reference datasets (assumed predefined)
pred.scRNA <- SingleR(
  test = scRNA@assays$RNA@data,
  ref = list(ref1 = ref_Hematopoietic, ref2 = ref_Human_all, ref3 = ref_Monaco),
  labels = list(ref_Hematopoietic$label.fine, ref_Human_all$label.fine, ref_Monaco$label.fine),
  clusters = scRNA@meta.data$clusters,
  fine.tune = TRUE
)

# Pruned labels
pred.scRNA$pruned.labels
heatmap <- plotScoreHeatmap(pred.scRNA, clusters = pred.scRNA@rownames, fontsize.row = 9, show_colnames = TRUE)
ggsave("singleR_heatmap.pdf", plot = heatmap, width = 18, height = 12)

# Rename clusters
new.cluster.ids1 <- pred.scRNA$pruned.labels
Idents(scRNA) <- scRNA$clusters
names(new.cluster.ids1) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids1)
scRNA@meta.data$celltype11 <- scRNA@meta.data$clusters

# UMAP & tSNE plot with celltype labels
DimPlot(scRNA, reduction = "umap", pt.size = 0.5, label = TRUE, raster = FALSE)
DimPlot(scRNA, reduction = "tsne", pt.size = 0.5, label = TRUE, raster = FALSE)
ggsave("celltype_umap.pdf", width = 8.5, height = 6)

table(scRNA$celltype)

#### Manual annotation results ####
scRNA1 <- scRNA
Idents(scRNA1) <- scRNA1$clusters
levels(scRNA1)

# Manually assigned cell types (must match number of clusters)
new.cluster.ids <- c("Monocyte", "Macrophage", "Endothelial_cell", "Hepatocyte", "B_cell",
                     "Endothelial_cell", "Hepatocyte", "Hepatocyte", "Hepatocyte", "T_cell",
                     "B_cell", "Hepatocyte", "iPS_cell", "Hepatocyte", "Dendritic_cell", "Hepatocyte",
                     "Epithelial_cell", "Tissue_stem_cell", "Hepatocyte", "T_cell",
                     "Endothelial_cell", "Hepatocyte", "Macrophage", "Hepatocyte", "NK_cell",
                     "T_cell", "Hepatocyte", "Endothelial_cell", "Tissue_stem_cell", "Tissue_stem_cell")

Idents(scRNA) <- scRNA$clusters
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$celltype <- Idents(scRNA)

saveRDS(scRNA, "scRNA_anno_liger.rds")

# Celltype UMAP
DimPlot(scRNA, reduction = "umap", label = TRUE, label.size = 2, group.by = "celltype", repel = TRUE, raster = FALSE)

# stage_group tSNE
DimPlot(scRNA_TN, reduction = "tsne", label = TRUE, label.size = 2, group.by = "stage_group", repel = TRUE, raster = FALSE) +
  theme(plot.title = element_blank(),
        axis.line = element_line(size = 0.75),
        axis.title = element_text(size = 6, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank())
ggsave("celltype_umap.pdf", width = 12, height = 9, units = "cm")
dev.off()

table(sce$celltype)

#### Find marker genes for each celltype ####
Idents(scRNA) <- "celltype"
sce.markers1 <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE)
write.csv(sce.markers1, "DEG_per_celltype.csv")

Idents(scRNA) <- "stage_group"
sce.markers2 <- FindAllMarkers(scRNA, logfc.threshold = 0.5, min.pct = 0.1, only.pos = FALSE)
write.csv(Scissor_diff, "scissor_diff.csv")

sce.markers_celltype <- sce.markers1
sce.markers_stage_group <- sce.markers2
save(sce.markers_celltype, sce.markers_stage_group, file = "findallmarkers.RData")

# Top 5 genes per cluster
sce.markers1$cluster
top5 <- sce.markers1 %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC)

# DotPlot
DotPlot(scRNA, features = unique(top2$gene), group.by = "celltype") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

DotPlot(scRNA, features = markers, group.by = "orig.ident") + RotatedAxis()
ggsave("celltype_DEG_dot.pdf", width = 6, height = 6)

### GO enrichment barplot for top5 marker genes ###
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert top 5 gene symbols to ENTREZ IDs
genelist <- bitr(top5$gene[41:45],
                 fromType = "SYMBOL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis (Biological Process)
ego_BP <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Format result
ego_result_BP <- as.data.frame(ego_BP)[1:22, ]
ego_result_BP$log_p.adjust <- (-log(ego_result_BP$p.adjust))  # Higher = more significant

go_enrich_df <- data.frame(
  ID = ego_result_BP$ID,
  Description = ego_result_BP$Description,
  log_p.adjust = ego_result_BP$log_p.adjust,
  type = factor(rep("biological process", 22), levels = "biological process"))

go_enrich_df$type_order <- factor(rev(as.integer(rownames(go_enrich_df))), labels = rev(go_enrich_df$Description))

# Plot GO enrichment barplot
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62", "#FFCC33", "#00EE76", "#FF6347", "#8A2BE2", "#CDAF95", "#EE1289")
ggplot(data = go_enrich_df, aes(x = type_order, y = log_p.adjust, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = COLS[9]) +
  coord_flip() +
  xlab("GO term") +
  ylab("-log10(p_adj)") +
  labs(title = "The Most Enriched GO Terms") +
  theme_bw()

ggsave("Fibroblasts_GO.pdf", width = 10, height = 5)
write.csv(ego_result_BP, file = "ego_result_BP_Fibroblasts.csv", row.names = TRUE)

# Violin plot for cell type marker genes
VlnPlot(scRNA, features = unique(top2$gene), stack = TRUE, pt.size = 0)
ggsave("celltype_marker_vln.pdf", width = 10, height = 6)

### GO Enrichment visualization using ClusterGVis ###
library(clusterProfiler)
# devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(org.Hs.eg.db)

head(top5)  # View top DEGs

# Prepare data for ClusterGVis
st.data <- prepareDataFromscRNA(object = scRNA,
                                diffData = top5,
                                showAverage = TRUE)

# Enrichment analysis
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.1,
                        topn = 5,
                        seed = 5314)

# Randomly select marker genes for annotation
markGenes <- unique(top5$gene)[sample(1:length(unique(top5$gene)), 35, replace = FALSE)]

# Line plot
visCluster(object = st.data, plot.type = "line")

# Full visualization with GO annotation
pdf('ClusterGVis_top5_celltype_GO.pdf', height = 10, width = 14, onefile = FALSE)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 35,
           show_row_dend = FALSE,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:7),
           go.col = rep(jjAnno::useMyCol("stallion", n = 7), each = 5),
           add.bar = TRUE)
dev.off()

# Violin plot again (legend removed)
VlnPlot(scRNA, features = unique(top5$gene), stack = TRUE, pt.size = 0) + NoLegend()
ggsave("celltype_marker_vln.pdf", width = 10, height = 6)

# Sankey / Donut / Circle plots of cell type proportions
devtools::install_github('DynamicBiosystems/cellPCT')
library(cellPCT)

# Use annotated Seurat object
scRNA <- scRNA_TN
save(scRNA_TN, file = "scRNA_TN.RData")

plot_sankey(scRNA, coord_flip = TRUE, group_by = "stage_group", cell_type = "celltype")
plot_multicircle(scRNA, group_by = "stage_group", cell_type = "celltype", outdir = '.')

# Single donut plot
pdf("ledia.pdf", width = 8, height = 5)
plot_singledonut(scRNA_TN, group_by = "stage_group", cell_type = "celltype", min_percent = 5, lab = TRUE)
dev.off()

# Subset Healthy vs HCC only
scRNA_TN <- subset(scRNA, stage_group == "HCC" | stage_group == "Healthy")

# TSNE plot of subset by cell type
DimPlot(scRNA_TN, label = TRUE, label.size = 2,
        reduction = "tsne", group.by = "celltype", repel = TRUE, raster = FALSE) +
  theme(plot.title = element_blank(),
        axis.line = element_line(size = 0.75),
        axis.title = element_text(size = 6, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank())
ggsave("celltype_umap.pdf", width = 12, height = 9, units = "cm")

#### Heatmap by cell type and condition ####
library(Seurat)
library(pheatmap)

# Correlation by stage group
Idents(scRNA_TN) <- scRNA_TN$stage_group
av.exp <- AverageExpression(scRNA_TN)$RNA
features <- names(tail(sort(apply(av.exp, 1, sd)), 2000))
av.exp <- av.exp[which(row.names(av.exp) %in% features), ]
av.exp <- as.data.frame(av.exp)
av.exp <- cor(av.exp, method = "spearman")
pheatmap::pheatmap(av.exp)

# Correlation by stage+celltype
library(tidyr)
scRNA_TN@meta.data <- unite(scRNA_TN@meta.data, "stage_celltype", stage_group, celltype, remove = FALSE)
Idents(scRNA_TN) <- scRNA_TN$stage_celltype
exp <- AverageExpression(scRNA_TN)$RNA
features <- names(tail(sort(apply(exp, 1, sd)), 2000))
exp <- exp[which(row.names(exp) %in% features), ]
exp <- as.data.frame(exp)
exp <- cor(exp, method = "spearman")

# Annotations for heatmap
annotation_col <- data.frame(
  celltype = c("Monocyte", "Macrophage", "Endothelial_cell", "Hepatocyte", "B_cell", "T_cell",
               "iPS_cell", "Dendritic_cell", "Epithelial_cell", "Tissue_stem_cell", "NK_cell"),
  stage_group = c(rep("HCC", 11), rep("Healthy", 11))
)
row.names(annotation_col) <- colnames(exp)

annotation_row <- annotation_col
row.names(annotation_row) <- rownames(exp)

pheatmap::pheatmap(exp, annotation_col = annotation_col, annotation_row = annotation_row,
                   color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
library(ClusterGVis)
Idents(scRNA_TN) <- scRNA_TN$celltype

# DEG for each celltype
pbmc.markers.all <- Seurat::FindAllMarkers(scRNA_TN,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

# Top 5 markers per cluster
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC)

# Prepare for ClusterGVis
st.data <- prepareDataFromscRNA(object = scRNA_TN,
                                diffData = pbmc.markers,
                                showAverage = TRUE)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.05,
                        topn = 3,
                        seed = 5201314)

head(enrich)

markGenes <- unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)), 55, replace = FALSE)]

# Visualization
visCluster(object = st.data, plot.type = "line")

pdf("stage.pdf", width = 13, height = 11)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 10,
           show_row_dend = FALSE,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:11),
           go.col = rep(jjAnno::useMyCol("stallion", n = 11), each = 3),
           add.bar = TRUE)
dev.off()
