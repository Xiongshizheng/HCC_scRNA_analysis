# Scissor&SCENIC.R
# This script integrates bulk RNA-seq survival data from TCGA-LIHC with single-cell RNA-seq data,
# identifies survival-related cell populations using Scissor, then performs SCENIC analysis
# to infer transcription factor regulatory networks in these populations.

# Load necessary libraries
library(Scissor)
library(Seurat)
library(preprocessCore)   # for quantile normalization
library(ggplot2)
library(dplyr)

# Load bulk RNA-seq data and survival data from TCGA-LIHC
surv <- read.table("TCGA-LIHC.survival.tsv", sep = "\t", row.names = 1, check.names = FALSE, header = TRUE)
exp <- read.table("TCGA-LIHC.gene_expression_TPM.tsv", sep = "\t", row.names = 1, check.names = FALSE, header = TRUE)

# Trim sample IDs to 16 characters to match survival data sample IDs
colnames(exp) <- substr(colnames(exp), 1, 16)

# Remove unwanted column from survival data and rename columns
surv <- surv[, -2]
colnames(surv) <- c("status", "time")

# Prepare data: keep intersecting samples between expression and survival data
comgene <- intersect(colnames(exp), rownames(surv))
exp <- exp[, comgene]
surv <- surv[comgene, ]

# Log-transform and quantile normalize expression matrix
tpm <- log2(exp + 1)
tpm <- as.matrix(tpm)
tpm <- normalize.quantiles(tpm)  # Quantile normalization
rownames(tpm) <- rownames(exp)
colnames(tpm) <- colnames(exp)

# Scissor analysis to associate single-cell data with survival outcome
set.seed(1234)
infos2 <- Scissor(tpm, scRNA, surv, alpha = 0.005, family = "cox", Save_file = 'Scissor_LIHC_survival.RData')

# Classify cells into Scissor+ (poor survival) and Scissor- (good survival) groups
Scissor_select <- rep("other", ncol(scRNA))
names(Scissor_select) <- colnames(scRNA)
Scissor_select[infos2$Scissor_pos] <- "Scissor+"
Scissor_select[infos2$Scissor_neg] <- "Scissor-"

# Add classification as metadata to Seurat object
scRNA <- AddMetaData(scRNA, metadata = Scissor_select, col.name = "scissor")

# Visualize Scissor cell groups with t-SNE plot
DimPlot(scRNA, reduction = 'tsne', group.by = 'scissor',
        cols = c('royalblue','indianred1','grey'), pt.size = 1.2, order = c(2,1))
ggsave("scissor.PDF", width = 6, height = 4)

# Further assign Scissor labels to celltype metadata
scRNA <- subset(scRNA, stage_group == "HCC")
x <- as.character(scRNA$celltype)
names(x) <- colnames(scRNA)
x[infos1$Scissor_pos] <- "Scissor+"
x[infos1$Scissor_neg] <- "Scissor-"
scRNA <- AddMetaData(scRNA, metadata = x, col.name = "scissor_celltype")

# Visualize scissor_celltype groups
DimPlot(scRNA, reduction = 'tsne', group.by = 'scissor_celltype', raster=FALSE, pt.size = 1.2, order = c(2,1))
ggsave("scissor_celltype.PDF", width = 6, height = 4)

# Extract Scissor+ and Scissor- cells for differential expression and SCENIC analysis
Idents(scRNA) <- scRNA$scissor_celltype
Scissor_data <- subset(scRNA, idents = c("Scissor-", "Scissor+"))
table(Scissor_data$scissor_celltype)
save(Scissor_data, file = "Scissor_data.RData")

# Find marker genes for Scissor+ and Scissor- populations
Scissor_diff <- FindAllMarkers(Scissor_data, logfc.threshold = 1.5, min.pct = 0.1, only.pos = TRUE)
Idents(Scissor_data) <- Scissor_data$scissor_celltype
Scissor_diff2 <- FindMarkers(Scissor_data, ident.1 = "Scissor-", ident.2 = "Scissor+",
                            logfc.threshold = 0.5, only.pos = TRUE,
                            test.use = "wilcox", min.pct = 0.1)

# Select top 50 markers per cluster based on average log2 fold change
Scissor_diff <- Scissor_diff %>% group_by(cluster) %>% top_n(50, wt = avg_log2FC)

write.csv(Scissor_diff, file = "ml_dt.csv")
save.image(file = "include_Scissor.RData")

# SCENIC analysis to identify transcription factor regulatory networks

library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)

# Load motif annotation database
data(list = "motifAnnotations_hgnc_v9", package = "RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

# Initialize SCENIC options with human gene database
scenicOptions <- initializeScenic(org = "hgnc",
                                  dbDir = "E:\\rdk\\Stage\\scenic\\cisTarget_databases",
                                  nCores = 1)

# Prepare expression matrix from Scissor_data (log-normalized counts)
exprMat <- as.matrix(Scissor_data@assays$RNA@data)
cellInfo <- Scissor_data@meta.data[, c("scissor_celltype", "nCount_RNA", "nFeature_RNA")]
colnames(cellInfo) <- c('CellType', 'nGene', 'nUMI')

# Filter genes based on expression criteria for SCENIC
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)

# Log-transform filtered expression matrix for GENIE3
exprMat_filtered_log <- log2(exprMat_filtered + 1)
runGenie3(exprMat_filtered_log, scenicOptions)

# Run SCENIC workflow steps
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)        # Step 1: Identify TF modules by co-expression
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod = c("top5perTarget"))  # Step 2: Motif enrichment and regulon detection
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)                # Step 3: Score regulon activity per cell
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)                                # Step 4: Binarize regulon activity for clustering

# Save SCENIC results
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")
scenicOptions <- readRDS("int/scenicOptions.Rds")

# Load regulon activity matrix and process names for compatibility
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- data.frame(t(AUCmatrix@assays@data@listData$AUC), check.names = FALSE)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(', '_', RegulonName_AUC)
RegulonName_AUC <- gsub('\\)', '', RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC

# Add AUC matrix to Seurat object metadata for visualization
pbmcauc <- AddMetaData(Scissor_data, AUCmatrix)
pbmcauc@assays$integrated <- NULL
saveRDS(pbmcauc, 'pbmcauc.rds')

# Load binary regulon activity matrix for visualization
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names = FALSE)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(', '_', RegulonName_BIN)
RegulonName_BIN <- gsub('\\)', '', RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN

pbmcbin <- AddMetaData(Scissor_data, BINmatrix)
pbmcbin@assays$integrated <- NULL
saveRDS(pbmcbin, 'pbmcbin.rds')

# Visualization of regulon activities with Seurat FeaturePlot
library(ggplot2)
dir.create('xsz_07/01_scenic_pdf')

GRNs <- intersect(colnames(AUCmatrix), colnames(BINmatrix))

for (i in 1:length(GRNs)) {
  p1 <- FeaturePlot(pbmcauc, features = GRNs[i], label = TRUE, reduction = 'tsne')
  p2 <- FeaturePlot(pbmcbin, features = GRNs[i], label = TRUE, reduction = 'tsne')
  plotc <- p1 | p2
  ggsave(paste0("xsz_07/01_scenic_pdf/", GRNs[i], ".pdf"), plotc, width = 8, height = 4)
}

