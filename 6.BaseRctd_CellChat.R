### 6.BaseRctd_CellChat.R ###
### Workflow for integrating spatial transcriptomics (ST) and CellChat analysis --- example using HCC1 ###

# Load RCTD proportion data (HCC5 here corresponds to HCC-1N; HCC1 would correspond to HCC1-T)
raw_data <- HCC5  

# Extract proportion matrix from RCTD assay data and rename columns (replace spaces with underscores)
proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
colnames(proportion) <- gsub(' ', '_', colnames(proportion))

# Calculate the number of cells per cell type by multiplying average proportion by number of spots
cells <- round(colMeans(proportion) * nrow(proportion))

# For each cell type, mark top spots (cells) where this cell type is dominant (set to 1, else 0)
tmp <- proportion %>% top_n(cells[1], get(names(cells[1])))
raw_data$'Monocyte' <- 0
raw_data$'Monocyte'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[2], get(names(cells[2])))
raw_data$Macrophage <- 0
raw_data$Macrophage[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[3], get(names(cells[3])))
raw_data$'Endothelial-cell' <- 0
raw_data$'Endothelial-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[4], get(names(cells[4])))
raw_data$Hepatocyte <- 0
raw_data$Hepatocyte[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[5], get(names(cells[5])))
raw_data$'B-cell' <- 0
raw_data$'B-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[6], get(names(cells[6])))
raw_data$'T-cell' <- 0
raw_data$'T-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[7], get(names(cells[7])))
raw_data$'iPS-cell' <- 0
raw_data$'iPS-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[8], get(names(cells[8])))
raw_data$'Dendritic-cell' <- 0
raw_data$'Dendritic-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[9], get(names(cells[9])))
raw_data$'Epithelial-cell' <- 0
raw_data$'Epithelial-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[10], get(names(cells[10])))
raw_data$'Tissue-stem-cell' <- 0
raw_data$'Tissue-stem-cell'[rownames(tmp)] <- 1

tmp <- proportion %>% top_n(cells[11], get(names(cells[11])))
raw_data$'NK-cell' <- 0
raw_data$'NK-cell'[rownames(tmp)] <- 1

# Sum binary indicators to get number of cell types per spot
raw_data$cells_no <- (raw_data$Monocyte + raw_data$Macrophage + raw_data$`Endothelial-cell` + raw_data$Hepatocyte + raw_data$`B-cell`
                      + raw_data$`T-cell` + raw_data$`iPS-cell` + raw_data$`Dendritic-cell`
                      + raw_data$`Tissue-stem-cell` + raw_data$`NK-cell`)

# Convert any spot with more than 1 cell type to 1 (binary presence)
raw_data$cells <- raw_data$cells_no
raw_data$cells[raw_data$cells > 1] <- 1

# Keep only spots with exactly one assigned cell type (remove ambiguous spots)
data <- raw_data[, raw_data$cells == '1']

# Assign cell type labels for these spots
data$cell_types <- NA
data$cell_types[data$Monocyte == 1] <- 'Monocyte'
data$cell_types[data$Macrophage == 1] <- 'Macrophage'
data$cell_types[data$`Endothelial-cell` == 1] <- 'Endothelial-cell'
data$cell_types[data$Hepatocyte == 1] <- 'Hepatocyte'
data$cell_types[data$`B-cell` == 1] <- 'B-cell'
data$cell_types[data$`T-cell` == 1] <- 'T-cell'
data$cell_types[data$`iPS-cell` == 1] <- 'iPS-cell'
data$cell_types[data$`Dendritic-cell` == 1] <- 'Dendritic-cell'
data$cell_types[data$`Epithelial-cell` == 1] <- 'Epithelial-cell'
data$cell_types[data$`Tissue-stem-cell` == 1] <- 'Tissue-stem-cell'
data$cell_types[data$`NK-cell` == 1] <- 'NK-cell'

# Convert cell_types to factor with defined order for plotting
data$cell_types <- factor(data$cell_types, levels = c('Monocyte', 'Macrophage', 'Endothelial-cell', 'Hepatocyte', 'B-cell', 'T-cell', 'iPS-cell', 'Dendritic-cell','Epithelial-cell','Tissue-stem-cell', 'NK-cell'))

# Set identity classes and default assay for Seurat object
Idents(data) <- 'cell_types'
DefaultAssay(data) <- "SCT"

# Define colors for cell types for visualization
my_cols <- c("Monocyte" = "#927A66", "Macrophage" = "#DBAEA4", "Endothelial-cell" = "#97A4AB","Hepatocyte"="#21A69A","B-cell"="#3A6688","T-cell"="#C98882","iPS-cell" = "#A593A7","Dendritic-cell" = "#CE8662","Epithelial-cell" = "#927A86","Tissue-stem-cell" = "#927A76","NK-cell" = "#227A66")

# Plot spatial distribution of cell types and save figure
SpatialDimPlot(data, group.by = "cell_types", cols = my_cols)
ggsave("HCC-4N_RCTD.pdf", width = 6, height = 6)

### CellChat analysis preparation ###
# Prepare Seurat object HCC2 for CellChat analysis
HCC2 <- data
Idents(HCC2) <- "cell_types"

# Define color palette using scPalette function (assumed user-defined or from a package)
color.use <- scPalette(nlevels(HCC2))
names(color.use) <- levels(HCC2)

# Plot spatial plot with labels
SpatialDimPlot(HCC2, label = TRUE, label.size = 3, cols = color.use)

# Extract normalized expression matrix for CellChat input
data.input <- Seurat::GetAssayData(HCC2, slot = "data", assay = "SCT")

# Create meta data frame with cell type labels
meta <- data.frame(labels = Idents(HCC2), row.names = names(Idents(HCC2)))

# Get spatial coordinates for each spot
spatial.locs <- Seurat::GetTissueCoordinates(HCC2, scale = NULL, cols = c("imagerow", "imagecol"))

# Read scale factors from spatial data JSON file
scale.factors <- jsonlite::fromJSON(txt = file.path("E:/R/xsz/ST/RAW_DATA/HCC-1N/spatial", 'scalefactors_json.json'))

# Define spot diameters and scale factors for CellChat
scale.factors <- list(
  spot.diameter = 65, 
  spot = scale.factors$spot_diameter_fullres, 
  fiducial = scale.factors$fiducial_diameter_fullres, 
  hires = scale.factors$tissue_hires_scalef, 
  lowres = scale.factors$tissue_lowres_scalef
)

# Create CellChat object with spatial transcriptomics data
cellchat <- createCellChat(
  object = data.input, 
  meta = meta, 
  group.by = "labels", 
  datatype = "spatial", 
  coordinates = spatial.locs, 
  scale.factors = scale.factors
)

# Load CellChatDB human database (can switch to mouse DB if needed)
CellChatDB <- CellChatDB.human 

# Use full CellChat database for analysis
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# Subset data to signaling genes to reduce computation time
cellchat <- subsetData(cellchat) 

# Use parallel processing with 8 workers
future::plan("multisession", workers = 8) 

# Identify overexpressed genes and ligand-receptor interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability considering spatial distances
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, scale.distance = 0.01)

# Filter communications to those with minimum cells involved
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute communication probability at pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Visualize number and strength of interactions using circle plots
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = TRUE, label.edge= FALSE, title.name = "Interaction weights/strength")

# Heatmap visualization of communication counts and weights
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
pdf("hcc1N-pheatmap.pdf", width = 10, height = 5)
p1 + p2
dev.off()

# Select pathways to visualize
cellchat@netP$pathways
pathways.show <- c("MK")
levels(cellchat@idents)

# Choose receiver vertex indices for network visualization
vertex.receiver <- c(1, 2, 3)

# Hierarchy layout plot for selected signaling pathway
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "hierarchy")

# Circle plot for same pathway
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot of communication network
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

# Compute network centrality scores for pathway-level communication
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize centrality scores (heatmap) to identify major signaling roles
par(mfrow = c(1, 1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Bubble plot for specific source and target cell types in MK signaling pathway
pdf("hcc-1-n-C1_C0_MK.pdf", width = 5, height = 5)
netVisual_bubble(cellchat, sources.use = c(1, 5), targets.use = c(1, 4, 3), remove.isolate = FALSE, signaling = c("MK"))
dev.off()

# Violin plots for specific gene expression in single-cell data (scRNA_TN assumed Seurat object)
VlnPlot(scRNA_TN, features = c("SDC1", "SDC2", "SDC4", "MDK"), pt.size = 0, ncol = 2) +  
  stat_compare_means(method = "t.test")  # t-test for comparison

### Plot ligands and receptors expression on spatial data ###
FeaturePlot(data, features = c("SDC1", "SDC2", "SDC4", "MDK"), ncol = 4)
