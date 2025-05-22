# Load necessary data
load("E:/rdk/Scissor_scenic/Scissor_data.RData")

# Define gene sets
MHC <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F",
         "HLA-DMA","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
         "HLA-DRA","HLA-DRB1","HLA-DRB5","MKI67","AURKA","TOP2A","PCNA")

# Alternative MHC gene set for some analysis (commented out original)
MHC <- c("MKI67","AURKA","TOP2A","PCNA") 

# Define hub genes for analysis
hub_gene <- c("JUNB","CREM","JUND","KLF4","MAFB","RAD21","KLF5",
              "ARNT","TRIM28","MYBL2","BMI1","EZH2","NUCKS1",
              "CEBPA","DNAJC2","FOXM1")

# Visualize expression of MHC genes using violin plot with Wilcoxon test significance
VlnPlot(Scissor_data, features = MHC) + 
  theme_minimal() +
  stat_compare_means(method = "wilcox.test", label = "p.signif")

# Filtering genes by expression threshold (commented out, for reference)
# min_cells <- round(0.005 * ncol(Scissor_data))
# Scissor_data <- subset(Scissor_data, features = rownames(Scissor_data)[
#   Matrix::rowSums(Scissor_data@assays$RNA@counts > 0) >= min_cells])

# DotPlot visualization of hub gene expression by scissor cell types
DotPlot(Scissor_data, features = hub_gene, group.by = "scissor_celltype") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0,1,0.2), 
                        colours = c('#330066','#336699','#66CC66','#FFCC33'))

dev.off()  # Close any open plotting device

# Extract metadata and set cell identities to scissor_celltype
df <- Scissor_data@meta.data
Idents(Scissor_data) <- Scissor_data$scissor_celltype

# Fetch expression data of hub genes
exprs <- data.frame(FetchData(object = Scissor_data, vars = hub_gene))

# Check expression range
range(exprs)

# Transpose data for filtering genes with too many zeros
data <- t(exprs)

# Filter genes expressed in >10% cells (for MHC genes, threshold is 40%)
non_zero_genes <- colSums(data != 0) / nrow(data) > 0.1
table(non_zero_genes)

# Keep only filtered genes
filtered_data <- data[, non_zero_genes]

# Transpose back
exprs <- t(filtered_data)

# Subset metadata for filtered expression data
df <- df[rownames(exprs),]

# Combine metadata and expression data
all <- cbind(df, exprs)

# Select columns 16 to 33 for further analysis (likely hub genes expression)
xx <- all[, c(16:33)]

# Check range of selected expression data
range(xx[, 2:14])

# Prepare data frame for multigroup visualization, removing first column
set.seed(29)
for_multigroup_vis <- data.frame(xx, 'sample_type' = xx$scissor_celltype, stringsAsFactors = FALSE)
for_multigroup_vis <- for_multigroup_vis[, -1]  # remove first column as noted

# Reshape data for ggplot2 visualization
melted_var_table <- reshape2::melt(for_multigroup_vis, 
                                  id.vars = 'sample_type', 
                                  variable.name = 'genes', 
                                  value.name = 'vst_counts')

# Set order of sample types for plotting
melted_var_table$sample_type <- factor(melted_var_table$sample_type, 
                                      levels = c('Scissor+', 'Scissor-'), 
                                      ordered = TRUE)

# Scale vst_counts for normalization (important for MHC genes)
melted_var_table$vst_counts <- scale(melted_var_table$vst_counts)

# Check scaled value range
range(melted_var_table$vst_counts)

library(dplyr)

# Perform Wilcoxon test for each gene between groups
test_results <- melted_var_table %>%
  group_by(genes) %>%
  summarize(p_value = wilcox.test(vst_counts ~ sample_type)$p.value)

# Adjust p-values for multiple testing using FDR method
test_results <- test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

# Calculate median expression per gene for sorting
summary_stats <- melted_var_table %>%
  group_by(genes) %>%
  summarize(median_expr = median(vst_counts, na.rm = TRUE)) %>%
  arrange(desc(median_expr))

# Reorder factor levels for genes based on median expression
melted_var_table$genes <- factor(melted_var_table$genes, levels = summary_stats$genes)

# Plot violin plot of gene expression by sample type with significance
multiple_vio <- ggplot(melted_var_table, aes(x = genes, y = vst_counts, fill = sample_type)) +
  geom_violin(scale = "width", alpha = 0.5, width = 0.5, size = 0.8, trim = FALSE) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.5), 
               width = 0.1, show.legend = TRUE) +
  scale_fill_manual(values = c('Scissor+' = 'orange', 'Scissor-' = '#4D85BD')) +
  stat_compare_means(aes(group = sample_type), method = "t.test", paired = FALSE,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif", label.y = 3.3, size = 4.5) +
  theme_bw() +
  theme(panel.grid = element_line(color = 'white'),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, size = 6)) +
  ylab('Expression (z-score)') +
  xlab('Genes') +
  ylim(-1, 4) +
  ggtitle('Violin plot of multiple genes')

pdf("ki67.pdf", width = 6, height = 4)
print(multiple_vio)
dev.off()

### Draw Venn Diagram
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    scissor_plus = rownames(differential_genes)[differential_genes$regulation == "Scissor+"],  # Genes up in Scissor+ group
    scissor_minus = rownames(differential_genes)[differential_genes$regulation == "Scissor-"], # Genes up in Scissor- group
    gene = genes  # NetHub Genes list
  ),
  category.names = c("Scissor+ Genes", "Scissor- Genes", "NetHub Genes"), # Labels for each set
  filename = NULL,  # Do not save file yet, draw on grid
  output = TRUE,
  fill = c("orange", "#4D85BD", "#8FBC8F"),  # Colors for each set
  alpha = 0.5,  # Transparency
  cat.col = c("orange", "#4D85BD", "#8FBC8F"), # Category label colors
  cat.cex = 1.5,  # Category label size
  cat.fontface = "bold",  # Category label font weight
  cat.pos = 0,  # Position category labels at default
  cex = 1.5,  # Set text size for gene counts
  fontfamily = "sans"
)
grid.draw(venn.plot)  # Draw Venn diagram on grid device


#### Scissor group enrichment analysis
library(BioEnricher)

# Extract Scissor+ genes with strong upregulation (log2FC > 2)
Scissorplus.genes <- differential_genes$gene[differential_genes$avg_log2FC > 2]

# Extract Scissor- genes with strong downregulation (log2FC < -2)
Scissorminus.genes <- differential_genes$gene[differential_genes$avg_log2FC < -2]

# Perform integrated Over-Representation Analysis (ORA) for Scissor+ genes
up.enrich <- lzq_ORA.integrated(
  genes = Scissorplus.genes,
  background.genes = NULL,
  GO.ont = 'ALL',
  perform.WikiPathways = TRUE,
  perform.Reactome = TRUE,
  perform.MsigDB = TRUE,
  MsigDB.category = 'H',
  perform.Cancer.Gene.Network = TRUE,
  perform.disease.ontoloty = TRUE,
  perform.DisGeNET = TRUE,
  perform.CellMarker = TRUE,
  perform.CMAP = TRUE,
  min.Geneset.Size = 3
)

# Perform integrated ORA for Scissor- genes
down.enrich <- lzq_ORA.integrated(
  genes = Scissorminus.genes,
  background.genes = NULL,
  GO.ont = 'ALL',
  perform.WikiPathways = TRUE,
  perform.Reactome = TRUE,
  perform.MsigDB = TRUE,
  MsigDB.category = 'H',
  perform.Cancer.Gene.Network = TRUE,
  perform.disease.ontoloty = TRUE,
  perform.DisGeNET = TRUE,
  perform.CellMarker = TRUE,
  perform.CMAP = TRUE,
  min.Geneset.Size = 3
)

# Save MsigDB Hallmark enrichment results as PDF with side-by-side barplots
pdf("MsigDB_H_Enrich.results.pdf", width = 6, height = 4)
lzq_ORA.barplot2(
  enrich.obj1 = up.enrich$MsigDB.H$Enrich.results,
  enrich.obj2 = down.enrich$MsigDB.H$Enrich.results,
  obj.types = c('Scissor+', 'Scissor-')
)
dev.off()


# Prepare ranked gene list for GSEA, sorted decreasingly by log2 fold change
grlist <- differential_genes$avg_log2FC
names(grlist) <- differential_genes$gene
grlist <- sort(grlist, decreasing = TRUE)

# Plot GSEA results (example)
lzq_GSEA.barplot1(enrich.obj = fit2$simplyGO)
lzq_GSEA.barplot1(enrich.obj = fit2$simplyGO, type = 'neg')

lzq_GSEA.barplot2(enrich.obj = fit2$CellMarker$Enrich.results)

lzq_gseaplot(
  fit2$simplyGO,
  Pathway.ID = 'GO:0002399',
  rank = FALSE,
  statistic.position = c(0.71, 0.85),
  rel.heights = c(1, 0.4)
)


##### Draw Inflammatory score violin + boxplot
# Must be done after AddModuleScore calculation in Seurat object

# Extract Inflammatory score and metadata
data <- FetchData(NPC.tmp, vars = c("score.Cytokines", "new_label"))

# Define order of groups so that Scissor+ and Scissor- corresponding cells are adjacent
order <- c(
  "Scissor-/B_cell", "Scissor+/B_cell",
  "Scissor-/Dendritic_cell", "Scissor+/Dendritic_cell",
  "Scissor-/Endothelial_cell", "Scissor+/Endothelial_cell",
  "Scissor-/Epithelial_cell", "Scissor+/Epithelial_cell",
  "Scissor-/Hepatocyte", "Scissor+/Hepatocyte",
  "Scissor-/iPS_cell", "Scissor+/iPS_cell",
  "Scissor-/Macrophage", "Scissor+/Macrophage",
  "Scissor-/Monocyte", "Scissor+/Monocyte",
  "Scissor-/NK_cell", "Scissor+/NK_cell",
  "Scissor-/T_cell", "Scissor+/T_cell",
  "Scissor-/Tissue_stem_cell", "Scissor+/Tissue_stem_cell"
)

# Define colors for groups: blue for Scissor-, orange for Scissor+
color_map <- c(
  "Scissor-/B_cell" = "#4D85BD", "Scissor+/B_cell" = "orange",
  "Scissor-/Dendritic_cell" = "#4D85BD", "Scissor+/Dendritic_cell" = "orange",
  "Scissor-/Endothelial_cell" = "#4D85BD", "Scissor+/Endothelial_cell" = "orange",
  "Scissor-/Epithelial_cell" = "#4D85BD", "Scissor+/Epithelial_cell" = "orange",
  "Scissor-/Hepatocyte" = "#4D85BD", "Scissor+/Hepatocyte" = "orange",
  "Scissor-/iPS_cell" = "#4D85BD", "Scissor+/iPS_cell" = "orange",
  "Scissor-/Macrophage" = "#4D85BD", "Scissor+/Macrophage" = "orange",
  "Scissor-/Monocyte" = "#4D85BD", "Scissor+/Monocyte" = "orange",
  "Scissor-/NK_cell" = "#4D85BD", "Scissor+/NK_cell" = "orange",
  "Scissor-/T_cell" = "#4D85BD", "Scissor+/T_cell" = "orange",
  "Scissor-/Tissue_stem_cell" = "#4D85BD", "Scissor+/Tissue_stem_cell" = "orange"
)

# Convert new_label to factor with defined order
data$new_label <- factor(data$new_label, levels = order)

library(ggplot2)
library(ggpubr)

# Plot violin + boxplot of score.Cytokines grouped by new_label
plot <- ggplot(data, aes(x = new_label, y = score.Cytokines, fill = new_label)) +
  geom_violin(scale = "width", alpha = 0.5, width = 0.5, size = 0.8, trim = FALSE) +  # Violin plot
  geom_boxplot(width = 0.2, position = position_dodge(0.75), outliers = FALSE, alpha = 0.7) +  # Boxplot overlay
  labs(title = "Inflammatory Score by Cell Type", x = "", y = "Score Inflammatory") +
  scale_fill_manual(values = color_map) +  # Custom colors
  theme(
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black", face = "plain", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black", face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_blank()
  ) +
  stat_compare_means(
    comparisons = list(
      c("Scissor-/B_cell", "Scissor+/B_cell"),
      c("Scissor-/Dendritic_cell", "Scissor+/Dendritic_cell"),
      c("Scissor-/Endothelial_cell", "Scissor+/Endothelial_cell"),
      c("Scissor-/Epithelial_cell", "Scissor+/Epithelial_cell"),
      c("Scissor-/Hepatocyte", "Scissor+/Hepatocyte"),
      c("Scissor-/iPS_cell", "Scissor+/iPS_cell"),
      c("Scissor-/Macrophage", "Scissor+/Macrophage"),
      c("Scissor-/Monocyte", "Scissor+/Monocyte"),
      c("Scissor-/NK_cell", "Scissor+/NK_cell"),
      c("Scissor-/T_cell", "Scissor+/T_cell"),
      c("Scissor-/Tissue_stem_cell", "Scissor+/Tissue_stem_cell")
    ),
    method = "wilcox.test",
    label = "p.signif",
    size = 4,
    hide.ns = TRUE,
    vjust = 0.1,
    bracket.size = 0.3
  )

print(plot)
ggsave("score.Inflammatory_boxplot.pdf", width = 12, height = 6)
