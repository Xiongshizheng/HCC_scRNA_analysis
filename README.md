# HCC_scRNA_analysis
Hepatocellular Carcinoma (HCC) Single-Cell RNA-seq Analysis Pipeline
This project contains 7 R scripts for analyzing single-cell RNA sequencing data of hepatocellular carcinoma (HCC). The pipeline integrates various computational methods to dissect tumor microenvironment heterogeneity, cellular evolution, and key gene regulatory mechanisms.

Script Descriptions
1. Data Preprocessing and Integration
Performs quality control, filtering of low-quality and doublet cells, normalization, and integration of single-cell data to generate a Seurat object ready for downstream analysis.

2. Cell Type Annotation and Marker Gene Identification
Annotates major cell types based on canonical marker genes combined with Seurat clustering results, and identifies cell type-specific marker genes.

3. Scissor and SCENIC Analysis
Uses the Scissor method to link single-cell data with TCGA survival information, identifying cell subpopulations associated with prognosis. Then performs SCENIC analysis to detect key transcription factors and regulatory networks.

4. Cell-Cell Communication Network Construction
Builds cell communication networks based on ligand-receptor expression to reveal signaling interactions between immune cells and tumor cells in the tumor microenvironment.

5. Spatial Transcriptomics Integration
Integrates spatial transcriptomics data to visualize spatial distribution of key cell subpopulations, exploring their organization within the tumor microenvironment.

6. Trajectory and Pseudotime Analysis
Reconstructs dynamic evolutionary trajectories of tumor and immune cells using pseudotime analysis, identifying gene modules and regulatory mechanisms involved in cellular state transitions.

7. GeneTrajectory-Based Gene Activity Trajectories and Cell Classification
Applies the GeneTrajectory method to build gene activity trajectories and uses Garnett classifier for refined cell type annotation to further characterize hepatocyte heterogeneity.
