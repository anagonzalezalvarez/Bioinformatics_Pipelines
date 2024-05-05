#### load libraries & utility function 
library(Seurat)
library(ggplot2)

# source utility functions
source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
object_path <- snakemake@input[[1]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/scrnaseq_processing_seurat/KOcall_NonTargeting/NORMALIZED_object.rds"

# outputs
umap_path <- snakemake@output[["umap_plot"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_results.csv" 
tsne_path <- snakemake@output[["tsne_plot"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_results.csv" 


### load data
data <- readRDS(file = file.path(object_path))
DefaultAssay(object = data) <- "RNA"
Idents(object = data) <- "cytokine.condition"

# plots dim
width_panel = 10
height_panel = 10


# Run non-linear dimensional reduction

# Scaling data
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes, verbose = FALSE)
seurat_object <- data

#Run PCA
# Find variable features if not already done
seurat_object <- FindVariableFeatures(seurat_object)
variable_features <- VariableFeatures(seurat_object)

seurat_object <- RunPCA(seurat_object, features = variable_features)

#seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction="pca")
#seurat_object <- FindClusters(seurat_object, resolution = 0.3)

# Run UMAP
your_seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Plot UMAP
umap_plot <- DimPlot(your_seurat_object, reduction = "umap")

# save plot
# ggsave("UMAP_plot.png", umap_path)
ggsave_new(filename = "UMAP_plot", 
           results_path=dirname(umap_path), 
           plot=umap_plot,
           width=width_panel, 
           height=height_panel)

# Run tSNE
your_seurat_object <- RunTSNE(your_seurat_object, dims = 1:10)

# #Plot tSNE
tSNE_plot <- DimPlot(your_seurat_object, reduction = "tsne")

# ggsave("tSNE_plot.png", tSNE_path)

ggsave_new(filename = "tSNE_plot", 
            results_path=dirname(tsne_path), 
            plot=tSNE_plot,
            width=width_panel, 
            height=height_panel)

