options(Seurat.object.assay.version = 'v5')
ls() 
library(dplyr)
library(Seurat)
library(patchwork)

setwd("/Users/kurowsaa/Desktop/dea_seurat-main/config/")
# Load the Memory cell data

Memory_Tcells <- readRDS("/Users/kurowsaa/Desktop/dea_seurat-main/test_data/Memory_Tcells_counts.rds")
meta_Memory <- read.csv("/Users/kurowsaa/Desktop/dea_seurat-main/test_data/Memory_Tcells_metadata.csv")
Memory_Tcells <- CreateSeuratObject(counts = Memory_Tcells, meta.data= meta_Memory, project = "Memory_Tcells", min.cells = 3, min.features = 200)
Memory_Tcells <- NormalizeData(Memory_Tcells, normalization.method = "LogNormalize", scale.factor = 10000) #normalization
Idents(Memory_Tcells) = Memory_Tcells$cytokine.condition
saveRDS(Memory_Tcells, "Memory_Tcells_updated.rds")

# Load the Naive cell data

Naive_Tcells <- readRDS("/Users/kurowsaa/Desktop/dea_seurat-main/test_data/Naive_Tcells_counts.rds")
meta_Naive <- read.csv("/Users/kurowsaa/Desktop/dea_seurat-main/test_data/Naive_Tcells_metadata.csv")
Naive_Tcells <- CreateSeuratObject(counts = Naive_Tcells, meta.data= meta_Naive, project = "Naive_Tcell", min.cells = 3, min.features = 200)
Naive_Tcells <- NormalizeData(Naive_Tcells, normalization.method = "LogNormalize", scale.factor = 10000) #normalization
Idents(Naive_Tcells) = Naive_Tcells$cytokine.condition
saveRDS(Naive_Tcells, "Naive_Tcells_updated.rds")


# Scaling data
all.genes <- rownames(Naive_Tcells)
data <- ScaleData(Naive_Tcells, features = all.genes, verbose = FALSE)

seurat_object <- data

# Find variable features if not already done
if (length(VariableFeatures(seurat_object)) == 0) {
  seurat_object <- FindVariableFeatures(seurat_object)
}

variable_features <- VariableFeatures(seurat_object)

# Ensure there are variable features identified
if (length(variable_features) > 0) {
  # Proceed with PCA
  seurat_object <- RunPCA(seurat_object, features = variable_features)
} else {
  stop("No variable features identified. PCA cannot be performed.")
}

# Now, check if UMAP has been run and proceed if not
if (!"umap" %in% Reductions(seurat_object)) {
  seurat_object <- RunUMAP(seurat_object, dims = 1:10) # Adjust dimensions as necessary
}

# Generate a basic UMAP plot without any additional modifications
umap_plot <- DimPlot(seurat_object, reduction = "umap")

# Try displaying the plot
print(umap_plot)

# Proceed to generate and save your UMAP plot
ggsave(umap_plot, filename = "UMAP_plot.png", width = 10, height = 8)

