
library(Matrix)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(clustree)
library(RColorBrewer)


basepath <- "/ibex/scratch/projects/c2169/Cecilia/CellToCell"


### SET WORKING SPACE -------------------------------------------------------------------------------------------------
analysis_name <- "1_scRNA" 
path_imag <- paste0(basepath,"/code/",analysis_name,"/imag/")
path_dir <- paste0(basepath,"/code/",analysis_name,"/")

### Create it if it doesn't exist
dir.create(path_dir, showWarnings = FALSE)
dir.create(path_imag, showWarnings = FALSE)
setwd(path_imag)


### Load the data
seurat_obj <- readRDS(paste0(path_dir,"/results/2_seurat_filtered.rds"))
head(seurat_obj@meta.data)
table(seurat_obj$sample)

### 1. Normalize Data ------------------------------------------------
seurat_obj <- NormalizeData(seurat_obj)

### 2. Identify the most variable genes -----------------------------------------------------------------
seurat_obj <- FindVariableFeatures(seurat_obj, 
                    selection.method = "vst",
                    nfeatures = 2000, 
                    verbose = FALSE)

### 3. Scale the data ------------------------------------------------
seurat_obj <- ScaleData(seurat_obj)

### 4. Run PCA ------------------------------------------------
seurat_obj <- RunPCA(seurat_obj)

### 6. Cluster ------------------------------------------------
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction='pca')
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2)) ## calculate all resolutions to be able to do treecluster

head(seurat_obj@meta.data)
### Create a new column with the chosen resolution s6_resolution and set it as the idents
seurat_obj@meta.data[['seurat_clusters']] <- seurat_obj@meta.data[['RNA_snn_res.0.8']]
Idents(object = seurat_obj) <- 'seurat_clusters'

### 7. UMAP ------------------------------------------------
seurat_obj <- RunUMAP(
  seurat_obj,
  dims=1:30,
  min.dist=0.3,
  reduction='pca',
  reduction.name = 'umap'
)


### 8. AZIMUTH ----------------------------------------------
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
seurat_obj <- Azimuth::RunAzimuth(seurat_obj, 
          reference = "mousecortexref")
head(seurat_obj@meta.data)
saveRDS(seurat_obj,paste0(path_dir,"results/3_seurat_annotated.rds"))
seurat_obj <- readRDS(paste0(path_dir,"results/3_seurat_annotated.rds"))


### 8. Plot UMAPS WITH AND without annot AZIMUTH------------------------------------------------
source("/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/1_scRNA/scripts/umaps.R")

nCount_name <- "nCount_RNA"
nFeature_name <- "nFeature_RNA"
wkg_assay <- "RNA"
s7_umap_name <- "umap"
s6_col_name_for_clusters <- "seurat_clusters"
p8_path_umap_plots <- paste0(path_imag,"3_umaps.pdf")

plotUMAPS(seurat_obj, nCount_name,nFeature_name,wkg_assay, s7_umap_name, s6_col_name_for_clusters,p8_path_umap_plots)

list <- c("P10","P28")
for (i in list){
  p8_path_umap_plots <- paste0(path_imag,"3_umaps_",i,".pdf")
  plotUMAPS(subset(seurat_obj, sample == i), nCount_name,nFeature_name,wkg_assay, s7_umap_name, s6_col_name_for_clusters,p8_path_umap_plots)
}

### 9. Save object and total cells ------------------------------------------------
head(seurat_obj@meta.data)

seurat_obj_p10 <- subset(seurat_obj, sample == "P10")
seurat_obj_p28 <- subset(seurat_obj, sample == "P28")

total_cells_per_subclass <- table(seurat_obj$predicted.subclass)
p10_total_cells_per_subclass <- table(seurat_obj_p10$predicted.subclass)
p28_total_cells_per_subclass <- table(seurat_obj_p28$predicted.subclass)

write.csv(total_cells_per_subclass, "total_cells_per_subclass_all.csv")
write.csv(p10_total_cells_per_subclass, "total_cells_per_subclass_p10.csv")
write.csv(p28_total_cells_per_subclass, "total_cells_per_subclass_p28.csv")



