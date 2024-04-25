library(Matrix)
library(Seurat)
library(stringr)
library(wesanderson)
basepath <- "/ibex/scratch/projects/c2169/Cecilia/CellToCell"


### SET WORKING SPACE -------------------------------------------------------------------------------------------------
analysis_name <- "1_scRNA" 
path_imag <- paste0(basepath,"/code/",analysis_name,"/imag/")
path_dir <- paste0(basepath,"/code/",analysis_name,"/")

### Create it if it doesn't exist
dir.create(path_dir, showWarnings = FALSE)
dir.create(path_imag, showWarnings = FALSE)
setwd(path_imag)

source(paste0(basepath,"/code/1_scRNA/scripts/qualitycontrol.R"))

### 1. Create a Seurat object for each sample ----------------------------------------------------------------------------------

### Path to my directories samples P10,P28
path <- '/ibex/scratch/projects/c2169/Cecilia/CellToCell/data/raw/'

dirs <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)
dirs
#[1] "GSM5014306_P10" "GSM5014308_P28"

### Create a Seurat object for each sample named after samples_ID
for (sample in dirs){
	### To read the mtx, features and barcodes files I need the full path
    cts <- ReadMtx(mtx = paste0(path, sample, "/", sample, "_V1_Dlxpos_RNA_matrix.mtx"), 
            features = paste0(path, sample, "/", sample, "_V1_Dlxpos_RNA_features.tsv"), 
            cells = paste0(path, sample, "/", sample, "_V1_Dlxpos_RNA_barcodes.tsv"))

    ### Create a seurat object
    assign(sample, CreateSeuratObject(counts = cts))
	}


summary(GSM5014306_P10@meta.data)
summary(GSM5014308_P28@meta.data)

head(GSM5014306_P10@meta.data)



### 2. Merge all individual samples ----------------------------------------------------------------------------------
merged_seurat <- merge(x = GSM5014306_P10, 
                      y = c(GSM5014308_P28),
                      add.cell.id = c("GSM5014306_P10" ,"GSM5014308_P28"))
merged_seurat
head(merged_seurat@meta.data)

### 3. Add more important info to metadata ----------------------------------------------------------------------------------------

# Add log10GenesPerUMI
merged_seurat$log10GenesPerUMI_RNA <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Add MT percentage
merged_seurat$mitoRatio <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create sample identifier replicate column
merged_seurat$sample <- NA
merged_seurat$sample[which(str_detect(rownames(merged_seurat@meta.data), "^GSM5014306_P10"))] <- "P10"
merged_seurat$sample[which(str_detect(rownames(merged_seurat@meta.data), "^GSM5014308_P28"))] <- "P28"
unique(merged_seurat$sample)
                 
# Create .RData object to load at any time
saveRDS(merged_seurat, paste0(path_dir,"results/1_merged_seurat.rds"))



#########################################################################################
#########################################################################################


### 4. Visualize the data without filtering  ---------------------------------

wkg_assay <- "RNA" 
col_name_to_split <- "sample"
mitoRatio_filt <- 0.01
nCount_name <- paste0("nCount_",wkg_assay)
nFeature_name <-  paste0("nFeature_",wkg_assay)
log10GenesPerUMI_name <-  paste0("log10GenesPerUMI_",wkg_assay)

colLib_samples <- wes_palette(2, name = "Chevalier1")
names(colLib_samples) = c("P10","P28")

### Check Quality Control -----------------------------------------------
path_imag_wkg <- paste0(path_imag,"1_quality_control.pdf")

pdf(path_imag_wkg)
  # violin plots for different QC stats
  features <- c(nCount_name, nFeature_name, 'mitoRatio', log10GenesPerUMI_name)

  plot_list <- lapply(features, function(x){ VlnPlot(
    merged_seurat,
    features = x,
    group.by = 'sample',
    pt.size=0) +
    RotatedAxis() +
    NoLegend() +
    geom_boxplot(notch=TRUE, fill=NA, outlier.shape=NA) +
    xlab('') +
    scale_fill_manual(values = colLib_samples)+
    theme(plot.title = element_text(size=10, hjust=0.5))
  })

  # assemble plots with patchwork
  print(wrap_plots(plot_list, ncol=2), nrow=2)
dev.off()

### 5. Filtering ------------------------------------------------

summary(merged_seurat@meta.data)


nCount_min_filt <- 2000
nCount_max_filt <- 70000
nFeature_min_filt <- 1000
nFeature_max_filt <- 8000
mitoRatio_filt <- 0.1


filtered_seurat <- subset(merged_seurat, nCount_RNA >= nCount_min_filt&
                            nCount_RNA <= nCount_max_filt & 
                            nFeature_RNA >= nFeature_min_filt & nFeature_RNA <= nFeature_max_filt &
                             mitoRatio < mitoRatio_filt)
table(filtered_seurat$sample)

### Check Quality Control Again :)
path_imag_wkg <- paste0(path_imag,"2_quality_control_1stfilter.pdf")

pdf(path_imag_wkg)
  # violin plots for different QC stats
  features <- c(nCount_name, nFeature_name, 'mitoRatio', log10GenesPerUMI_name)

  plot_list <- lapply(features, function(x){ VlnPlot(
    filtered_seurat,
    features = x,
    group.by = 'sample',
    pt.size=0) +
    RotatedAxis() +
    NoLegend() +
    geom_boxplot(notch=TRUE, fill=NA, outlier.shape=NA) +
    xlab('') +
    scale_fill_manual(values = colLib_samples)+
    theme(plot.title = element_text(size=10, hjust=0.5))
  })

  # assemble plots with patchwork
  print(wrap_plots(plot_list, ncol=2), nrow=2)
dev.off()


saveRDS(filtered_seurat, paste0(path_dir,"results/2_seurat_filtered.rds"))


unique(filtered_seurat@meta.data$sample)
summary(filtered_seurat@meta.data)

table(filtered_seurat$sample)