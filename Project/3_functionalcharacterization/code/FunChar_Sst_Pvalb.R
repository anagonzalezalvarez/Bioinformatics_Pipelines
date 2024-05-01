# libraries
# install.packages("ggVennDiagram")
# BiocManager::install("pathview")

library(ggVennDiagram)
library(DOSE)
library(clusterProfiler)
library("org.Mm.eg.db")
library("pathview")
library(ggplot2)


### SET WORKING SPACE ------------------------------------------------------------------------------------------------------------------------------------------------------------
basepath <- "/Users/gonzalac/Desktop/PhD_2nd/BESE398_Pipelines/Bioinformatics_Pipelines/Project/"
analysis_name <- "3_functionalcharacterization" 
path_imag <- paste0(basepath,analysis_name,"/imag/")
path_dir <- paste0(basepath,analysis_name,"/")

### Create it if it doesn't exist
dir.create(path_dir, showWarnings = FALSE)
dir.create(path_imag, showWarnings = FALSE)
setwd(path_imag)

### 1. Load LIANA outputs of final (robust) interactions  ------------------------------------------------------------------------------------------------------------------------
robust_all_df_p10 <- read.csv(paste0(basepath, "2_liana/results/robust_all_df_p10.csv")); robust_all_df_p10 <- robust_all_df_p10[-1]
robust_all_df_p28 <- read.csv(paste0(basepath, "2_liana/results/robust_all_df_p28.csv")); robust_all_df_p28 <- robust_all_df_p28[-1]

head(robust_all_df_p10)

### 2. Create separate objects for interactions of sub-types of interest, considering the direction of interactions ---------------------------------------------------------------
SSTtoPVALB_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Sst" & robust_all_df_p10$target == "Pvalb")
SSTtoPVALB_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Sst" & robust_all_df_p28$target == "Pvalb")
PVALBtoSST_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Pvalb" & robust_all_df_p10$target == "Sst")
PVALBtoSST_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Pvalb" & robust_all_df_p28$target == "Sst")



### 3. Plotting VennDiagram to visualize the number of unique and shared interactions (across the two stages of developemnt; p10 and p28)  -------------------------------------------------------------------------------------------
#--- Sst -> Pvalb 
pdf("1_vendiagram_interactions_SSTtoPVALB.pdf")
x <- list(p10 = SSTtoPVALB_p10,
          p28 = SSTtoPVALB_p28)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()

#--- Pvalb -> Sst
x <- list(p10 = SSTtoPVALB_p10,
          p28 = SSTtoPVALB_p28)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()
dev.off()


### 4. Subseting LIANA output for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
#--- SST to PVALB
unique_pairs_SSTtoPVALB_p10_inter <- Reduce(setdiff, list(A = SSTtoPVALB_p10$interaction, B = SSTtoPVALB_p28$interaction))
unique_pairs_SSTtoPVALB_p10_df <- subset(SSTtoPVALB_p10, SSTtoPVALB_p10$interaction %in% unique_pairs_SSTtoPVALB_p10_inter)
unique_pairs_SSTtoPVALB_p28_inter <- Reduce(setdiff, list(A = SSTtoPVALB_p28$interaction, B = SSTtoPVALB_p10$interaction))
unique_pairs_SSTtoPVALB_p28_df <- subset(SSTtoPVALB_p28, SSTtoPVALB_p28$interaction %in% unique_pairs_SSTtoPVALB_p28_inter)
unique_pairs_SSTtoPVALB_shared_inter <- intersect(SSTtoPVALB_p10$interaction, SSTtoPVALB_p28$interaction)
unique_pairs_SSTtoPVALB_shared_df <- subset(SSTtoPVALB_p28, SSTtoPVALB_p28$interaction %in% unique_pairs_SSTtoPVALB_shared_inter)
#--- PVALB to SST 
unique_pairs_PVALBtoSST_p10_inter <- Reduce(setdiff, list(A = PVALBtoSST_p10$interaction, B = PVALBtoSST_p28$interaction))
unique_pairs_PVALBtoSST_p10_df <- subset(PVALBtoSST_p10, PVALBtoSST_p10$interaction %in% unique_pairs_PVALBtoSST_p10_inter)
unique_pairs_PVALBtoSST_p28_inter <- Reduce(setdiff, list(A = PVALBtoSST_p28$interaction, B = PVALBtoSST_p10$interaction))
unique_pairs_PVALBtoSST_p28_df <- subset(PVALBtoSST_p28, PVALBtoSST_p28$interaction %in% unique_pairs_PVALBtoSST_p28_inter)
unique_pairs_PVALBtoSST_shared_inter <- intersect(PVALBtoSST_p10$interaction, PVALBtoSST_p28$interaction)
unique_pairs_PVALBtoSST_shared_df <- subset(PVALBtoSST_p28, PVALBtoSST_p28$interaction %in% unique_pairs_PVALBtoSST_shared_inter)



### 5. Plotting dotplots for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------

#--- SST to Pvalb
names(unique_pairs_SSTtoPVALB_p10_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoPVALB_p10_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoPVALB_p10_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Pvalb"), target_groups =  c("Sst", "Pvalb"))

names(unique_pairs_SSTtoPVALB_shared_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoPVALB_shared_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoPVALB_shared_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Pvalb"), target_groups =  c("Sst", "Pvalb"))

names(unique_pairs_SSTtoPVALB_p28_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoPVALB_p28_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoPVALB_p28_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Pvalb"), target_groups =  c("Sst", "Pvalb"), ntop = 20)

#--- Pvalb to SST 
names(unique_pairs_PVALBtoSST_p10_df)[3] <- paste("ligand.complex")
names(unique_pairs_PVALBtoSST_p10_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_PVALBtoSST_p10_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Pvalb"), target_groups =  c("Sst", "Pvalb")) 

names(unique_pairs_PVALBtoSST_shared_df)[3] <- paste("ligand.complex")
names(unique_pairs_PVALBtoSST_shared_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_PVALBtoSST_shared_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Pvalb"), target_groups =  c("Sst", "Pvalb")) 

names(unique_pairs_PVALBtoSST_p28_df)[3] <- paste("ligand.complex")
names(unique_pairs_PVALBtoSST_p28_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_PVALBtoSST_p28_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Pvalb"), target_groups =  c("Sst", "Pvalb"))




### 6. Creating vectors/objects for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
SSTtoPVALB_p10_vectors <- paste(SSTtoPVALB_p10$ligand, SSTtoPVALB_p10$receptor, sep = "-")
SSTtoPVALB_p28_vectors <- paste(SSTtoPVALB_p28$ligand, SSTtoPVALB_p28$receptor, sep = "-")
PVALBtoSST_p10_vectors <- paste(PVALBtoSST_p10$ligand, PVALBtoSST_p10$receptor, sep = "-")
PVALBtoSST_p28_vectors <- paste(PVALBtoSST_p28$ligand, PVALBtoSST_p28$receptor, sep = "-")

#--- Sst -> Pvalb 
unique_pairs_SSTtoPVALB_p10 <- Reduce(setdiff, list(A = SSTtoPVALB_p10_vectors, B = SSTtoPVALB_p28_vectors))
unique_pairs_SSTtoPVALB_p28 <- Reduce(setdiff,list(A = SSTtoPVALB_p28_vectors, B = SSTtoPVALB_p10_vectors))
shared_pairs_SSTtoPVALB <- intersect(SSTtoPVALB_p10_vectors, SSTtoPVALB_p28_vectors)
#--- Pvalb -> Sst
unique_pairs_PVALBtoSST_p10 <- Reduce(setdiff, list(A = PVALBtoSST_p10_vectors, B = PVALBtoSST_p28_vectors))
unique_pairs_PVALBtoSST_p28 <- Reduce(setdiff,list(A = PVALBtoSST_p28_vectors, B = PVALBtoSST_p10_vectors))
shared_pairs_PVALBtoSSTP <- intersect(PVALBtoSST_p10_vectors, PVALBtoSST_p28_vectors)



### 7. Over-representation analysis (ORA) using ligands and receptors as geneset ----------------------------------------------------------------------------------------------------------------------------------
# defining background genes 
allgenes <- unique(c(SSTtoPVALB_p10$ligand, SSTtoPVALB_p10$receptor,
                     SSTtoPVALB_p28$ligand, SSTtoPVALB_p28$receptor,
                     PVALBtoSST_p10$ligand, PVALBtoSST_p10$receptor,
                     PVALBtoSST_p28$ligand, PVALBtoSST_p28$receptor))

# defining genesets 
genes <- unique(unlist(strsplit(unique_pairs_SSTtoPVALB_p10, "-"))) # change to geneset of interest 
# ORA (enrichPathway())
genes_entrez <- mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID") 
allgenes_entrez <- mapIds(org.Mm.eg.db, keys = allgenes, keytype = "SYMBOL", column = "ENTREZID") 
ePath <- enrichPathway(genes_entrez, universe = allgenes_entrez, organism = "mouse", pvalueCutoff = 1, pAdjustMethod = "none", qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500, readable = FALSE)
ePath_df <- data.frame(ePath)
barplot(ePath, showCategory=15) 
dotplot(ePath, showCategory=15) 
write.xlsx(ePath_df, file = "~/Library/Mobile Documents/com~apple~CloudDocs/Riney/SC/PC/Cell2Cell_2/All_PCs_vs_Niche/Interactions_PC_MSC/FunctionalAnnotation/MM_PCtoMSC_ePath.xlsx")

