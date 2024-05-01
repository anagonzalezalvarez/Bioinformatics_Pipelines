# libraries
library(ggVennDiagram)
library(DOSE)
library(clusterProfiler)
library("org.Mm.eg.db")
library("pathview")


### SET WORKING SPACE ------------------------------------------------------------------------------------------------------------------------------------------------------------

analysis_name <- "functional_charcterization_Sst_Vip" 
path_imag <- paste0(basepath,"/code/",analysis_name,"/imag/")
path_dir <- paste0(basepath,"/code/",analysis_name,"/")

### Create it if it doesn't exist
dir.create(path_dir, showWarnings = FALSE)
dir.create(path_imag, showWarnings = FALSE)
setwd(path_imag)


### 1. Load LIANA outputs of final (robust) interactions  ------------------------------------------------------------------------------------------------------------------------
robust_all_df_p10 <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Final/Data/robust_all_df_p10.csv"); robust_all_df_p10 <- robust_all_df_p10[-1]
robust_all_df_p28 <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Final/Data/robust_all_df_p28.csv"); robust_all_df_p28 <- robust_all_df_p28[-1]



### 2. Create separate objects for interactions of sub-types of interest, considering the direction of interactions ---------------------------------------------------------------
SSTtoVIP_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Sst" & robust_all_df_p10$target == "Vip")
SSTtoVIP_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Sst" & robust_all_df_p28$target == "Vip")
VIPtoSST_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Vip" & robust_all_df_p10$target == "Sst")
VIPtoSST_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Vip" & robust_all_df_p28$target == "Sst")



### 3. Plotting VennDiagram to visualize the number of unique and shared interactions (across the two stages of developemnt; p10 and p28)  -------------------------------------------------------------------------------------------
#--- Sst -> Vip 
x <- list(p10 = SSTtoVIP_p10$interaction,
          p28 = SSTtoVIP_p28$interaction)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()
#----  Vip -> Sst
x <- list(p10 = VIPtoSST_p10$interaction,
          p28 = VIPtoSST_p28$interaction)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()



### 4. Subseting LIANA output for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
#--- SST to VIP
unique_pairs_SSTtoVIP_p10_inter <- Reduce(setdiff, list(A = SSTtoVIP_p10$interaction, B = SSTtoVIP_p28$interaction))
unique_pairs_SSTtoVIP_p10_df <- subset(SSTtoVIP_p10, SSTtoVIP_p10$interaction %in% unique_pairs_SSTtoVIP_p10_inter)
unique_pairs_SSTtoVIP_p28_inter <- Reduce(setdiff, list(A = SSTtoVIP_p28$interaction, B = SSTtoVIP_p10$interaction))
unique_pairs_SSTtoVIP_p28_df <- subset(SSTtoVIP_p28, SSTtoVIP_p28$interaction %in% unique_pairs_SSTtoVIP_p28_inter)
unique_pairs_SSTtoVIP_shared_inter <- intersect(SSTtoVIP_p10$interaction, SSTtoVIP_p28$interaction)
unique_pairs_SSTtoVIP_shared_df <- subset(SSTtoVIP_p28, SSTtoVIP_p28$interaction %in% unique_pairs_SSTtoVIP_shared_inter)
#--- VIP to SST 
unique_pairs_VIPtoSST_p10_inter <- Reduce(setdiff, list(A = VIPtoSST_p10$interaction, B = VIPtoSST_p28$interaction))
unique_pairs_VIPtoSST_p10_df <- subset(VIPtoSST_p10, VIPtoSST_p10$interaction %in% unique_pairs_VIPtoSST_p10_inter)
unique_pairs_VIPtoSST_p28_inter <- Reduce(setdiff, list(A = VIPtoSST_p28$interaction, B = VIPtoSST_p10$interaction))
unique_pairs_VIPtoSST_p28_df <- subset(VIPtoSST_p28, VIPtoSST_p28$interaction %in% unique_pairs_VIPtoSST_p28_inter)



### 5. Plotting dotplots for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
#--- SST to VIP
names(unique_pairs_SSTtoVIP_p10_df)[3] <- paste("ligand.complex") 
names(unique_pairs_SSTtoVIP_p10_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoVIP_p10_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))

names(unique_pairs_SSTtoVIP_shared_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoVIP_shared_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoVIP_shared_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))

names(unique_pairs_SSTtoVIP_p28_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoVIP_p28_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoVIP_p28_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"), ntop = 20)

#--- VIP to SST 
names(unique_pairs_VIPtoSST_p10_df)[3] <- paste("ligand.complex")
names(unique_pairs_VIPtoSST_p10_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_VIPtoSST_p10_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))

names(unique_pairs_VIPtoSST_p28_df)[3] <- paste("ligand.complex")
names(unique_pairs_VIPtoSST_p28_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_VIPtoSST_p28_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))



### 6. Creating vectors/objects for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
SSTtoVIP_p10_vectors <- paste(SSTtoVIP_p10$ligand, SSTtoVIP_p10$receptor, sep = "-")
SSTtoVIP_p28_vectors <- paste(SSTtoVIP_p28$ligand, SSTtoVIP_p28$receptor, sep = "-")
VIPtoSST_p10_vectors <- paste(VIPtoSST_p10$ligand, VIPtoSST_p10$receptor, sep = "-")
VIPtoSST_p28_vectors <- paste(VIPtoSST_p28$ligand, VIPtoSST_p28$receptor, sep = "-")
#--- SST to VIP
unique_pairs_SSTtoVIP_p10 <- Reduce(setdiff, list(A = SSTtoVIP_p10_vectors, B = SSTtoVIP_p28_vectors))
unique_pairs_SSTtoVIP_p28 <- Reduce(setdiff,list(A = SSTtoVIP_p28_vectors, B = SSTtoVIP_p10_vectors))
shared_pairs_SSTtoVIP <- intersect(SSTtoVIP_p10_vectors, SSTtoVIP_p28_vectors)
#--- VIP to SST 
unique_pairs_VIPtoSST_p10 <- Reduce(setdiff, list(A = VIPtoSST_p10_vectors, B = VIPtoSST_p28_vectors))
unique_pairs_VIPtoSST_p28 <- Reduce(setdiff,list(A = VIPtoSST_p28_vectors, B = VIPtoSST_p10_vectors))
shared_pairs_VIPtoSSTP <- intersect(VIPtoSST_p10_vectors, VIPtoSST_p28_vectors)



### 7. Over-representation analysis (ORA) using ligands and receptors as geneset ----------------------------------------------------------------------------------------------------------------------------------
# defining background genes 
allgenes <- unique(c(SSTtoVIP_p10$ligand, SSTtoVIP_p10$receptor,
                     SSTtoVIP_p28$ligand, SSTtoVIP_p28$receptor,
                     VIPtoSST_p10$ligand, VIPtoSST_p10$receptor,
                     VIPtoSST_p28$ligand, VIPtoSST_p28$receptor))
# defining genesets 
genes <- unique(unlist(strsplit(unique_pairs_SSTtoVIP_p10, "-"))) # change to geneset of interest 
# ORA (enrichPathway())
genes_entrez <- mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID") 
allgenes_entrez <- mapIds(org.Mm.eg.db, keys = allgenes, keytype = "SYMBOL", column = "ENTREZID") 
ePath <- enrichPathway(genes_entrez, universe = allgenes_entrez, organism = "mouse", pvalueCutoff = 1, pAdjustMethod = "none", qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500, readable = FALSE)
ePath_df <- data.frame(ePath)
barplot(ePath, showCategory=15) 
dotplot(ePath, showCategory=15) 
write.xlsx(ePath_df, file = "~/Library/Mobile Documents/com~apple~CloudDocs/Riney/SC/PC/Cell2Cell_2/All_PCs_vs_Niche/Interactions_PC_MSC/FunctionalAnnotation/MM_PCtoMSC_ePath.xlsx")

