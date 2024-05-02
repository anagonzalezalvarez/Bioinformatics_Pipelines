# libraries
# install.packages("ggVennDiagram")
BiocManager::install("pathview")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
BiocManager::install("ReactomePA")

remotes::install_github('saezlab/liana')
library(ggVennDiagram)
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(ggplot2)
library(liana)
library(ReactomePA)

### SET WORKING SPACE ------------------------------------------------------------------------------------------------------------------------------------------------------------
basepath <- "/ibex/scratch/projects/c2169/Cecilia/Bioinformatics_Pipelines/Project/"
analysis_name <- "3_functionalcharacterization" 
path_imag <- paste0(basepath,analysis_name,"/imag/SST_VIP")
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
SSTtoVIP_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Sst" & robust_all_df_p10$target == "Vip")
SSTtoVIP_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Sst" & robust_all_df_p28$target == "Vip")
VIPtoSST_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Vip" & robust_all_df_p10$target == "Sst")
VIPtoSST_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Vip" & robust_all_df_p28$target == "Sst")


### 3. Plotting VennDiagram to visualize the number of unique and shared interactions (across the two stages of developemnt; p10 and p28)  -------------------------------------------------------------------------------------------
#--- Sst -> Vip 
x <- list(p10 = SSTtoVIP_p10$interaction,
          p28 = SSTtoVIP_p28$interaction)
plot1 <- ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()
#----  Vip -> Sst
x <- list(p10 = VIPtoSST_p10$interaction,
          p28 = VIPtoSST_p28$interaction)
plot2 <- ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()

# Save the plots to PDF
ggsave("1_vendiagram_interactions_SSTtoVIP_plot1.pdf", plot1, width = 8, height = 6)
ggsave("1_vendiagram_interactions_VIPtoSST_plot2.pdf", plot2, width = 8, height = 6)


### 4. Subseting LIANA output for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
#--- SST to VIP

## UNIQUE to p10 (left side of ven plot = 2)
unique_pairs_SSTtoVIP_p10_inter <- Reduce(setdiff, list(A = SSTtoVIP_p10$interaction, B = SSTtoVIP_p28$interaction))
unique_pairs_SSTtoVIP_p10_df <- subset(SSTtoVIP_p10, SSTtoVIP_p10$interaction %in% unique_pairs_SSTtoVIP_p10_inter)

## UNIQUE to p28 (right side of ven = 31)
unique_pairs_SSTtoVIP_p28_inter <- Reduce(setdiff, list(A = SSTtoVIP_p28$interaction, B = SSTtoVIP_p10$interaction))
unique_pairs_SSTtoVIP_p28_df <- subset(SSTtoVIP_p28, SSTtoVIP_p28$interaction %in% unique_pairs_SSTtoVIP_p28_inter)

## SHARED (middle of ven = 3)
unique_pairs_SSTtoVIP_shared_inter <- intersect(SSTtoVIP_p10$interaction, SSTtoVIP_p28$interaction)
unique_pairs_SSTtoVIP_shared_df <- subset(SSTtoVIP_p28, SSTtoVIP_p28$interaction %in% unique_pairs_SSTtoVIP_shared_inter)

#--- VIP to SST 
## UNIQUE to p10 (left side of ven plot = 3)
unique_pairs_VIPtoSST_p10_inter <- Reduce(setdiff, list(A = VIPtoSST_p10$interaction, B = VIPtoSST_p28$interaction))
unique_pairs_VIPtoSST_p10_df <- subset(VIPtoSST_p10, VIPtoSST_p10$interaction %in% unique_pairs_VIPtoSST_p10_inter)

## UNIQUE to p28 (right side of ven plot = 18)
unique_pairs_VIPtoSST_p28_inter <- Reduce(setdiff, list(A = VIPtoSST_p28$interaction, B = VIPtoSST_p10$interaction))
unique_pairs_VIPtoSST_p28_df <- subset(VIPtoSST_p28, VIPtoSST_p28$interaction %in% unique_pairs_VIPtoSST_p28_inter)

## in this direction no SHARED interactions


### 5. Plotting dotplots for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
#--- SST to VIP

pdf("2_dotplot_interactions_SSTtoVIP.pdf")

names(unique_pairs_SSTtoVIP_p10_df)[3] <- paste("ligand.complex")  ### rename cols ligand to ligand.complex
names(unique_pairs_SSTtoVIP_p10_df)[4] <- paste("receptor.complex") 
liana_dotplot(unique_pairs_SSTtoVIP_p10_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))

names(unique_pairs_SSTtoVIP_shared_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoVIP_shared_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoVIP_shared_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))

names(unique_pairs_SSTtoVIP_p28_df)[3] <- paste("ligand.complex")
names(unique_pairs_SSTtoVIP_p28_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_SSTtoVIP_p28_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"), ntop = 20)

dev.off()
#--- VIP to SST 

pdf("2_dotplot_interactions_VIPtoSST.pdf")
names(unique_pairs_VIPtoSST_p10_df)[3] <- paste("ligand.complex")
names(unique_pairs_VIPtoSST_p10_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_VIPtoSST_p10_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))

names(unique_pairs_VIPtoSST_p28_df)[3] <- paste("ligand.complex")
names(unique_pairs_VIPtoSST_p28_df)[4] <- paste("receptor.complex")
liana_dotplot(unique_pairs_VIPtoSST_p28_df, specificity = "cellphonedb.pvalue", source_groups = c("Sst", "Vip"), target_groups =  c("Sst", "Vip"))
dev.off()




### 7. Over-representation analysis (ORA) using ligands and receptors as geneset ----------------------------------------------------------------------------------------------------------------------------------
# defining background genes 
allgenes <- unique(c(SSTtoVIP_p10$ligand, SSTtoVIP_p10$receptor,
                    SSTtoVIP_p28$ligand, SSTtoVIP_p28$receptor,
                    VIPtoSST_p10$ligand, VIPtoSST_p10$receptor,
                    VIPtoSST_p28$ligand, VIPtoSST_p28$receptor))
allgenes_entrez <- mapIds(org.Mm.eg.db, keys = allgenes, keytype = "SYMBOL", column = "ENTREZID") 



# Function to perform ORA on a set of gene pairs
perform_ORA <- function(gene_pairs, file_name) {
  genes <- unique(unlist(strsplit(gene_pairs, "-")))
  genes_entrez <- mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  
  ePath <- enrichPathway(genes_entrez, universe = allgenes_entrez, organism = "mouse", pvalueCutoff = 1, pAdjustMethod = "none", qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500, readable = FALSE)
  ePath_df <- data.frame(ePath)
  
  barplot <- barplot(ePath, showCategory = 15)
  dotplot <- dotplot(ePath, showCategory = 15)

  # Save plots using ggsave
  ggsave(paste0("3_", file_name, "_bar_plot.pdf"), barplot, width = 10, height = 8)
  ggsave(paste0("3_", file_name, "_dot_plot.pdf"), dotplot, width = 10, height = 8)
}



### Creating vectors/objects for unique and shared interactions (across the two stages of developemnt; p10 and p28) -------------------------------------------------------------------------------------------
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



# Applying the function to each gene pair set
perform_ORA(unique_pairs_SSTtoVIP_p10, "Unique_SSTtoVIP_p10_ORA")
perform_ORA(unique_pairs_SSTtoVIP_p28, "Unique_SSTtoVIP_p28_ORA")
perform_ORA(shared_pairs_SSTtoVIP, "Shared_SSTtoVIP_ORA")
perform_ORA(unique_pairs_VIPtoSST_p10, "Unique_VIPtoSST_p10_ORA")
perform_ORA(unique_pairs_VIPtoSST_p28, "Unique_VIPtoSST_p28_ORA")
perform_ORA(shared_pairs_VIPtoSST, "Shared_VIPtoSST_ORA")