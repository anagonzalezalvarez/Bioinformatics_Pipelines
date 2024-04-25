# libraries 
#remotes::install_github('saezlab/liana', force = T)
#remotes::install_github("sqjin/CellChat")


library(CellChat)
library(liana)
library(igraph)
library(Seurat)
library(dplyr)
library(ggVennDiagram)
library(ggplot2)
#library(xlsx)

source("/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/bootstrapping_functions.R")

### SET WORKING SPACE -------------------------------------------------------------------------------------------------
basepath <- "/ibex/scratch/projects/c2169/Cecilia/CellToCell"

analysis_name <- "2_liana" 
path_imag <- paste0(basepath,"/code/",analysis_name,"/imag/")
path_dir <- paste0(basepath,"/code/",analysis_name,"/")

### Create it if it doesn't exist
dir.create(path_dir, showWarnings = FALSE)
dir.create(path_imag, showWarnings = FALSE)
setwd(path_imag)


## 1. Load annotated seurat object -------------------------------------------------------------------------------------------------
data <- readRDS("/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/1_scRNA/results/3_seurat_annotated.rds")
data@meta.data$seurat_annotations <- data@meta.data$predicted.subclass
head(data@meta.data)
p28 <- subset(data, sample == "P28")
table(p28$predicted.subclass)

## 2. Create separate objects for cell of interest -------------------------------------------------------------------------------------------------

Lamp5 <- subset(p28, predicted.subclass == "Lamp5")
Pvalb <- subset(p28, predicted.subclass == "Pvalb")
Vip <- subset(p28, predicted.subclass == "Vip")
Sst <- subset(p28, predicted.subclass == "Sst")


## 3. Run LIANA 100 times (bootstrapping) -------------------------------------------------------------------------------------------------
LIANA_bootstrap <- bootstrap_ccc(Lamp5, Pvalb, Vip, Sst, 500, 100, "predicted.subclass")
saveRDS(LIANA_bootstrap, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/LIANA_bootstrap_p28.rds")
LIANA_bootstrap <- readRDS(file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/LIANA_bootstrap_p28.rds")

## 4. Remove self interactions -------------------------------------------------------------------------------------------------
Lamp5_Sst_bootstrap_dir <- bootstrap_keep_dir(LIANA_bootstrap, "Lamp5", "Sst", 3)
Lamp5_Pvalb_bootstrap_dir <- bootstrap_keep_dir(LIANA_bootstrap, "Lamp5", "Pvalb", 3)
Lamp5_Vip_bootstrap_dir <- bootstrap_keep_dir(LIANA_bootstrap, "Lamp5", "Vip", 3)
Pvalb_Sst_bootstrap_dir <- bootstrap_keep_dir(LIANA_bootstrap, "Pvalb", "Sst", 3)
Pvalb_Vip_bootstrap_dir <- bootstrap_keep_dir(LIANA_bootstrap, "Pvalb", "Vip", 3)
Vip_Sst_bootstrap_dir <- bootstrap_keep_dir(LIANA_bootstrap, "Vip", "Sst", 3)

# Getting totals
list_of_lists <- list(Lamp5_Sst = Lamp5_Sst_bootstrap_dir,
  Lamp5_Pvalb = Lamp5_Pvalb_bootstrap_dir,
  Lamp5_Vip = Lamp5_Vip_bootstrap_dir,
  Pvalb_Sst = Pvalb_Sst_bootstrap_dir,
  Pvalb_Vip = Pvalb_Vip_bootstrap_dir,
  Vip_Sst = Vip_Sst_bootstrap_dir
)

total_rows_results <- sapply(list_of_lists, function(bootstrap_dir) {
  sum(sapply(bootstrap_dir, function(list_item) {
    sum(sapply(list_item, nrow))
  }))
})
print(total_rows_results)

## 5. Remove duplicates -------------------------------------------------------------------------------------------------

# DUPLICATED ERROR FIX: Loop through each list in the list of lists
remove_duplicates_from_lists <- function(data_structure) {
  # Loop through each list in the list of lists
  for (i in seq_along(data_structure)) {
    current_list <- data_structure[[i]]
    
    # Loop through each tibble in the current list
    for (j in names(current_list)) {
      current_tibble <- current_list[[j]]
      
      # Create a unique identifier for each row
      current_tibble$unique_id <- with(current_tibble, paste(source, target, ligand, receptor, sep = "⊎"))
      
      # Find all duplicates based on the unique identifier
      duplicate_indices <- duplicated(current_tibble$unique_id) | duplicated(current_tibble$unique_id, fromLast = TRUE)
      
      # Check if there are any duplicates and remove them
      if (any(duplicate_indices)) {        
        # Remove all duplicates by keeping rows where duplicate_indices is FALSE
        data_structure[[i]][[j]] <- current_tibble[!duplicate_indices, ]
      }
      
      # Remove the unique identifier column to clean up the tibble
      data_structure[[i]][[j]]$unique_id <- NULL
    }
  }
  
  return(data_structure)
}

Lamp5_Sst_bootstrap_dir_nodup <- remove_duplicates_from_lists(Lamp5_Sst_bootstrap_dir)
Lamp5_Pvalb_bootstrap_dir_nodup <- remove_duplicates_from_lists(Lamp5_Pvalb_bootstrap_dir)
Lamp5_Vip_bootstrap_dir_nodup <- remove_duplicates_from_lists(Lamp5_Vip_bootstrap_dir)
Pvalb_Sst_bootstrap_dir_nodup <- remove_duplicates_from_lists(Pvalb_Sst_bootstrap_dir)
Pvalb_Vip_bootstrap_dir_nodup <- remove_duplicates_from_lists(Pvalb_Vip_bootstrap_dir)
Vip_Sst_bootstrap_dir_nodup <- remove_duplicates_from_lists(Vip_Sst_bootstrap_dir)

# Getting totals
list_of_lists <- list(Lamp5_Sst = Lamp5_Sst_bootstrap_dir_nodup,
  Lamp5_Pvalb = Lamp5_Pvalb_bootstrap_dir_nodup,
  Lamp5_Vip = Lamp5_Vip_bootstrap_dir_nodup,
  Pvalb_Sst = Pvalb_Sst_bootstrap_dir_nodup,
  Pvalb_Vip = Pvalb_Vip_bootstrap_dir_nodup,
  Vip_Sst = Vip_Sst_bootstrap_dir
)

total_rows_results <- sapply(list_of_lists, function(bootstrap_dir) {
  sum(sapply(bootstrap_dir, function(list_item) {
    sum(sapply(list_item, nrow))
  }))
})
print(total_rows_results)

## 6. Aggregate -------------------------------------------------------------------------------------------------
Lamp5_Sst_bootstrap_agg <- aggregate_dir_inter(Lamp5_Sst_bootstrap_dir_nodup)
Lamp5_Pvalb_bootstrap_agg <- aggregate_dir_inter(Lamp5_Pvalb_bootstrap_dir_nodup)
Lamp5_Vip_bootstrap_agg <- aggregate_dir_inter(Lamp5_Vip_bootstrap_dir_nodup) ### manual modifiction
Pvalb_Sst_bootstrap_agg <- aggregate_dir_inter(Pvalb_Sst_bootstrap_dir_nodup)
Pvalb_Vip_bootstrap_agg <- aggregate_dir_inter(Pvalb_Vip_bootstrap_dir_nodup)
Vip_Sst_bootstrap_agg <- aggregate_dir_inter(Vip_Sst_bootstrap_dir_nodup)

## Visualize full object
Lamp5_Sst_bootstrap_agg[[1]] %>%
print(width=Inf)

# Getting totals
list_of_lists <- list(Lamp5_Sst = Lamp5_Sst_bootstrap_agg,
  Lamp5_Pvalb = Lamp5_Pvalb_bootstrap_agg,
  Lamp5_Vip = Lamp5_Vip_bootstrap_agg,
  Pvalb_Sst = Pvalb_Sst_bootstrap_agg,
  Pvalb_Vip = Pvalb_Vip_bootstrap_agg,
  Vip_Sst = Vip_Sst_bootstrap_agg
)
total_rows_by_group <- setNames(numeric(length(list_of_lists)), names(list_of_lists))
for (group in names(list_of_lists)) {
  total_rows <- 0
  for (i in seq_along(list_of_lists[[group]])) {
    total_rows <- total_rows + nrow(list_of_lists[[group]][[i]])
  }
  total_rows_by_group[group] <- total_rows
}
total_rows_by_group

## 7. Add interaction string T-S-L-R column -------------------------------------------------------------------------------------------------
Lamp5_Sst_bootstrap_agg <- add_inter_col(Lamp5_Sst_bootstrap_agg)
Lamp5_Pvalb_bootstrap_agg <- add_inter_col(Lamp5_Pvalb_bootstrap_agg)
Lamp5_Vip_bootstrap_agg <- add_inter_col(Lamp5_Vip_bootstrap_agg) ### NEED TO CHANGE
Pvalb_Sst_bootstrap_agg <- add_inter_col(Pvalb_Sst_bootstrap_agg)
Pvalb_Vip_bootstrap_agg <- add_inter_col(Pvalb_Vip_bootstrap_agg)
Vip_Sst_bootstrap_agg <- add_inter_col(Vip_Sst_bootstrap_agg)

## Visualize full object
Lamp5_Sst_bootstrap_agg[[1]] %>%
print(width=Inf)

# SAVING 
saveRDS(Lamp5_Sst_bootstrap_agg, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/Lamp5_Sst_bootstrap_agg_28.rds")
saveRDS(Lamp5_Pvalb_bootstrap_agg, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/Lamp5_Pvalb_bootstrap_agg_28.rds")
saveRDS(Lamp5_Vip_bootstrap_agg, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/Lamp5_Vip_bootstrap_agg_28.rds") ### NEED TO CHANGE
saveRDS(Pvalb_Sst_bootstrap_agg, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/Pvalb_Sst_bootstrap_agg_28.rds")
saveRDS(Pvalb_Vip_bootstrap_agg, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/Pvalb_Vip_bootstrap_agg_28.rds")
saveRDS(Vip_Sst_bootstrap_agg, file = "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/Vip_Sst_bootstrap_agg_28.rds")


## 8. Filter `aggregate_rank < 0.01`  -------------------------------------------------------------------------------------------------
Lamp5_Sst_bootstrap_agg_filt <- filter_aggRank(Lamp5_Sst_bootstrap_agg, 0.01)
Lamp5_Pvalb_bootstrap_agg_filt <- filter_aggRank(Lamp5_Pvalb_bootstrap_agg, 0.01)
Lamp5_Vip_bootstrap_agg_filt <- filter_aggRank(Lamp5_Vip_bootstrap_agg, 0.01)
Pvalb_Sst_bootstrap_agg_filt <- filter_aggRank(Pvalb_Sst_bootstrap_agg, 0.01)
Pvalb_Vip_bootstrap_agg_filt <- filter_aggRank(Pvalb_Vip_bootstrap_agg, 0.01)
Vip_Sst_bootstrap_agg_filt <- filter_aggRank(Vip_Sst_bootstrap_agg, 0.01)


# Getting totals 
list_of_lists <- list(Lamp5_Sst = Lamp5_Sst_bootstrap_agg_filt,
  Lamp5_Pvalb = Lamp5_Pvalb_bootstrap_agg_filt,
  Lamp5_Vip = Lamp5_Vip_bootstrap_agg_filt,
  Pvalb_Sst = Pvalb_Sst_bootstrap_agg_filt,
  Pvalb_Vip = Pvalb_Vip_bootstrap_agg_filt,
  Vip_Sst = Vip_Sst_bootstrap_agg_filt
)
total_rows_by_group <- setNames(numeric(length(list_of_lists)), names(list_of_lists))
for (group in names(list_of_lists)) {
  total_rows <- 0
  for (i in seq_along(list_of_lists[[group]])) {
    total_rows <- total_rows + nrow(list_of_lists[[group]][[i]])
  }
  total_rows_by_group[group] <- total_rows
}
total_rows_by_group


## 9. Get vector with all possible interactions predicted in the 100 runs and get the percentage of presence in the bootstrapped results -------------------------------------------------------------------------------------------------
Lamp5_Sst_interactions <- interactions_vector(Lamp5_Sst_bootstrap_agg_filt) # all interactions possible 
Lamp5_Pvalb_interactions <- interactions_vector(Lamp5_Pvalb_bootstrap_agg_filt) # all interactions possible 
Lamp5_Vip_interactions <- interactions_vector(Lamp5_Vip_bootstrap_agg_filt) # all interactions possible 
Pvalb_Sst_interactions <- interactions_vector(Pvalb_Sst_bootstrap_agg_filt) # all interactions possible 
Pvalb_Vip_interactions <- interactions_vector(Pvalb_Vip_bootstrap_agg_filt) # all interactions possible 
Vip_Sst_interactions <- interactions_vector(Vip_Sst_bootstrap_agg_filt) # all interactions possible 


robust_Lamp5_Sst <- check_inter(Lamp5_Sst_bootstrap_agg_filt, Lamp5_Sst_interactions)
robust_Lamp5_Pvalb <- check_inter(Lamp5_Pvalb_bootstrap_agg_filt, Lamp5_Pvalb_interactions)
robust_Lamp5_Vip <- check_inter(Lamp5_Vip_bootstrap_agg_filt, Lamp5_Vip_interactions)
robust_Pvalb_Sst <- check_inter(Pvalb_Sst_bootstrap_agg_filt, Pvalb_Sst_interactions)
robust_Pvalb_Vip <- check_inter(Pvalb_Vip_bootstrap_agg_filt, Pvalb_Vip_interactions)
robust_Vip_Sst <- check_inter(Vip_Sst_bootstrap_agg_filt, Vip_Sst_interactions)

## Separate for interactions c1 -> c2
robust_Lamp5toSst <- robust_Lamp5_Sst[startsWith(robust_Lamp5_Sst$interactions, "L"), ] 
robust_Lamp5toPvalb <- robust_Lamp5_Pvalb[startsWith(robust_Lamp5_Pvalb$interactions, "L"), ] 
robust_Lamp5toVip <- robust_Lamp5_Vip[startsWith(robust_Lamp5_Vip$interactions, "L"), ] 
robust_PvalbtoSst <- robust_Pvalb_Sst[startsWith(robust_Pvalb_Sst$interactions, "P"), ] 
robust_PvalbtoVip <- robust_Pvalb_Vip[startsWith(robust_Pvalb_Vip$interactions, "P"), ] 
robust_ViptoSst <- robust_Vip_Sst[startsWith(robust_Vip_Sst$interactions, "V"), ] 

## Separate for interactions c2 -> c1
robust_SsttoLamp5 <- robust_Lamp5_Sst[startsWith(robust_Lamp5_Sst$interactions, "S"), ] 
robust_PvalbtoLamp5 <- robust_Lamp5_Pvalb[startsWith(robust_Lamp5_Pvalb$interactions, "P"), ] 
robust_ViptoLamp5 <- robust_Lamp5_Vip[startsWith(robust_Lamp5_Vip$interactions, "V"), ] 
robust_SsttoPvalb <- robust_Pvalb_Sst[startsWith(robust_Pvalb_Sst$interactions, "S"), ] 
robust_ViptoPvalb <- robust_Pvalb_Vip[startsWith(robust_Pvalb_Vip$interactions, "V"), ] 
robust_SsttoVip <- robust_Vip_Sst[startsWith(robust_Vip_Sst$interactions, "S"), ] 


## 10 Plotting -------------------------------------------------------------------------------------------------

robust_Lamp5_Sst$dir <- "Lamp5-Sst"
robust_Lamp5_Pvalb$dir <- "Lamp5_Pvalb"
robust_Lamp5_Vip$dir <- "Lamp5_Vip"
robust_Pvalb_Sst$dir <- "Pvalb_Sst"
robust_Pvalb_Vip$dir <- "Pvalb_Vip"
robust_Vip_Sst$dir <- "Vip_Sst"

all <- rbind(robust_Lamp5_Sst,
             robust_Lamp5_Pvalb,
             robust_Lamp5_Vip,
             robust_Pvalb_Sst,
             robust_Pvalb_Vip,
             robust_Vip_Sst)

# ploting and checking desnity
all %>%
  ggplot(aes(x=percentage)) +
  geom_density(fill="indianred", color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept = 86.30, linetype = "dashed", color = "red") +
  theme_minimal()  # You can customize the theme as needed

pdf("1_all_interactions_density_p28.pdf")
# ploting and checking desnity
all %>%
  ggplot(aes(x=percentage)) +
  geom_density(fill="indianred", color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept = 86.30, linetype = "dashed", color = "red") +
  theme_minimal()  # You can customize the theme as needed

all 
ggplot(all, aes(x = interactions, y = percentage, fill = dir)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size =3.5))  # You can customize the theme as needed
dev.off()


## 11 Filtering for only interactions withing trshold  -------------------------------------------------------------------------------------------------

### 11.1 Keep only interactions above 75% of the 100 predictions

##interactions c1 -> c2
robust_Lamp5toSst_vecfilt <- subset(robust_Lamp5toSst, subset = robust_Lamp5toSst$percentage >= 75)
robust_Lamp5toPvalb_vecfilt <- subset(robust_Lamp5toPvalb, subset = robust_Lamp5toPvalb$percentage >= 75)
robust_Lamp5toVip_vecfilt <- subset(robust_Lamp5toVip, subset = robust_Lamp5toVip$percentage >= 75)
robust_PvalbtoSst_vecfilt <- subset(robust_PvalbtoSst, subset = robust_PvalbtoSst$percentage >= 75)
robust_PvalbtoVip_vecfilt <- subset(robust_PvalbtoVip, subset = robust_PvalbtoVip$percentage >= 75)
robust_ViptoSst_vecfilt <- subset(robust_ViptoSst, subset = robust_ViptoSst$percentage >= 75)

##interactions c2 -> c1
robust_SsttoLamp5_vecfilt <- subset(robust_SsttoLamp5, subset = robust_SsttoLamp5$percentage >= 75)
robust_PvalbtoLamp5_vecfilt <- subset(robust_PvalbtoLamp5, subset = robust_PvalbtoLamp5$percentage >= 75)
robust_ViptoLamp5_vecfilt <- subset(robust_ViptoLamp5, subset = robust_ViptoLamp5$percentage >= 75)
robust_SsttoPvalb_vecfilt <- subset(robust_SsttoPvalb, subset = robust_SsttoPvalb$percentage >= 75)
robust_ViptoPvalb_vecfilt <- subset(robust_ViptoPvalb, subset = robust_ViptoPvalb$percentage >= 75)
robust_SsttoVip_vecfilt <- subset(robust_SsttoVip, subset = robust_SsttoVip$percentage >= 75)

### 11.2 Subset the bootstrapped results (still a list of 100 items)

##interactions c1 -> c2
robust_Lamp5toSst_agg <- keep_specific_inter(Lamp5_Sst_bootstrap_agg_filt, robust_Lamp5toSst_vecfilt$interactions) 
robust_Lamp5toPvalb_agg <- keep_specific_inter(Lamp5_Pvalb_bootstrap_agg_filt, robust_Lamp5toPvalb_vecfilt$interactions) 
robust_Lamp5toVip_agg <- keep_specific_inter(Lamp5_Vip_bootstrap_agg_filt, robust_Lamp5toVip_vecfilt$interactions) 
robust_PvalbtoSst_agg <- keep_specific_inter(Pvalb_Sst_bootstrap_agg_filt, robust_PvalbtoSst_vecfilt$interactions) 
robust_PvalbtoVip_agg <- keep_specific_inter(Pvalb_Vip_bootstrap_agg_filt, robust_PvalbtoVip_vecfilt$interactions) 
robust_ViptoSst_agg <- keep_specific_inter(Vip_Sst_bootstrap_agg_filt, robust_ViptoSst_vecfilt$interactions) 

##interactions c2 -> c1
robust_SsttoLamp5_agg <- keep_specific_inter(Lamp5_Sst_bootstrap_agg_filt, robust_SsttoLamp5_vecfilt$interactions) 
robust_PvalbtoLamp5_agg <- keep_specific_inter(Lamp5_Pvalb_bootstrap_agg_filt, robust_PvalbtoLamp5_vecfilt$interactions) 
robust_ViptoLamp5_agg <- keep_specific_inter(Lamp5_Vip_bootstrap_agg_filt, robust_ViptoLamp5_vecfilt$interactions) 
robust_SsttoPvalb_agg <- keep_specific_inter(Pvalb_Sst_bootstrap_agg_filt, robust_SsttoPvalb_vecfilt$interactions) 
robust_ViptoPvalb_agg <- keep_specific_inter(Pvalb_Vip_bootstrap_agg_filt, robust_ViptoPvalb_vecfilt$interactions) 
robust_SsttoVip_agg <- keep_specific_inter(Vip_Sst_bootstrap_agg_filt, robust_SsttoVip_vecfilt$interactions) 


### 11.3 Randomly choose one set of interactions (a single df) from the bootstrapped resutls 
robust_Lamp5toSst_agg_df <- robust_Lamp5toSst_agg[[sample(1:100, 1)]]
robust_Lamp5toPvalb_agg_df <- robust_Lamp5toPvalb_agg[[sample(1:100, 1)]]
robust_Lamp5toVip_agg_df <- robust_Lamp5toVip_agg[[sample(1:100, 1)]]
robust_PvalbtoSst_agg_df <- robust_PvalbtoSst_agg[[sample(1:100, 1)]]
robust_PvalbtoVip_agg_df <- robust_PvalbtoVip_agg[[sample(1:100, 1)]]
robust_ViptoSst_agg_df <- robust_ViptoSst_agg[[sample(1:100, 1)]]

robust_SsttoLamp5_agg_df <- robust_SsttoLamp5_agg[[sample(1:100, 1)]]
robust_PvalbtoLamp5_agg_df <- robust_PvalbtoLamp5_agg[[sample(1:100, 1)]]
robust_ViptoLamp5_agg_df <- robust_ViptoLamp5_agg[[sample(1:100, 1)]]
robust_SsttoPvalb_agg_df <- robust_SsttoPvalb_agg[[sample(1:100, 1)]]
robust_ViptoPvalb_agg_df <- robust_ViptoPvalb_agg[[sample(1:100, 1)]]
robust_SsttoVip_agg_df <- robust_SsttoVip_agg[[sample(1:100, 1)]]

getwd()






## 12 Plotting  -------------------------------------------------------------------------------------------------
pdf("2_heatmaps_p28.pdf")
robust_all <- rbind(robust_Lamp5toSst_agg_df,robust_Lamp5toPvalb_agg_df,robust_Lamp5toVip_agg_df,robust_PvalbtoSst_agg_df,robust_PvalbtoVip_agg_df,robust_ViptoSst_agg_df,robust_SsttoLamp5_agg_df,robust_PvalbtoLamp5_agg_df,robust_ViptoLamp5_agg_df,robust_SsttoPvalb_agg_df,robust_ViptoPvalb_agg_df,robust_SsttoVip_agg_df)
heat_freq(robust_all)

robust_all_Lamp5 <- rbind(robust_Lamp5toSst_agg_df, 
                         robust_Lamp5toPvalb_agg_df,
                         #robust_Lamp5toVip_agg_df,
                         robust_SsttoLamp5_agg_df,
                         robust_PvalbtoLamp5_agg_df
                         #,robust_ViptoLamp5_agg_df
                         )
heat_freq(robust_all_Lamp5)

robust_all_Sst <- rbind(robust_Lamp5toSst_agg_df, 
                         robust_PvalbtoSst_agg_df,
                         robust_ViptoSst_agg_df,
                         robust_SsttoLamp5_agg_df,
                         robust_SsttoPvalb_agg_df,
                         robust_SsttoVip_agg_df
                         )
heat_freq(robust_all_Sst)

robust_all_Pvalb <- rbind(robust_Lamp5toPvalb_agg_df, 
                         robust_PvalbtoSst_agg_df,
                         robust_PvalbtoVip_agg_df,
                         robust_PvalbtoLamp5_agg_df,
                         robust_SsttoPvalb_agg_df,
                         robust_ViptoPvalb_agg_df
                         )
heat_freq(robust_all_Pvalb)

robust_all_Vip <- rbind(#robust_Lamp5toVip_agg_df, 
                         robust_PvalbtoVip_agg_df,
                         robust_ViptoSst_agg_df,
                         #robust_ViptoLamp5_agg_df,
                         robust_ViptoPvalb_agg_df,
                         robust_SsttoVip_agg_df
                         )
heat_freq(robust_all_Vip)
dev.off()


### Venn Diagrams

# VennDiagram plots 

### Save into df
library(dplyr)

robust_all_df <- bind_rows(
  robust_Lamp5toSst_agg_df,
  robust_Lamp5toPvalb_agg_df,
  robust_Lamp5toVip_agg_df,  # Omitted as it's commented out in your setup
  robust_PvalbtoSst_agg_df,
  robust_PvalbtoVip_agg_df,
  robust_ViptoSst_agg_df,
  robust_SsttoLamp5_agg_df,
  robust_PvalbtoLamp5_agg_df,
  robust_ViptoLamp5_agg_df,  # Omitted as it's commented out in your setup
  robust_SsttoPvalb_agg_df,
  robust_ViptoPvalb_agg_df,
  robust_SsttoVip_agg_df
)

## save as csv
write.csv(robust_all_df, "/ibex/scratch/projects/c2169/Cecilia/CellToCell/code/2_liana/results/robust_all_df_p28.csv")



pdf("2_totalinteractions_p28.pdf")

data_summary <- robust_all_df %>%
  group_by(source, target) %>%
  summarise(count = n(), .groups = 'drop')

ggplot(data_summary, aes(x = interaction(paste(source, target, sep="-")), y = count, fill = source)) +
  geom_bar(stat = "identity") +  # Use identity to use the counts as heights
  labs(title = "Count of Rows Per Source-Target Combination",
       x = "Source-Target Combination",
       y = "Count of Rows") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability



dev.off()


pdf("3_venn.pdf")
PvalbtoSst_vectors <- paste(robust_PvalbtoSst_agg_df$ligand, robust_PvalbtoSst_agg_df$receptor, sep = "-")
SsttoPvalb_vectors <- paste(robust_SsttoPvalb_agg_df$ligand, robust_SsttoPvalb_agg_df$receptor, sep = "-")

x <- list(PtoS = PvalbtoSst_vectors,
          StoP = SsttoPvalb_vectors)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom")

Lamp5toSst_vectors <- paste(robust_Lamp5toSst_agg_df$ligand, robust_Lamp5toSst_agg_df$receptor, sep = "-")
SsttoLamp5_vectors <- paste(robust_SsttoLamp5_agg_df$ligand, robust_SsttoLamp5_agg_df$receptor, sep = "-")

x <- list(LtoS = Lamp5toSst_vectors,
          StoL = SsttoLamp5_vectors)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom")


dev.off()
