
# Built-in Functions

#### 1) Sub-sampling objects and store in a list ####

## Function to sample 
SampleCells <- function(seurat_object, sample_size){
  set.seed(NULL)
  supsampled_seurat_object <- seurat_object[, sample(colnames(seurat_object), size = sample_size, replace=F)]
  return(supsampled_seurat_object)
}

## Function to combine three sub-sampled datasets
combine_subsampled <- function(se1, se2, se3, se4){
  merged_se <-  merge(se1, y=c(se2, se3, se4))
  return(merged_se)
}


## Main bootstrapping
bootstrap_ccc <- function(seurat_object1, seurat_object2, seurat_object3, seurat_object4, sample_size, N, ident){
  num_bootstraps <- N   # Set the number of bootstraps
  bootstrap_results <- list()
  subsampled_seurat_1_list <- lapply(1:num_bootstraps, function(i) {SampleCells(seurat_object1, sample_size)})
  subsampled_seurat_2_list <- lapply(1:num_bootstraps, function(i) {SampleCells(seurat_object2, sample_size)})
  subsampled_seurat_3_list <- lapply(1:num_bootstraps, function(i) {SampleCells(seurat_object3, sample_size)})
  subsampled_seurat_4_list <- lapply(1:num_bootstraps, function(i) {SampleCells(seurat_object4, sample_size)})

  
  # perform bootstrap
  for (i in 1:num_bootstraps){
    subsampled_se1 <- subsampled_seurat_1_list[[i]]
    subsampled_se2 <- subsampled_seurat_2_list[[i]]
    subsampled_se3 <- subsampled_seurat_3_list[[i]]
    subsampled_se4 <- subsampled_seurat_4_list[[i]]
    merged_se <- combine_subsampled(subsampled_se1, subsampled_se2, subsampled_se3, subsampled_se4)
    merged_se[["RNA"]] <- JoinLayers(merged_se[["RNA"]])

    liana_result <- liana_wrap(merged_se, method = c('call_cellchat', 'cellphonedb', 'sca'),
                               resource = 'MouseConsensus', idents_col = ident, return_all = T)
    # Store the results in the list
    bootstrap_results[[i]] <- liana_result
  } 

  # Store the results in the list
  return(bootstrap_results)
}


#### 2) Function for keeping direct interaction #### 
bootstrap_results_dir <- list() # Initialize a list to store the results
bootstrap_keep_dir <- function(list, celltype1, celltype2, n_methods){
  for (n in 1:length(list)){
    liana_dir <- list()
    for (x in 1:n_methods){
      
      l <- list[[n]][[x]][which(list[[n]][[x]]$source == celltype1 & list[[n]][[x]]$target == celltype2 | 
                                  list[[n]][[x]]$source == celltype2 & list[[n]][[x]]$target == celltype1), ]
      liana_dir[[x]] <- l
    }
    names(liana_dir) <- c('call_cellchat', 'cellphonedb', 'sca')
    bootstrap_results_dir[[n]] <- liana_dir
  }
  return(bootstrap_results_dir)
}


#### 3) Function to aggregate results #### 
bootstrap_results_agg <- list() 
aggregate_dir_inter <- function(list){
  for (x in 1:length(list)){
    for (df in names(list[[x]])) {
      row.names(list[[x]][[df]]) <- NULL  # Remove existing row names
    }
    a <- liana_aggregate(list[[x]]) 
    bootstrap_results_agg[[x]] <- a
  }
  return(bootstrap_results_agg)
}


#### 4) Function to add metadata with interaction #### 
new_list <- list() 
add_inter_col <- function(list){
  for (x in 1:length(list)){
    vect <- paste(list[[x]]$source, list[[x]]$target, list[[x]]$ligand, list[[x]]$receptor, sep = "-")
    list[[x]]$interaction <- vect
    new_list[[x]] <- list[[x]]
  }
  return(new_list)
}


#### 5) Function to filter aggregate_rank #### 
fitered_list <- list() 
filter_aggRank <- function(list, pval_threshold){
  for (x in 1:length(list)){
    f <- list[[x]] %>% filter(aggregate_rank < pval_threshold)
    fitered_list[[x]] <- f
  }
  return(fitered_list)
}

#### 6) generate a vector of all interactions predicted #### 
all_interactions <- c()
interactions_vector <- function(list){
  for (x in 1:length(list)){
    inter <- list[[x]]$interaction
    all_interactions <- c(all_interactions, inter)
  }
  all_interactions <- unique(all_interactions) # to remove duplicates 
  return(all_interactions)
}


#### 7) To check for how man times(%) an interaction appears in the bootstrapped results #### 
interactions <- c()
percentage <- c()
check_inter <- function(list, inter_vec){
  for (n in 1:length(inter_vec)){
    total_count <- 0
    for (i in 1:length(list)){
      count <- sum(list[[i]]$interaction == inter_vec[n])
      total_count <- total_count + count
    }
    percentage <- c(percentage, total_count)
    interactions <- c(interactions, inter_vec[n])
  }
  df <- data.frame(interactions, percentage)
  return(df)
}

#### 8) Keeping interactions that are of interest #### 
new_list <- list() 
keep_specific_inter <- function(list, interactions_vec){
  for (x in 1:length(list)){
    all_inter <- list[[x]]$interaction
    index <- as.numeric(which(all_inter %in% interactions_vec))
    l <- list[[x]][c(index),]
    new_list[[x]] <- l
  }
  return(new_list)
}


