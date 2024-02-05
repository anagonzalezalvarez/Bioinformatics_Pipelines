
library(DESeq2)
library(GEOquery)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(magrittr)
library(limma)
library(edgeR)
set.seed(123)

# 0 Read counts -----------------------------------------------------------------
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
dim(GSE198256_count)

# 0 Read Meta data -----------------------------------------------------------------
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Group <- Meta_GSE198256[,c("disease state:ch1")]
Group

## Make sure row names in samples match column names in counts
### all columns in counts in varaibles rows?
all(colnames(GSE198256_count) %in% rownames(Meta_GSE198256)) ##T
### same order?
all(colnames(GSE198256_count) == rownames(Meta_GSE198256)) ##T

# Rename my conditions just so its easier to understand
Group[Group=="Healthy"] <- "Controls"
Group[Group=="Covid19: Acute infection"] <- "Acute"
Group[Group=="Covid19: Recovery 3Mo"] <- "EarlyRecovery"
Group[Group=="Covid19: Recovery 6Mo"] <- "LateRecovery"

# 1 DE WORKFLOW ----------------------------------------------------------------------------

### LIMMA VOOM  ADD  CONTRASTS 
###  If Contrasts need
###  1. colnames(design)
###  2. contrast.matrix
###  3. contrasts.fit


## Create DGEList
dge <- DGEList(counts=GSE198256_count)
nrow(dge)
# Design
design <- model.matrix(~ Group )

# Filter
keep <- filterByExpr(dge, design=design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
nrow(dge)

# Normalization
dge <- calcNormFactors(dge)
dge$samples$norm.factors 

## VOOM
v <- voom(dge, design, plot=TRUE)
colnames(design) <- c("Intercept","Acute","EarlyRecovery","LateRecovery")

### intercept is Control is what is missing
### Controls Acute EarlyRecovery LateRecovery
fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(Acute, EarlyRecovery, 
                                 LateRecovery,    
                                 levels=design) ### Contrasts comparing everything to the intercept=Control
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
voom_topTable_CONTRASTS<-topTable(fit2)
#1Acute 2EarlyRecovery 3LateRecovery

results <- decideTests(fit2, adjust = "BH", p.value = 0.05, lfc = log2(2))
summary(decideTests(fit2, adjust = "BH", p.value = 0.05, lfc = log2(2)))
# Acute EarlyRecovery LateRecovery
# Down     629          2086          326
# NotSig 15448         13137        15681
# Up       332          1186          402

vennDiagram(results)
vennDiagram(results, include = "up")  # Only upregulated
vennDiagram(results, include = "down")  # Or downregulated


# 2 ORA WORKFLOW ----------------------------------------------------------------------------

# 1. Genes filtered
diff_table <- topTable(fit2,coef=3,p.value=0.01,number=10000)  ### CHANGE COEF TO CHECK THE OTHER CONTRASTS
genes_dif<- rownames(diff_table)

# Step 2: determine background.
background_set <- unique(rownames(v))

# Step 3: Determine gene sets.
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
hs_kegg_df <- subset(hs_msigdb_df, gs_cat == "C2" & gs_subcat == "CP:KEGG")


# Step 4: conduct ORA.
# 
selected_hs_kegg_df <- hs_kegg_df[, c("gs_name", "human_entrez_gene")]
kegg_ora_results <- enricher(
  gene = genes_dif,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  universe = background_set,
  TERM2GENE = selected_hs_kegg_df
)

# Step 5: Visualize / explore
enrich_plot <- enrichplot::dotplot(kegg_ora_results)
enrich_plot

upset_plot <- enrichplot::upsetplot(kegg_ora_results)
upset_plot

# 2 GSEA WORKFLOW ----------------------------------------------------------------------------
# msigdbr_collections()
# 1 C1     ""                         299
# 2 C2     "CGP"                     3384
# 3 C2     "CP"                        29
# 4 C2     "CP:BIOCARTA"              292
# 5 C2     "CP:KEGG"                  186
# 6 C2     "CP:PID"                   196
# 7 C2     "CP:REACTOME"             1615
# 8 C2     "CP:WIKIPATHWAYS"          664
# 9 C3     "MIR:MIRDB"               2377
# 10 C3     "MIR:MIR_Legacy"           221

# 1. Genes ALL p.value=1, number=nrow(logCPM)
diff_table_all <- topTable(fit2,coef=3,p.value=1,number=nrow(v)) 

# Step 2: determine background. -> NO BACKGROUND

# Step 3: Determine gene sets.
#hs_msigdb_df <- msigdbr(species = "Homo sapiens")
hs_kegg_df <- subset(hs_msigdb_df, gs_cat == "C2" & gs_subcat == "CP:REACTOME")

# Step 4: conduct GSEA
list_ordered <- diff_table_all[,"B"]
names(list_ordered) <- rownames(diff_table_all)

selected_hs_kegg_df <- hs_kegg_df[, c("gs_name", "human_entrez_gene")]
gsea_results <- GSEA(
  geneList = list_ordered,
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = selected_hs_kegg_df
)

# Step 5: Visualize / explore
enrich_plot <- enrichplot::dotplot(gsea_results)
enrich_plot

# Step 5: Visualize / explore
head(gsea_results@result)



