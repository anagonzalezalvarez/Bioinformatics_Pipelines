### Packages ----------------------------------------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")
# BiocManager::install("NOISeq")


library("NOISeq")
library("tidyverse")
library("biomaRt")
library(DESeq2)


setwd("~/Desktop/PhD_2nd/BESE398_Pipelines/Bioinformatics_Pipelines/W1/1_ClassExcercise")

### Load Data ----------------------------------------------------------------------------------

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
ncol(GSE198256_count)
nrow(GSE198256_count) ##39376


### Metadata ----------------------------------------------------------------------------------

# load metadata table from GEO
library(GEOquery)
gds <- getGEO("GSE198256")
## https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
## Retrieve data
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)

## Filter for only the columns i am interested in
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]

Factors_GSE198256 <- Meta_GSE198256[,c("disease state:ch1")]



### NOISESeq ----------------------------------------------------------------------------------
# We need length, gc, biotype, chromosome

##### BIOMART

### Get rownames which contain entrez id as csv
entrez_id <- rownames(GSE198256_count)

# Save to  file
write.table(entrez_id, file="data/biomart/entrez_id.csv",col.names = FALSE,row.names = FALSE,quote=F)

#### Go to biomart and upload the file to retrive length, gc, biotype, chromosome

annotgene <- read.table("data/biomart/gene_info_mart_export.txt",  sep = "\t", header = TRUE)
colnames(annotgene)

#### Rename cols
colnames(annotgene)[colnames(annotgene) == 'NCBI.gene..formerly.Entrezgene..ID'] <- 'Entrezgene'
colnames(annotgene)[colnames(annotgene) == 'Chromosome.scaffold.name'] <- 'Chromosome'
colnames(annotgene)[colnames(annotgene) == 'Gene.start..bp.'] <- 'start'
colnames(annotgene)[colnames(annotgene) == 'Gene.end..bp.'] <- 'end'
colnames(annotgene)[colnames(annotgene) == 'Gene...GC.content'] <- 'GC'
colnames(annotgene)[colnames(annotgene) == 'Transcript.type'] <- 'type'

## create matrix for summary of filtered counts
summary <- matrix(NA, nrow = 4, ncol = 4)
colnames(summary) <- c("filter","n_rows", "n_unique", "n_in_GSE198256_count")

summary[1,1] <- "original"
summary[1,2] <- nrow(annotgene)
summary[1,3] <- length(unique(annotgene$Entrezgene))
summary[1,4] <- sum(rownames(GSE198256_count) %in% annotgene$Entrezgene) # How many genes do I get annotated?


# Filter the information
annotgene <- annotgene[annotgene$Chromosome %in% c(as.character(1:22) ,"X","Y"),]

summary[2,1] <- "chromosome filter"
summary[2,2] <- nrow(annotgene)
summary[2,3] <- length(unique(annotgene$Entrezgene))
summary[2,4] <- sum(rownames(GSE198256_count) %in% annotgene$Entrezgene) # How many genes do I get annotated?


## Multiple entries for one entrezid -> need to keep only 1 to be able to set it as rownames
annotgene_filt <- annotgene[!duplicated(annotgene$Entrezgene),] ### remove duplicates

summary[3,1] <- "duplicate removal"
summary[3,2] <- nrow(annotgene_filt)
summary[3,3] <- length(unique(annotgene_filt$Entrezgene))
summary[3,4] <- sum(rownames(GSE198256_count) %in% annotgene_filt$Entrezgene) # How many genes do I get annotated?


## Remove nas
annotgene_filt <- na.omit(annotgene_filt) ##27,234 remove na's in annotgene_filt$Entrezgene

summary[4,1] <- "nas removal"
summary[4,2] <- nrow(annotgene_filt)
summary[4,3] <- length(unique(annotgene_filt$Entrezgene))
summary[4,4] <- sum(rownames(GSE198256_count) %in% annotgene_filt$Entrezgene) # How many genes do I get annotated?


## Match annotation annotgene_filt$Entrezgene with counts GSE198256_count
rownames(annotgene_filt) <- as.character(annotgene_filt$Entrezgene) ### set annotgene_filt$Entrezgene as rowname 

##  Filter GSE198256_count for only the ones present in annotgene_filt (the ones i retrieved the biomart info)
GSE198256_count_filt <- GSE198256_count[rownames(GSE198256_count) %in% rownames(annotgene_filt),]
nrow(GSE198256_count_filt)
## Also save the excluded ones :)
GSE198256_count_exc <-GSE198256_count[!(rownames(GSE198256_count) %in% rownames(annotgene_filt)),]
## reorder annotgene_filt to the order in GSE198256_count_filt
annotgene_ord <- annotgene_filt[rownames(GSE198256_count_filt ),]
## now both annotgene_ord and GSE198256_count_filt contain same entrez id in the same order

### NOISEQ WORKFLOW ---------------------------------------------------------------------------

### Get factors order of the condition col
Factors_GSE198256 <- data.frame(Meta_GSE198256 [ colnames(GSE198256_count_filt),c("disease state:ch1")])
colnames(Factors_GSE198256)[1]<- "Group" ### RENAME 

### Extract column for lenght, GC, biotype and chromosome and create a vector with the entrez id's as names

lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
names(lengthuse) <- rownames(annotgene_ord)
head(lengthuse)

gc <- annotgene_ord$GC ##exctract col
names(gc) <- rownames(annotgene_ord) ## add names

biotype <-annotgene_ord$type ##exctract col
names(biotype) <- rownames(annotgene_ord) ## add names

chromosome <- annotgene_ord[,c("Chromosome","start","end")]

### RUN NOISEQ --------------------------------------------------------------
data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)

myexplodata <- dat(data_NOISEQ, type = "biodetection")

### EXPLORE RESULTS --------------------------------------------------------------

explo.plot(myexplodata, plottype = "persample")

par(mfrow = c(1, 2))
explo.plot(myexplodata, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")


mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")

mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:2, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 30:34)

explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")

explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")

mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd,samples = 1:12)

myPCA = dat(data_NOISEQ, type = "PCA")
explo.plot(myPCA, factor = "Group")

QCreport(data_NOISEQ, samples = NULL, factor = "Group", norm = FALSE)

save(data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="GSE198256_step1.Rda")

### NORMALIZATION  --------------------------------------------------------------

myRPKM = rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA = uqua(assayData(data_NOISEQ)$exprs, long = lengthuse, lc = 0.5, k = 0)
myTMM = tmm(assayData(data_NOISEQ)$exprs, long = 1000, lc = 0)

### DESEQ --------------------------------------------------------------
## Need:
## 1. Counts: GSE198256_count_filt
## 2. Sample info: data_NOISEQ@phenoData@data

# 1 prepare data ----------------------------------------------------------------------------
## Make sure row names in samples match column names in counts
### all columns in counts in varaibles rows?
all(colnames(GSE198256_count_filt) %in% rownames(data_NOISEQ@phenoData@data)) ##T
### same order?
all(colnames(GSE198256_count_filt) == rownames(data_NOISEQ@phenoData@data)) ##T

# Rename my conditions just so its easier to understand
pDataUSE <- data_NOISEQ@phenoData@data
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="Healthy"] <- "Controls"
pDataUSE[pDataUSE=="Covid19: Acute infection"] <- "Acute"
pDataUSE[pDataUSE=="Covid19: Recovery 3Mo"] <- "EarlyRecovery"
pDataUSE[pDataUSE=="Covid19: Recovery 6Mo"] <- "LateRecovery"

#pDataUSE[pDataUSE=="Healthy"] <- "Covid19AI"
#pDataUSE[pDataUSE=="Covid19: Acute infection"] <- "Covid19AI"
#pDataUSE[pDataUSE=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
#pDataUSE[pDataUSE=="Covid19: Recovery 6Mo"] <- "Covid196Mo"

pDataUSE[,1] <- as.factor(pDataUSE[,1]) ## as factors


# 2 Create DESeqDataSet object --------------------------------------------------------------------------
GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ Group) ## in colData what is the name of the col that contains condition


### In case we do not want to use intercept in the design use: design = ~ -1 + Group

# 3 QC and Filtering  ---------------------------------------------------------------------------
# Keeping only rows with at least 10 reads
smallestGroupSize <- 6
keep <- rowSums(counts(GSE198256_DESeq2) >= 10) >= smallestGroupSize
GSE198256_DESeq2_F <- GSE198256_DESeq2[keep,]
nrow(GSE198256_DESeq2)
nrow(GSE198256_DESeq2_F) ###14,055

# 4 Factors  ---------------------------------------------------------------------------
GSE198256_DESeq2_F$Group ## Levels: Covid19: Acute infection Covid19: Recovery 3Mo Covid19: Recovery 6Mo Healthy
GSE198256_DESeq2_F$Group <-  relevel(GSE198256_DESeq2_F$Group, ref = 'Controls') ## set healthy as ref

# 5 Run DESeq ---------------------------------------------------------------------------
GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
resultsNames(GSE198256_DESeq2_F) ##[1] "Intercept" "Group_Acute_vs_Controls"   "Group_EarlyRecovery_vs_Controls" "Group_LateRecovery_vs_Controls" 

GSE198256_res <- results(GSE198256_DESeq2_F)
GSE198256_res

# 6 Explore Results ---------------------------------------------------------------------------
summary(results(GSE198256_DESeq2_F, contrast =c('Group','Acute','Controls')))
#LFC > 0 (up)       : 1093, 7.8%
#LFC < 0 (down)     : 1442, 10%

summary(results(GSE198256_DESeq2_F, contrast =c('Group','EarlyRecovery','Controls')))
##LFC > 0 (up)       : 1184, 8.4%
##LFC < 0 (down)     : 1246, 8.9%
summary(results(GSE198256_DESeq2_F, contrast =c('Group','LateRecovery','Controls')))
#LFC > 0 (up)       : 98, 0.7%
#LFC < 0 (down)     : 40, 0.28%


### 6.1 Contrasts ---------------------------------------------------------------------------
res <- results(GSE198256_DESeq2_F, contrast=c('Group','Acute','Controls'))
filtered_res <- res[!is.na(res$log2FoldChange) & 
                      !is.na(res$pvalue) &                            
                      (abs(res$log2FoldChange) > 1) &                            
                      res$padj < 0.05, ]
summary(filtered_res)
##LFC > 0 (up)       : 369, 37%
##LFC < 0 (down)     : 639, 63%

res <- results(GSE198256_DESeq2_F, contrast=c('Group','EarlyRecovery','Controls'))
filtered_res <- res[!is.na(res$log2FoldChange) & 
                      !is.na(res$pvalue) &                            
                      (abs(res$log2FoldChange) > 1) &                            
                      res$padj < 0.05, ]
summary(filtered_res)
##LFC > 0 (up)       : 419, 51%
##LFC < 0 (down)     : 396, 49%

res <- results(GSE198256_DESeq2_F, contrast=c('Group','LateRecovery','Controls'))
filtered_res <- res[!is.na(res$log2FoldChange) & 
                      !is.na(res$pvalue) &  !is.na(res$padj)   &                       
                      (abs(res$log2FoldChange) > 1) &                            
                      res$padj < 0.05, ]
summary(filtered_res)
##LFC > 0 (up)       : 20, 83%
##LFC < 0 (down)     : 4, 17%


### 7. Shrinkage ---------------------------------------------------------------------------
res_lfcShrink <- lfcShrink(GSE198256_DESeq2_F,coef=c("Group_EarlyRecovery_vs_Controls"))
res_lfcShrink <- lfcShrink(GSE198256_DESeq2_F,coef=c("Group_Covid196Mo_vs_Covid193Mo"))
plotMA(res_lfcShrink, ylim=c(-2,2))


#GSE198256_DESeq2_F <- DESeq(GSE198256_DESeq2_F, test="LRT", reduced=~1)
#GSE198256_DESeq2_res_LRT <- results(GSE198256_DESeq2_F)
#GSE198256_DESeq2_res_LRT
#res <- results(GSE198256_DESeq2_res_LRT)

# Technical replicates?

# How to interpret the results?

plotCounts(GSE198256_DESeq2_F, gene="100287102", intgroup="Group")

## STEP 3.1.4: QC??

# How do visualize?
vsd <- vst(GSE198256_DESeq2_F, blind=FALSE)
rld <- rlog(GSE198256_DESeq2_F, blind=FALSE)
head(assay(vsd), 3)

# heatmap
library("pheatmap")
select <- order(rowMeans(counts(GSE198256_DESeq2_F,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(GSE198256_DESeq2_F)[,c("Group")])
colnames(df) <- "Group"

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE,annotation_col=df)

# PCA
plotPCA(vsd, intgroup=c("Group"))






############
# STEP 3. BIS: NORMALIZATION AND DIFFERENTIAL EXPRESSION BASED ON LIMMA
############

# Starts at Page 70
# https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

# In the limma approach to RNA-seq, read counts are converted to log2-counts-per-million
# (logCPM) and the mean-variance relationship is modeled either with precision weights
# or with an empirical Bayes prior trend. The precision weights approach is called 
# “voom” and the prior trend approach is called “limma-trend”

require(limma)
require(edgeR)
dge <- DGEList(counts=GSE198256_count)
design <- model.matrix(~ pDataUSE[,1] )

# Filter
keep <- filterByExpr(dge, design=design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Normalization
dge <- calcNormFactors(dge)

############
# STEP 3.1 LIMMA: TREND

logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

############
# STEP 3.2 LIMMA: VOOM

v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))


###########
#### EXERCISE: COMPARE
###########

# How do we compare limma trend vs limma voom?


# How do we compare DESeq2 vs limma trend?


# How do we compare DESeq2 vs limma voom?



############
# STEP 4: BIOLOGICAL INTERPRETATION
############

# Gene Set Enrichment Analysis.
#    a. ORA.
#    b. GSEA
# How do we bring all the information together at once?


library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library("org.Hs.eg.db")


# 2 APPROACH 1: Gene Enrichment (set treshold)-------------------------------------------------------

res <- results(GSE198256_DESeq2_F, contrast=c('Group','Controls','EarlyRecovery'))##Acute ##LateRecovery ##EarlyRecovery
filtered_res <- res[!is.na(res$log2FoldChange) & 
                      !is.na(res$pvalue) &  !is.na(res$padj)   &                       
                      (abs(res$log2FoldChange) > 1) &                            
                      res$padj < 0.05, ]

filtered_res <- as.data.frame(filtered_res)
signGenes <- rownames(filtered_res)

signGenes_df <- as.data.frame(filtered_res)

ego <- enrichGO(gene = signGenes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                minGSSize = 2,
                maxGSSize = 200,
                readable = TRUE)
dim(ego) ##9
filteredEgo <- simplify(ego) ##simplify so similar bio process eg aging get cluster together
dim(filteredEgo) ##9
dotplot(filteredEgo)


# 3 APPROACH 2: Gene SET Enrichment (no treshold)-------------------------------------------------------
res <- results(GSE198256_DESeq2_F, contrast=c('Group','Controls','EarlyRecovery')) ##EarlyRecovery  ##Acute
               
#filtered_res <- res[!is.na(res$log2FoldChange) & 
#                     !is.na(res$pvalue) &  !is.na(res$padj)   &                       
#                     (abs(res$log2FoldChange) > 1) &                            
#                     res$padj < 0.05, ]
filtered_res <- as.data.frame(res)

geneList <- filtered_res$log2FoldChange ## do a list of all the fold values
geneList
names(geneList) <- as.character(rownames(filtered_res)) ###set name of each item in the list as the entrezid
geneList <- sort(geneList, decreasing = TRUE) ## order so biggest first

keggGO <- gseGO(geneList =geneList,
                OrgDb =  org.Hs.eg.db,
                ont = "ALL",
                nPerm = 1000,
                minGSSize = 20,
                maxGSSize = 200,
                pAdjustMethod = 'fdr',
                pvalueCutoff = 1,
                verbose = FALSE)
dotplot(keggGO)

keggGO

# Review Example:
# https://rpubs.com/jrgonzalezISGlobal/enrichment






