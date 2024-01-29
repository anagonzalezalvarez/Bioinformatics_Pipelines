

library(NOISeq)
library(tidyverse)
library(biomaRt)
library(DESeq2)
library(GEOquery)

setwd("~/Desktop/PhD_2nd/BESE398_Pipelines/Bioinformatics_Pipelines/W1/2_NewDataset")
### Load Data ----------------------------------------------------------------------------------

## Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249911

GSE132408_counts <- read.csv("data/raw/GSE132408_TET2_RawCounts.csv")
ncol(GSE132408_counts)
colnames(GSE132408_counts)
nrow(GSE132408_counts)


### Metadata ----------------------------------------------------------------------------------

# load metadata table from GEO
gds <- getGEO("GSE132408")
## https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
## Retrieve data
Meta_GSE132408 <- pData(gds$GSE132408_series_matrix.txt.gz@phenoData)
colnames(Meta_GSE132408)
## Filter for only the columns i am interested in
Meta_GSE132408 <- Meta_GSE132408[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1")]
## arranging data in description in the right way
Meta_GSE132408$description_rev <- stri_reverse(Meta_GSE132408$description)
Meta_GSE132408$description_rev_split <- sapply(strsplit(Meta_GSE132408$description_rev, "\\."), function(x) ifelse(length(x) > 2, x[3], ""))
Meta_GSE132408$description_rev_split <- sapply(strsplit(Meta_GSE132408$description_rev, "\\."), function(x) ifelse(length(x) > 2, paste(x[-(1:2)], collapse = "."), ""))
Meta_GSE132408$mydescription <- stri_reverse(Meta_GSE132408$description_rev_split)


### NOISESeq ----------------------------------------------------------------------------------
# We need length, gc, biotype, chromosome

##### BIOMART

### Set X as rownames
rownames(GSE132408_counts) <- GSE132408_counts[,1] 
GSE132408_counts$X <- NULL ## remove x col (we already have it as rownames)

### Get rownames which contain entrez id as csv
entrez_id <- rownames(GSE132408_counts)

# Save to  file
write.table(entrez_id, file="data/biomart/entrez_id_GSE132408.csv",col.names = FALSE,row.names = FALSE,quote=F)

#### Go to biomart and upload the file to retrive length, gc, biotype, chromosome

annotgene <- read.table("data/biomart/mart_export_GSE132408.txt",  sep = "\t", header = TRUE)
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
colnames(summary) <- c("filter","n_rows", "n_unique", "n_in_GSE132408_counts")

summary[1,1] <- "original"
summary[1,2] <- nrow(annotgene)
summary[1,3] <- length(unique(annotgene$Entrezgene))
summary[1,4] <- sum(rownames(GSE132408_counts) %in% annotgene$Entrezgene) # How many genes do I get annotated?


# Filter the information
annotgene <- annotgene[annotgene$Chromosome %in% c(as.character(1:22) ,"X","Y"),]

summary[2,1] <- "chromosome filter"
summary[2,2] <- nrow(annotgene)
summary[2,3] <- length(unique(annotgene$Entrezgene))
summary[2,4] <- sum(rownames(GSE132408_counts) %in% annotgene$Entrezgene) # How many genes do I get annotated?


## Multiple entries for one entrezid -> need to keep only 1 to be able to set it as rownames
annotgene_filt <- annotgene[!duplicated(annotgene$Entrezgene),] ### remove duplicates

summary[3,1] <- "duplicate removal"
summary[3,2] <- nrow(annotgene_filt)
summary[3,3] <- length(unique(annotgene_filt$Entrezgene))
summary[3,4] <- sum(rownames(GSE132408_counts) %in% annotgene_filt$Entrezgene) # How many genes do I get annotated?


## Remove nas
annotgene_filt <- na.omit(annotgene_filt) ##27,234 remove na's in annotgene_filt$Entrezgene

summary[4,1] <- "nas removal"
summary[4,2] <- nrow(annotgene_filt)
summary[4,3] <- length(unique(annotgene_filt$Entrezgene))
summary[4,4] <- sum(rownames(GSE132408_counts) %in% annotgene_filt$Entrezgene) # How many genes do I get annotated?


## Match annotation annotgene_filt$Entrezgene with counts GSE132408_counts
rownames(annotgene_filt) <- as.character(annotgene_filt$Entrezgene) ### set annotgene_filt$Entrezgene as rowname 

##  Filter GSE132408_counts for only the ones present in annotgene_filt (the ones i retrieved the biomart info)
GSE132408_counts_filt <- GSE132408_counts[rownames(GSE132408_counts) %in% rownames(annotgene_filt),]
nrow(GSE132408_counts_filt)

## Also save the excluded ones :)
GSE132408_counts_exc <-GSE132408_counts[!(rownames(GSE132408_counts) %in% rownames(annotgene_filt)),]
## reorder annotgene_filt to the order in GSE132408_counts_filt
annotgene_ord <- annotgene_filt[rownames(GSE132408_counts_filt ),]
## now both annotgene_ord and GSE132408_counts_filt contain same entrez id in the same order

### NOISEQ WORKFLOW ---------------------------------------------------------------------------

### Get factors order of the condition col
Factors_GSE132408 <- data.frame(Meta_GSE132408[,c("mydescription")])

#Factors_GSE132408 <- data.frame(Meta_GSE132408 [ colnames(GSE132408_counts_filt),c("description")])
colnames(Factors_GSE132408)[1]<- "Group" ### RENAME 

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
data_NOISEQ <- readData(data = GSE132408_counts_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE132408)

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

save(data_NOISEQ,GSE132408_counts_filt,annotgene_ord,file="GSE132408_step1.Rda")

### NORMALIZATION  --------------------------------------------------------------

myRPKM = rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA = uqua(assayData(data_NOISEQ)$exprs, long = lengthuse, lc = 0.5, k = 0)
myTMM = tmm(assayData(data_NOISEQ)$exprs, long = 1000, lc = 0)

### DESEQ --------------------------------------------------------------
## Need:
## 1. Counts: GSE132408_counts_filt
## 2. Sample info: data_NOISEQ@phenoData@data

# 1 prepare data ----------------------------------------------------------------------------
## Make sure row names in samples match column names in counts
### all columns in counts in varaibles rows?
all(colnames(GSE132408_counts_filt) %in% rownames(data_NOISEQ@phenoData@data)) ##T
### same order?
all(colnames(GSE132408_counts_filt) == rownames(data_NOISEQ@phenoData@data)) ##T

# Rename my conditions just so its easier to understand
pDataUSE <- data_NOISEQ@phenoData@data
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="CL"] <- "CTRL"
pDataUSE[pDataUSE=="CL.IFNG"] <- "CTRL.IFNγ"
pDataUSE[pDataUSE=="KO"] <- "TET2.KO"
pDataUSE[pDataUSE=="KO.IFNG"] <- "TET2.KO.IFNγ"


pDataUSE[,1] <- as.factor(pDataUSE[,1]) ## as factors


# 2 Create DESeqDataSet object --------------------------------------------------------------------------
GSE132408_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE132408_counts_filt,
                                           colData = pDataUSE,
                                           design = ~ Group) ## in colData what is the name of the col that contains condition


### In case we do not want to use intercept in the design use: design = ~ -1 + Group

# 3 QC and Filtering  ---------------------------------------------------------------------------
# Keeping only rows with at least 10 reads
smallestGroupSize <- 6
keep <- rowSums(counts(GSE132408_DESeq2) >= 10) >= smallestGroupSize
GSE132408_DESeq2_F <- GSE132408_DESeq2[keep,]
nrow(GSE132408_DESeq2) ##18062
nrow(GSE132408_DESeq2_F) ###11192

# 4 Factors  ---------------------------------------------------------------------------
GSE132408_DESeq2_F$Group ## Levels: CTRL CTRL + IFNγ TET2 KO TET2 KO + IFNγ
GSE132408_DESeq2_F$Group <-  relevel(GSE132408_DESeq2_F$Group, ref = 'CTRL') ## set CTRL as ref

# 5 Run DESeq ---------------------------------------------------------------------------
GSE132408_DESeq2_F<- DESeq(GSE132408_DESeq2_F)
resultsNames(GSE132408_DESeq2_F) ###"Group_CTRL.IFNγ_vs_CTRL"    "Group_TET2.KO_vs_CTRL"      "Group_TET2.KO.IFNγ_vs_CTRL"


# 6 Explore Results ALL ---------------------------------------------------------------------------
summary(results(GSE132408_DESeq2_F, contrast =c('Group','CTRL.IFNγ','CTRL')))
#LFC > 0 (up)       : 1442, 10%
#LFC < 0 (down)     : 1093, 7.8%


# 6 FILTER alpha = 0.05 ---------------------------------------------------------------------------

res_1 <- results(GSE132408_DESeq2_F, contrast =c('Group','CTRL.IFNγ','CTRL'))
res_1_df <- as.data.frame(res_1)
res_1_df_filt <- res_1_df[!is.na(res_1_df$log2FoldChange) & !is.na(res_1_df$pvalue) &
                            (abs(res_1_df$log2FoldChange) > 2) &
                            res_1_df$padj < 0.05, ]
summary(res_1_df_filt)
nrow(res_1_df_filt) ###342
# VOLCANO PLOT ---------------------------------------------------------------------------

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
res <- results(GSE132408_DESeq2_F, contrast =c('Group','CTRL','CTRL.IFNγ'))

topT <- as.data.frame(res)


#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)


# 5 reRun DESeq ---------------------------------------------------------------------------


GSE132408_DESeq2_F$Group <-  relevel(GSE132408_DESeq2_F$Group, ref = 'TET2.KO') ## set TET2.KO as ref

GSE132408_DESeq2_F<- DESeq(GSE132408_DESeq2_F)
resultsNames(GSE132408_DESeq2_F) ###"Group_CTRL.IFNγ_vs_CTRL"    "Group_TET2.KO_vs_CTRL"      "Group_TET2.KO.IFNγ_vs_CTRL"

## "CTRL""CTRL.IFNγ""TET2.KO""TET2.KO.IFNγ"


summary(results(GSE132408_DESeq2_F, contrast =c('Group','TET2.KO.IFNγ','TET2.KO')))
res_2 <- results(GSE132408_DESeq2_F, contrast =c('Group','TET2.KO.IFNγ','TET2.KO'))
res_2_df <- as.data.frame(res_4)
res_2_df_filt <- res_2_df[!is.na(res_2_df$log2FoldChange) & !is.na(res_2_df$pvalue) &
                            (abs(res_2_df$log2FoldChange) > 2) &
                            res_2_df$padj < 0.05, ]
nrow(res_2_df_filt)
###12



###











### 6.1 Contrasts ---------------------------------------------------------------------------
resultsNames(GSE132408_DESeq2_F) ##[1] "Intercept"                    "Group_CTRL...IFNγ_vs_CTRL"    "Group_TET2.KO_vs_CTRL"       [4] "Group_TET2.KO...IFNγ_vs_CTRL"
plotMA(GSE132408_res, ylim=c(-2,2))


### 7. Shrinkage ---------------------------------------------------------------------------
res_lfcShrink <- lfcShrink(GSE132408_DESeq2_F,coef=c("Group_CTRL...IFNγ_vs_CTRL"))
res_lfcShrink <- lfcShrink(GSE132408_DESeq2_F,coef=c("Group_Covid196Mo_vs_Covid193Mo"))
plotMA(res_lfcShrink, ylim=c(-2,2))


#GSE132408_DESeq2_F <- DESeq(GSE132408_DESeq2_F, test="LRT", reduced=~1)
#GSE132408_DESeq2_res_LRT <- results(GSE132408_DESeq2_F)
#GSE132408_DESeq2_res_LRT
#res <- results(GSE132408_DESeq2_res_LRT)

# Technical replicates?

# How to interpret the results?

plotCounts(GSE132408_DESeq2_F, gene="100287102", intgroup="Group")

## STEP 3.1.4: QC??

# How do visualize?
vsd <- vst(GSE132408_DESeq2_F, blind=FALSE)
rld <- rlog(GSE132408_DESeq2_F, blind=FALSE)
head(assay(vsd), 3)

# heatmap
library("pheatmap")
select <- order(rowMeans(counts(GSE132408_DESeq2_F,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(GSE132408_DESeq2_F)[,c("Group")])
colnames(df) <- "Group"

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)

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
dge <- DGEList(counts=GSE132408_counts)
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

BiocManager::install("topGO")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Review Example:
# https://rpubs.com/jrgonzalezISGlobal/enrichment






