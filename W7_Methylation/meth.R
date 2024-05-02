
## GEO ID: GSE188573
## Met data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188573

## Tutorial: https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

setwd('/Users/gonzalac/Desktop/PhD_2nd/BESE398_Pipelines/Bioinformatics_Pipelines/W7_Methylation/results/imag')
### DMA WORKFLOW -------------------------
## 1. Obtaining the data


# read in the sample sheet for the experiment
dataDirectory <- "data/GSE188573_RAW"

targets <- read.metharray.sheet(dataDirectory, pattern="Data_sheet.csv")
targets$Basename <- unique(sub('data/GSE188573_RAW 2/', '', targets$Basename))
head(targets)

# Read the data
rgSet <- read.metharray.exp(base = dataDirectory, targets = targets)

## Change the sample names to be more informative
targets$Sample_Group <- gsub("Healthy Donor", "HD", targets$Sample_Group)
targets$ID <- paste0(targets$Sample_Group,".", targets$Sample_Name)
sampleNames(rgSet) <- targets$ID

## 3. Quality Control

### 3.1 pvals 
detP <- detectionP(rgSet, type = "m+u")
head(detP)

pdf("1_detection_pvals.pdf")
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")
dev.off()

### 3.2 report
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")


### 3.3 Filtering
# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

## 4. Normalization
mSetSq <- preprocessQuantile(rgSet) 

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation
pdf("2_normalization.pdf")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()


## 5. Data Exploration
pdf("3_data_exploration.pdf")
# MDS plots to look at largest sources of variation
plotMDS(getM(mSetSq), top=1000, gene.selection="common", labels = NULL,
        col=pal[factor(targets$Sample_Group)], cex = 0.5)
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

dev.off()

## 6. Filtering

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt


# Remove probes on the sex chromosomes
# creating the annotation
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
table(keep)
#keep
#FALSE   TRUE 
#18683 811890 
mSetSqFlt <- mSetSqFlt[keep,]

# plot MDS after filtering
pdf("4_afterfilt_MDS.pdf")
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", labels = NULL,
        col=pal[factor(targets$Sample_Group)], cex = 0.5)
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

## 7. Probe-wise differential methylation analysis
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

pdf("5_m-b-values.pdf")
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()


cellType <- factor(targets$Sample_Group, levels = c('HD', 'COVID19'))
design <- model.matrix(~ 0 + cellType, data = targets)
colnames(design) <- gsub("cellType", "", colnames(design))


# linear model
fit <- lmFit(mVals, design)

contMatrix <- makeContrasts(COVID19-HD,
                            levels=design)
contMatrix       

#          Contrasts
# Levels    COVID19 - HD
#   HD                -1
#   COVID19            1

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))
#        COVID19 - HD
# Down          14334
# NotSig       776658
# Up            20898


# get the table of results for the first contrast (COVID19 - HD)
annSub <- annEPIC[match(rownames(mVals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]
# DMP results
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annSub)
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# plot the top 4 most significantly differentially methylated CpGs 
pdf("6_top4_cpgs.pdf")
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})
dev.off()


## 8. Region-wise differential methylation analysis
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "COVID19 - HD", arraytype = "EPIC")
str(myAnnotation)

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
png("7_dvp_top_dmr.png", width = 200, height = 400, units = "mm", res = 600)
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")
dev.off()

#### 8. Gene Ontology Testing


sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
sigCpGs[1:10]# First 10 significant CpGs

length(sigCpGs) ##35232

all <- DMPs$Name
length(all)

pdf("7_DMP.pdf")
par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)
dev.off()

topGSA(gst, number=10)


## Differential Variability ------------------------
detP <- detectionP(rgSet, type = "m+u")
mSetSq <- preprocessQuantile(rgSet)


targets$Sex <- getSex(mSetSq)$predictedSex
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

keep <- rowSums(detP < 0.01) == ncol(detP) 
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))


pdf('8_diff_variability.pdf')
pal <- brewer.pal(8,"Set1")
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], labels=targets$Sex, 
        main="With Sex CHR Probes")
legend("topleft", legend=levels(factor(targets$Sample_Group)), 
       text.col=pal)

plotMDS(getM(mSetSqFlt[keep,]), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], labels=targets$Sex, 
        main="Without Sex CHR Probes")
legend("top", legend=levels(factor(targets$Sample_Group)),
       text.col=pal)
dev.off()



keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]



# get M-values for analysis
mVals <- getM(mSetSqFlt)

design <- model.matrix(~factor(targets$Sample_Group)) 
fitvar <- varFit(mVals, design = design, coef = c(1,2))

# Summary of differential variability
summary(decideTests(fitvar))

# > summary(decideTests(fitvar))
#        (Intercept) factor(targets$Sample_Group)HD
# Down             0                              4
# NotSig          19                         811571
# Up          811871                            315


topDV <- topVar(fitvar, coef=2)
topDV


bVals <- getBeta(mSetSqFlt)


pdf("9_top4_dv_cpgs.pdf")
par(mfrow=c(2,2))
sapply(rownames(topDV)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, 
          ylab = "Beta values")
})
dev.off()








