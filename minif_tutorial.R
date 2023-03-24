## https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
## GSE68777


library("minfi")

## Read in the raw .idat data
rgSet <- read.metharray.exp("./example_file/")
rgSet
## Design of 450K probes
manifest <- getManifest(rgSet)
manifest
head(getProbeInfo(manifest))

## preprocess

### Methylation Signal
MSet <- preprocessRaw(rgSet)
MSet
head(getMeth(MSet))
head(getUnmeth(MSet))     
     
### Beta/M value
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet

beta <- getBeta(RSet)


### Map to genome
GRset <- mapToGenome(RSet)
GRset

### Retrieve information
beta <- getBeta(GRset)
M <- getM(GRset)
sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)

### Probe info
annotation <- getAnnotation(GRset)
names(annotation)
write.csv(annotation, '1.csv')


### QC
qc <- getQC(MSet)
plotQC(qc)
densityBeanPlot(MSet)
#### Control probe
controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")


### Sex Prediction
predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
predictedSex

### Normalization
MSet.swan <- preprocessSWAN(rgSet)




 
### DMP WTF is this result???
group <- c(rep('c', 20), rep('t', 20))
beta <- getBeta(MSet.swan)
dmp <- dmpFinder(beta, pheno = group  , type = "categorical")

### Try limma
library(limma)
library(edgeR)

group_info <- factor(group)
design_limma <- model.matrix(~group_info)
limma_object <- DGEList(counts = beta)
limma_object <- calcNormFactors(limma_object)
limma_object <- voom(limma_object,design_limma)
fit <- lmFit(limma_object,design_limma)
fit <- eBayes(fit)
limma_results <- topTable(fit, n = Inf)

### Get info
### Specific bata value
beta[which(rownames(beta)%in%'cg27160931'),]

### Specific probe info 
annotation[which(annotation$Name %in% 'cg27160931'),]



