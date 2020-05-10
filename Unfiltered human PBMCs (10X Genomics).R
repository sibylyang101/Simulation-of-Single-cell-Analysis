### Orchestrating Single Cell Analysis
### Chapter 25 Unfiltered human PBMCs (10X Genomics)
### peripheral blood mononuclear cell (PBMC) dataset from 10X Genomics website

## data loading
library(BiocFileCache) # manage file across sessions
bfc <- BiocFileCache("raw_data", ask = FALSE) #represents the location of files stored on disk
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE) #Creates a SingleCellExperiment

library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol) #e.g combine gene symbol with Ensembl

library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID") #gets the mapped ids (column) for a set of keys
View(location)


##quality control
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc)) #drop empty droplet
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

unfiltered <- sce.pbmc

#Compute per-cell quality control metrics 
#remove cells with large mitochondrial proportions ( = cell damage)
#risk of removing cell types with low RNA content
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT"))) 
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher") # one-tail
sce.pbmc <- sce.pbmc[,!high.mito]

summary(high.mito)

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- high.mito

gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="subsets_Mito_percent",
              colour_by="discard") + ggtitle("Mito percent"),
  ncol=2
) #Set up a gtable layout to place multiple grobs on a page

#Proportion of mitochondrial reads in each cell of the PBMC dataset compared to its total count. 
plotColData(unfiltered, x="sum", y="subsets_Mito_percent",
            colour_by="discard") + scale_x_log10()


#normalization
library(scran)
set.seed(1000)
#cluster using either log-expression values or ranks
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters) #Normalization by deconvolution of size factor
sce.pbmc <- logNormCounts(sce.pbmc)

summary(sizeFactors(sce.pbmc)) #sets the size factors for all cells

par("mar")
par(mar=c(1,1,1,1))
plot(librarySizeFactors(sce.pbmc), sizeFactors(sce.pbmc), pch=16,
     xlab="Library size factors", ylab="Deconvolution factors", log="xy")


## variance modeling
set.seed(1001)
#Model the per-gene variance with Poisson noise
#decomposing it into technical and biological components
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

plot(dec.pbmc$mean, dec.pbmc$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")
curfit <- metadata(dec.pbmc)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)

# Dimensionality reduction
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

ncol(reducedDim(sce.pbmc, "PCA"))

## clustering
##？？？？ Here's error
g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
library(tables)
colLabels(sce.pbmc) <- factor(clust)

table(colLabels(sce.pbmc))

plotTSNE(sce.pbmc, colour_by="label")

## interpretation
markers <- findMarkers(sce.pbmc, pval.type="some", direction="up")
marker.set <- markers[["7"]]
as.data.frame(marker.set[1:30,1:3])
plotExpression(sce.pbmc, features=c("CD14", "CD68",
                                    "MNDA", "FCGR3A"), x="label", colour_by="label")