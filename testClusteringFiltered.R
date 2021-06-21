#This is a test clustering file for data retrieved excluding the '-I' parameter in the fast q dump in retrieveDataSRA.py 
#Filtered Cell Ranger output 
#set up environment and load packages
library(dplyr) 
library(Seurat)
library(patchwork)

#run this line once before running code below 
#reticulate::py_install(packages ='umap-learn')

### SETUP THE SEURAT OBJECT ###

#get full directory to mouse heart GEO data 
current_path<-getwd()

#run in Rstudio (add repository to paths at all locations in code):
Combined_path<-paste(current_path, "/test_combined_cellranger_output/outs/filtered_feature_bc_matrix", sep="")    #path to Combined data 
SAN_path<-paste(current_path, "/test_SAN_cellranger_output/outs/filtered_feature_bc_matrix", sep="")              #path to SAN data 
AVN_path<-paste(current_path, "/test_AVN_cellranger_output/outs/filtered_feature_bc_matrix", sep="")              #path to AVN data 
LPF_path<-paste(current_path, "/test_LPF_cellranger_output/outs/filtered_feature_bc_matrix", sep="")              #path to LPF data 
RPF_path<-paste(current_path, "/test_RPF_cellranger_output/outs/filtered_feature_bc_matrix", sep="")              #path to RPF data 

#load data
Combined.data<-Read10X(Combined_path)
SAN.data<-Read10X(SAN_path)
AVN.data<-Read10X(AVN_path)
LPF.data<-Read10X(LPF_path)
RPF.data<-Read10X(RPF_path)


#initialize the Seurat objects with the raw (non-normalized data)
#combined zones object 
zoneALL<-CreateSeuratObject(counts = Combined.data, project = "zone ALL", min.cells = 3, min.features = 200)

#zone I object 
zoneI<-CreateSeuratObject(counts = SAN.data, project = "zone I", min.cells = 3, min.features = 200)

#zone II object 
zoneII<-CreateSeuratObject(counts = AVN.data, project = "zone II", min.cells = 3, min.features = 200)

#Zone III object 
zoneIIILPF<-CreateSeuratObject(counts = LPF.data, project = "zone IIILPF",min.cells = 3, min.features = 200)
zoneIIIRPF<-CreateSeuratObject(counts = RPF.data, project = "zone IIIRPF", min.cells = 3, min.features = 200)
zoneIII.combined <- merge(zoneIIILPF, y = zoneIIIRPF, add.cell.ids = c("LPF", "RPF"), project = "zone III")




### STANDARD PRE-PROCESSING WORKFLOW ###

#QC and selecting cells for further analysis
zoneALL[["percent.mt"]] <- PercentageFeatureSet(zoneALL, pattern = "^MT-")
zoneI[["percent.mt"]] <- PercentageFeatureSet(zoneI, pattern = "^MT-")
zoneII[["percent.mt"]] <- PercentageFeatureSet(zoneII, pattern = "^MT-")
zoneIII.combined[["percent.mt"]] <- PercentageFeatureSet(zoneIII.combined, pattern = "^MT-")

#visualize QC metrics as a violin plot
VlnPlot(zoneALL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zoneI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zoneII, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zoneIII.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#use featureScatter to visualize feature-feature relationships, also can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(zoneALL, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneALL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(zoneI, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(zoneII, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneII, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(zoneIII.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zoneIII.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

zoneALL <- subset(zoneALL, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
zoneI <- subset(zoneI, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
zoneII <- subset(zoneII, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
zoneIII.combined <- subset(zoneIII.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



### NORMALIZING THE DATA ###

#normalize the data by employing a global-scaling normalization method “LogNormalize”
zoneALL <- NormalizeData(zoneALL, normalization.method = "LogNormalize", scale.factor = 10000)
zoneI <- NormalizeData(zoneI, normalization.method = "LogNormalize", scale.factor = 10000)
zoneII <- NormalizeData(zoneII, normalization.method = "LogNormalize", scale.factor = 10000)
zoneIII.combined <- NormalizeData(zoneIII.combined, normalization.method = "LogNormalize", scale.factor = 10000)




### IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION) ###

#calculate a subset of features that exhibit high cell-to-cell variation in the dataset by directly modeling the mean-variance relationship inherent in single-cell data
zoneALL <- FindVariableFeatures(zoneALL, selection.method = "vst", nfeatures = 2000)
zoneI <- FindVariableFeatures(zoneI, selection.method = "vst", nfeatures = 2000)
zoneII <- FindVariableFeatures(zoneII, selection.method = "vst", nfeatures = 2000)
zoneIII.combined <- FindVariableFeatures(zoneIII.combined, selection.method = "vst", nfeatures = 2000)

#identify the 10 most highly variable genes
top10zoneALL <- head(VariableFeatures(zoneALL), 10)
top10zoneI <- head(VariableFeatures(zoneI), 10)
top10zoneII <- head(VariableFeatures(zoneII), 10)
top10zoneIII <- head(VariableFeatures(zoneIII.combined), 10)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(zoneALL)
plot2 <- LabelPoints(plot = plot1, points = top10zoneALL, repel = FALSE)
plot1 + plot2

plot1 <- VariableFeaturePlot(zoneI)
plot2 <- LabelPoints(plot = plot1, points = top10zoneI, repel = FALSE)
plot1 + plot2

plot1 <- VariableFeaturePlot(zoneII)
plot2 <- LabelPoints(plot = plot1, points = top10zoneII, repel = FALSE)
plot1 + plot2

plot1 <- VariableFeaturePlot(zoneIII.combined)
plot2 <- LabelPoints(plot = plot1, points = top10zoneIII, repel = FALSE)
plot1 + plot2



### SCALING THE DATA ###

#apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
#Shifts the expression of each gene, so that the mean expression across cells is 0, Scales the expression of each gene, so that the variance across cells is 1
all.genes <- rownames(zoneALL)
zoneALL <- ScaleData(zoneALL, features = all.genes)

all.genes <- rownames(zoneI)
zoneI <- ScaleData(zoneI, features = all.genes)

all.genes <- rownames(zoneII)
zoneII <- ScaleData(zoneII, features = all.genes)

all.genes <- rownames(zoneIII.combined)
zoneIII.combined <- ScaleData(zoneIII.combined, features = all.genes)

###  PERFORM LINEAR DIMENSIONAL REDUCTION ###

#perform PCA on the scaled data
zoneALL <- RunPCA(zoneALL, features = VariableFeatures(object = zoneALL))
zoneI <- RunPCA(zoneI, features = VariableFeatures(object = zoneI))
zoneII <- RunPCA(zoneII, features = VariableFeatures(object = zoneII))
zoneIII.combined <- RunPCA(zoneIII.combined, features = VariableFeatures(object = zoneIII.combined))


### DETERMINE THE 'DIMENSIONALITY' OF THE DATASET ###

#NOTE: This process can take a long time for big datasets, comment out for expediency. More
#approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
zoneALL <- JackStraw(zoneALL, num.replicate = 100)
zoneALL <- ScoreJackStraw(zoneALL, dims = 1:15)
JackStrawPlot(zoneALL, dims = 1:15)

zoneI <- JackStraw(zoneI, num.replicate = 100)
zoneI <- ScoreJackStraw(zoneI, dims = 1:15)
JackStrawPlot(zoneI, dims = 1:15)

zoneII <- JackStraw(zoneII, num.replicate = 100)
zoneII <- ScoreJackStraw(zoneII, dims = 1:13)
JackStrawPlot(zoneII, dims = 1:13)

zoneIII.combined <- JackStraw(zoneIII.combined, num.replicate = 100)
zoneIII.combined <- ScoreJackStraw(zoneIII.combined, dims = 1:14)
JackStrawPlot(zoneIII.combined, dims = 1:14)

#heuristic alternative to JackStraw  
ElbowPlot(zoneALL)
ElbowPlot(zoneI)
ElbowPlot(zoneII)
ElbowPlot(zoneIII.combined)



### CLUSTER THE CELLS ###
zoneALL <- FindNeighbors(zoneALL, dims = 1:15)
zoneALL <- FindClusters(zoneALL, resolution = 0.5)
tabALL<-table(Idents(zoneALL))

zoneI <- FindNeighbors(zoneI, dims = 1:15)
zoneI <- FindClusters(zoneI, resolution = 0.5)
tabI<-table(Idents(zoneI))


zoneII <- FindNeighbors(zoneII, dims = 1:13)
zoneII <- FindClusters(zoneII, resolution = 0.5)
tabII<-table(Idents(zoneII))


zoneIII.combined <- FindNeighbors(zoneIII.combined, dims = 1:14)
zoneIII.combined <- FindClusters(zoneIII.combined, resolution = 0.5)
tabIII<-table(Idents(zoneIII.combined))

#look at cluster IDs of the first 5 cells
head(Idents(zoneALL), 5)
head(Idents(zoneI), 5)
head(Idents(zoneII), 5)
head(Idents(zoneIII.combined), 5)


### RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/t-SNE) ###
#Run t-distributed Stochastic Neighbor Embedding
zoneALL_tsne_path<-paste(current_path, "/seurat_output/filtered_test_zoneALL_tsne.jpeg", sep="")    #path to zone ALL tsne
zoneI_tsne_path<-paste(current_path, "/seurat_output/filtered_test_zoneI_tsne.jpeg", sep="")    #path to zone I tsne
zoneII_tsne_path<-paste(current_path, "/seurat_output/filtered_test_zoneII_tsne.jpeg", sep="")    #path to zone II tsne
zoneIII_tsne_path<-paste(current_path, "/seurat_output/filtered_test_zoneIII_tsne.jpeg", sep="")    #path to zone III tsne


zoneALL <- RunTSNE(zoneALL,dims.use = 1:15,reduction.use = "pca")
jpeg(file = zoneALL_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneALL, reduction = "tsne")
dev.off()


zoneI <- RunTSNE(zoneI,dims.use = 1:15,reduction.use = "pca")
jpeg(file = zoneI_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneI, reduction = "tsne")
dev.off()


zoneII <- RunTSNE(zoneII,dims.use = 1:13, reduction.use = "pca")
jpeg(file = zoneII_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneII, reduction = "tsne")
dev.off()

zoneIII.combined <- RunTSNE(zoneIII.combined,dims.use = 1:14, reduction.use = "pca")
jpeg(file = zoneIII_tsne_path)         #save tsne plot as jpeg 
DimPlot(zoneIII.combined, reduction = "tsne")
dev.off()

