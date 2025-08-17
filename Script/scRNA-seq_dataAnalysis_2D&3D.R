#####################################################
##                                                 ##
## SINGLE-CELL SEQUENCING ANALYSIS USING SEURAT V4 ##
##       apply scTransform normalization           ##
##        DEVELOPED BY MARCHA VERHEIJEN            ## 
##                                                 ##
#####################################################

# For function reference list see: https://satijalab.org/seurat/reference/index.html
# For a guide to seurat see: https://satijalab.org/seurat/articles/get_started.html

# normalization: scTransform [replaces NormalizeData(), ScaleData(), and FindVariableFeatures() of the basic workflow].
# Tutorial on integration: youtube https://www.youtube.com/watch?v=HrbeaEJqKcY
# Tutorial on obtaining DEGS: https://www.youtube.com/watch?v=1i6T9hpvwg0

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("dplyr", "patchwork", "glmGamPoi", "tidyr", "metap"))
#library('devtools')
#devtools::install_github("satijalab/seurat", ref = "develop")
#install_github("omicsCore/scTyper")
#BiocManager::install("infercnv")
#install.packages("rjags")
#devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
#remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)

#### Nessesary packages: (some might not be nessesary anymore for this version of the script)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(tidyr)
library(AnnotationHub)
library(ensembldb)
library(multtest)
library(tibble)
library(openxlsx)
library(HGNChelper)
library(tuple)
library(ggplot2)
library(ggridges)
library(SeuratWrappers)

####################################################
#### Settings which need to be adapted by user #####
####################################################

Run_ID="APAP_2D_3D_ControlsMerged"
Additional_info=""
tissue = c("Liver")

BaseDIR="~/Research/Data/Single_Cell/"
dir.create(BaseDIR)   
BaseOutputDIR=paste0(BaseDIR, Run_ID)
dir.create(BaseOutputDIR)   

DataDirStart="~/Research/Data/Single_Cell/" 
DataDirEnd="/outs/raw_feature_bc_matrix/"

##2D & 3D Control samples
Samples=c("1_Ctrl_2D_v2_5k", "4_Ctrl_3D_v2_5k") #control sample gives issues due to different name (no dose)
MergeName<-c("1_Ctrl_2D_ND", "4_Ctrl_3D_ND") #names should have equal amount of "_" to facilitate merge

COMPARE<-"Compound" #set column of Data.all@meta.data that will be used for comparisons
Condition1<-c("2D") #condition of interest
Condition2<-c("3D") #control


# Data filtering thresholds:
nFeature_RNA_Threshold = 500 # tutorial: 500
nCount_RNA_lowerThreshold = 800 # tutorial: 800
	
resolution = 2000 #for clustering 	
pAdjValue = 0.05
L2FC=0

#Celltype annotation database:
scTypeDIR<-paste0(BaseDIR,"/Script/scType/")

#possible db_files to choose:
#db_file<-"ScTypeDB_full.xlsx" #-> Original
#db_file<-"ScTypeDB_full_add_CellMarkerLiver.xlsx" #-> Original + added cell types from CellMarker2.0 (Human, liver, non-cancer, non-embryo)

db_file<-"~/Research/Data/Single_Cell/scType/Annotation_CellTypes_selected_ScTypeDB_CellMarker2Liver.xlsx" #-> subset expected types + MT +APAP

# CellMarker 2.0 DB files are downloaded from: http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html


########################
## Data loading and QC #
########################
OutputDIR=paste0(BaseOutputDIR, "QC/")
dir.create(OutputDIR)
setwd(OutputDIR)
UMAP_DIR=paste0(OutputDIR, "UMAP/")
dir.create(UMAP_DIR)
TSNE_DIR=paste0(OutputDIR, "TSNE/")
dir.create(TSNE_DIR)

nCellsTable<-matrix(NA, nrow=length(Samples), ncol= 2)
colnames(nCellsTable) <- c("BeforeFilter", "AfterBasicFilter")
row.names(nCellsTable)<- MergeName

for (x in 1:length(Samples)) {
        ###############
	## Load data ##
        ###############
		print(paste0("loading data: ", MergeName[x]))
		options(Seurat.object.assay.version="v3")
		count.data<- Read10X(data.dir = paste0(DataDirStart, Samples[x], DataDirEnd))
		seuratObject<-CreateSeuratObject(counts = count.data, project = Run_ID)
		
		nCellsTable[MergeName[x],"BeforeFilter"]<-length(Cells(seuratObject))

		## Adding sample info to metadata
		seuratObject$sample<-MergeName[x]
		seuratObject@meta.data<-separate(seuratObject@meta.data, col="sample", into = c("Sample_ID", "xD", "Compound", "Dose"), sep="_")
		DefaultAssay(seuratObject)<-ASSAY<-"RNA"

		## calculating MT% anf RP%
		seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")
		seuratObject[["percent.rp"]] <- PercentageFeatureSet(seuratObject, pattern = "^RP[SL]")

		assign(paste0("S", x, "_beforeFilter"),seuratObject) # save this object to visualize later

        ####################################
	## Filter data to improve quality ##
        ####################################
		print(paste0("Filter data: ", MergeName[x]))
		Data.filtered<-subset(seuratObject, subset=nCount_RNA >nCount_RNA_lowerThreshold & nFeature_RNA > nFeature_RNA_Threshold)
		
		print(paste0(MergeName[x],": ", length(Cells(Data.filtered)), "/", length(Cells(seuratObject)), " cells left after filtering"))
		seuratObject<-Data.filtered
		assign(paste0("S", x, "_afterBasicFilter"),seuratObject) # save this object to visualize later
		nCellsTable[MergeName[x],"AfterBasicFilter"]<-length(Cells(seuratObject))

	####################################
	## Normalize data for integration ##
        ####################################
		print(paste0("LogNormalize: ", MergeName[x]))

		seuratObject = NormalizeData(seuratObject, normalization.method = 'LogNormalize', scale.factor = 10000)
		assign(paste0("S", x, "_Norm"),seuratObject)

		print(paste0(MergeName[x], " individual processing DONE"))
}

################
## Merge data ##
################
	print(paste0("Merging all samples"))
	ASSAY<-"RNA"
	obj.list<-NULL
	for (x in 1:length(Samples)) {
		obj.list<- append(obj.list, get(paste0("S", x, "_Norm")))
	}

	# Initialize the Seurat object with the normalized data.
	if (length(Samples) == 2) {
		Data.merged = merge(S1_Norm, y = c(S2_Norm),
	 	add.cell.ids = MergeName, project = Run_ID) 
	}
	if (length(Samples) == 3) {
		Data.merged = merge(S1_Norm, y = c(S2_Norm, S3_Norm),
	 	add.cell.ids = MergeName, project = Run_ID) 
	}
		
	# scale and run PCA
	Data.merged<-JoinLayers(Data.merged)
	Data.merged = ScaleData(Data.merged, features = rownames(Data.merged))
	Data.merged = FindVariableFeatures(Data.merged, selection.method = 'vst', nfeatures = 2000)
	Data.merged = RunPCA(Data.merged, features = VariableFeatures(object = Data.merged))
	Data.processed<-Data.merged

	normalized_data <- GetAssayData(Data.processed, assay="RNA", slot="data")
	write.table(normalized_data, "NormData.csv", sep=",")

######################################
# Create plots of previous QC steps ##
######################################
write.table(nCellsTable, paste0(Run_ID, "Ncells_Table.csv"))

	print("Creating plots of previous QC steps")
	stage=c("_beforeFilter", "_afterBasicFilter", "_afterMerge") 
	ASSAY = "RNA"

	for (s in 1:length(stage)) {
		########################
		# Get the data to plot #
		########################
		if (stage[s]!="_afterMerge") {
			for (x in 1:length(Samples)) {
				assign(paste0("S", x), get(paste0("S", x, stage[s])))
			}
			## Merge seurat objects (this is not the same as integration)
			if (length(Samples)==2) {Data.plot<-merge(S1, S2, add.cell.ids=MergeName, project=Run_ID)}
			if (length(Samples)==3) {Data.plot<-merge(S1, c(S2, S3), add.cell.ids=MergeName, project=Run_ID)}
			if (length(Samples)==4) {Data.plot<-merge(S1, c(S2, S3, S4), add.cell.ids=MergeName, project=Run_ID)}
			if (length(Samples)==5) {Data.plot<-merge(S1, c(S2, S3, S4, S5), add.cell.ids=MergeName, project=Run_ID)}
			if (length(Samples)==6) {Data.plot<-merge(S1, c(S2, S3, S4, S5, S6), add.cell.ids=MergeName, project=Run_ID)}
			if (length(Samples)==7) {Data.plot<-merge(S1, c(S2, S3, S4, S5, S6, S7), add.cell.ids=MergeName, project=Run_ID)}
			if (length(Samples)==8) {Data.plot<-merge(S1, c(S2, S3, S4, S5, S6, S7, S8), add.cell.ids=MergeName, project=Run_ID)}
		} 
		if (stage[s]=="_afterMerge") { Data.plot<-Data.merged }



		####################
		## Generate plots ##
		####################
		## Generate plot of nCells per sample
		print(Data.plot@meta.data %>% 
		  	ggplot(aes(x=Dose, fill=Dose)) + 
		  	geom_bar() +
		  	theme_classic() +
		  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
		  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
		  	ggtitle("NCells"))
		savePlot(filename= paste0("QCmetrics", stage[s], "_1_nCellsPerSample.png"),type="png")


		## Visualize QC metrics (before filtering/normalization/integration)
		Data.plot@meta.data$percent.mt[is.nan(Data.plot@meta.data$percent.mt)]<-0
		Data.plot@meta.data$percent.rp[is.nan(Data.plot@meta.data$percent.rp)]<-0
		  print(VlnPlot(Data.plot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), assay=ASSAY, pt.size=0.01, ncol = 4))
			savePlot(filename=paste0("QCmetrics", stage[s],"_2_Vlnplot_AllSamples.png"),type="png")
		  print(VlnPlot(Data.plot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, group.by="Compound", pt.size=0.01, assay=ASSAY))
			savePlot(filename=paste0("QCmetrics", stage[s],"_3_Vlnplot_SplitSamples.png"),type="png")
  		  
		  print(FeatureScatter(Data.plot, feature1 = "nCount_RNA", feature2 = "percent.mt"))
			savePlot(filename=paste0("QCmetrics", stage[s],"_4_FeatureScatter_nCountRNA-percentMt.png"),type="png")
		  print(FeatureScatter(Data.plot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
			savePlot(filename=paste0("QCmetrics", stage[s],"_5_FeatureScatter_nFeatureRNA.png"),type="png")

		#FeatureScatter with MT% colored
		print(Data.plot@meta.data %>% 
		  	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
		  	geom_point() + 
			scale_colour_gradient(low = "gray90", high = "black") +
		  	stat_smooth(method=lm) +
		  	scale_x_log10() + 
		  	scale_y_log10() + 
		  	theme_classic() +
		  	geom_vline(xintercept = nCount_RNA_lowerThreshold) +
		  	geom_hline(yintercept = nFeature_RNA_Threshold) +
		  	facet_wrap(~Dose))
		savePlot(filename=paste0("QCmetrics", stage[s],"_6_FeatureScatter_nCountRNA-nFeatureRNA_colorMT.png"),type="png")

		#FeatureScatter with RP% colored
		print(Data.plot@meta.data %>% 
		  	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.rp)) + 
		  	geom_point() + 
			scale_colour_gradient(low = "gray90", high = "black") +
		  	stat_smooth(method=lm) +
		  	scale_x_log10() + 
		  	scale_y_log10() + 
		  	theme_classic() +
		  	geom_vline(xintercept = nCount_RNA_lowerThreshold) +
		  	geom_hline(yintercept = nFeature_RNA_Threshold) +
		  	facet_wrap(~Dose))
		savePlot(filename=paste0("QCmetrics", stage[s],"_7_FeatureScatter_nCountRNA-nFeatureRNA_colorRP.png"),type="png")
	}

##################################################################
# Identification of highly variable features (feature selection) #
##################################################################

	# Calculate a subset of features that exhibit high cell-to-cell variation in the dataset 
	# (i.e, they are highly expressed in some cells, and lowly expressed in others). 
	# Focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

	DefaultAssay(Data.processed)<-ASSAY<-"RNA"

	#### Identify the 10 most highly variable genes
	top10 <- head(VariableFeatures(Data.processed, assay=ASSAY), 10)
	top25 <- head(VariableFeatures(Data.processed, assay=ASSAY), 25)
	top50 <- head(VariableFeatures(Data.processed, assay=ASSAY), 50)
	write.table(top50, "1_Top50_variableGenes.csv", sep=",")

	#heatmap of top50
	DoHeatmap(Data.processed, top10, group.by="Compound")
		savePlot(filename=paste0("Heatmap_Top10_vaiable features.png"),type="png")
	DoHeatmap(Data.processed, top25, group.by="Compound")
		savePlot(filename=paste0("Heatmap_Top25_vaiable features.png"),type="png")
	DoHeatmap(Data.processed, top50, group.by="Compound")
		savePlot(filename=paste0("Heatmap_Top50_vaiable features.png"),type="png")


	#### Perform linear dimensional reduction for PCA 
	#Data.processed<-RunPCA(Data.processed, features = VariableFeatures(object = Data.processed)) #already done after merge
	VizDimLoadings(Data.processed, dims = 1:2, reduction = "pca")
		savePlot(filename="6a_Loadings.png",type="png")
	DimPlot(Data.processed, reduction = "pca")
		savePlot(filename="6b_PCA.png",type="png")

	#### Determine the 'dimensionality' of the dataset 
	ElbowPlot(Data.processed)
		savePlot(filename="7_ElbowPlot.png",type="png")

###################
# Cell Clustering #
###################

	# first construct a KNN graph based on the euclidean distance in PCA space, 
	# and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).
	# input: selected amount of PCs
	Data.processed <- FindNeighbors(Data.processed, dims = 1:25)

	# modularity optimization techniques: Louvain algorithm (default) or SLM
	# setting the resolution between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
	Data.processed <- FindClusters(Data.processed, resolution = resolution, assay=ASSAY) # determines the amount of clusters to use more precisely then reading the elbowplot
	clusters<-Idents(Data.processed)
	write.table(clusters, "8_CellClusters.csv", sep=",")

	# Non-linear dimensional reduction
	# The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space.

	###############
	# UMAP & TSNE #
	###############
	Data.processed$Dose <- factor(Data.processed$Dose, levels = c("ND", "LD", "HD")) # Set order of the doses

	## UMAP ##
	Data.processed <- RunUMAP(Data.processed, dims = 1:length(levels(clusters)), assay=ASSAY)
		setwd(UMAP_DIR)
		DimPlot(Data.processed, reduction = "umap", group.by= "Compound", cols=c("red", "green"))
			savePlot(filename="9a_UMAP_Clustering_2D3D.png",type="png")

		if(nrow(unique(Data.processed@meta.data["Dose"]))==2) {
		  	print(DimPlot(Data.processed, reduction = "umap", group.by= "Dose", cols=c("grey", "green")))
			savePlot(filename=paste0("9b_UMAP_Clustering_", "Dose", ".png"),type="png")
		}
		if(nrow(unique(Data.processed@meta.data["Dose"]))==3) {
		 	print(DimPlot(Data.processed, reduction = "umap", group.by= "Dose", cols=c("grey", "green", "red")))
			savePlot(filename=paste0("9b_UMAP_Clustering_", "Dose", ".png"),type="png")
		}

		DimPlot(Data.processed, reduction = "umap", group.by= "Sample_ID", cols=c("red", "purple", "blue", "green", "yellow"))
			savePlot(filename="9c_UMAP_Clustering_Sample_ID.png",type="png")
		DimPlot(Data.processed, reduction = "umap")
			savePlot(filename="9d_UMAP_Clustering_Clusters.png",type="png")
		DimPlot(Data.processed, reduction = "umap", label=TRUE)
			savePlot(filename="9e_UMAP_Clustering_Clusters_labels.png",type="png")
		
	## TSNE ##
	Data.processed <- RunTSNE(Data.processed, dims = 1:length(levels(clusters)), assay=ASSAY)
		setwd(TSNE_DIR)
	DimPlot(Data.processed, reduction = "tsne", group.by= "Compound", cols=c("red", "green"))
			savePlot(filename="9a_TSNE_Clustering_Condition.png",type="png")

		if(nrow(unique(Data.processed@meta.data["Dose"]))==2) {
		  	print(DimPlot(Data.processed, reduction = "tsne", group.by= "Dose", cols=c("grey", "green")))
			savePlot(filename=paste0("9b_TSNE_Clustering_Dose.png"),type="png")
		}
		if(nrow(unique(Data.processed@meta.data["Dose"]))==3) {
		 	print(DimPlot(Data.processed, reduction = "tsne", group.by= "Dose", cols=c("grey", "green", "red")))
			savePlot(filename=paste0("9b_TSNE_Clustering_Dose.png"),type="png")
		}

		DimPlot(Data.processed, reduction = "tsne", group.by= "Sample_ID", cols=c("red", "purple", "blue", "green", "yellow"))
			savePlot(filename="9c_TSNE_Clustering_Sample_ID.png",type="png")
		DimPlot(Data.processed, reduction = "tsne")
			savePlot(filename="9d_TSNE_Clustering_Clusters.png",type="png")
		DimPlot(Data.processed, reduction = "tsne", label=TRUE)
			savePlot(filename="9e_TSNE_Clustering_Clusters_labels.png",type="png")

	setwd(OutputDIR)

	#######################################################################
	# Saving Seurat object (Data) for later                               #
	saveRDS(Data.processed, file = paste0(Run_ID, "_processed_Seurat_object_R.rds"))
	####################################################################### 
	# Load the saved Seurat object (Data.processed)
	##Run_ID=Run_ID
	##COMPARE=COMPARE
	##BaseDIR=BaseDIR 
	##OutputDIR=paste0(BaseOutputDIR, QC/")
	##setwd(OutputDIR)
	#Data.processed<-readRDS(paste0(Run_ID, "_processed_Seurat_object_R.rds"))
	## inspect the imported Data.processed
		## str(Data.processed)
		## Data.processed@meta.data
	#######################################################################

#####################################
# Annotate clusters with cell types #
#####################################
	DefaultAssay(Data.processed)<-ASSAY<-"RNA" 

	OutputDIR=paste0(BaseOutputDIR, "Annotation/")
	dir.create(OutputDIR)
	setwd(OutputDIR)
	UMAP_DIR=paste0(OutputDIR, "UMAP/")
	dir.create(UMAP_DIR)
	TSNE_DIR=paste0(OutputDIR, "TSNE/")
	dir.create(TSNE_DIR)


	##########################
	##########################
	## scType for cluster annotation (annotate all samples in 1 go)

	# load gene set preparation function
	source(paste0(scTypeDIR, "/gene_sets_prepare.R"))
	# load cell type annotation function
	source(paste0(scTypeDIR, "/sctype_score_.R"))

	# DB file
	db_ = paste0(scTypeDIR, db_file)
	tissue = tissue

	# prepare gene sets
	gs_list = gene_sets_prepare(db_, tissue)

	# get cell-type by cell matrix
	es.max = sctype_score(scRNAseqData = Data.processed[[ASSAY]]$scale.data, scaled = TRUE, 
		              gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

	# merge by cluster
	cL_results = do.call("rbind", lapply(unique(Data.processed@meta.data$seurat_clusters), function(cl){
	    es.max.cl = sort(rowSums(es.max[ ,rownames(Data.processed@meta.data[Data.processed@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
	    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Data.processed@meta.data$seurat_clusters==cl)), 10)
	}))
	sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

	# set low-confident (low ScType score) clusters to "unknown"
	# lowConfident= sctype_scores$scores < (sctype_scores$ncells/4)
	# HighConfident= sctype_scores$scores < (sctype_scores$ncells/4) --> used for second round annotation with scMAGIC
	sctype_scores$name=sctype_scores$type
	sctype_scores$confidence="Confident"
	sctype_scores$test=sctype_scores$scores/sctype_scores$ncells
	for (r in 1:nrow(sctype_scores)) {
		if (sctype_scores$scores[r] < sctype_scores$ncells[r]/4) { 
			sctype_scores$confidence[r]="lowConfident" 
			sctype_scores$name[r]=sctype_scores$type[r]
	}}
	sctype_scores
	write.table(sctype_scores, paste0("scType_ClusterAnnotationScores.csv"), sep=",")

	# Add CellType name to @meta.data
	Data.processed@meta.data$CellType = ""
	Data.processed@meta.data$RefData = ""
	for(j in unique(sctype_scores$cluster)){
	  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
	  Data.processed@meta.data$CellType[Data.processed@meta.data$seurat_clusters == j] = as.character(cl_type$name[1])
		Data.processed@meta.data$RefData[Data.processed@meta.data$seurat_clusters == j] = as.character(cl_type$confidence[1])
	}
	Data.processed$Cell_Con <-paste0(Data.processed@meta.data[,"CellType"], "_", Data.processed@meta.data[,COMPARE])


##############################
# Plot cell type umap & tsne #
##############################

	setwd(UMAP_DIR)
		DimPlot(Data.processed, reduction = "umap", group.by= "CellType")
			savePlot(filename="11a_UMAP_annotated_cellTypes.png",type="png")
		DimPlot(Data.processed, reduction = "umap", group.by= "CellType", label=TRUE)
			savePlot(filename="11b_UMAP_annotated_cellTypes_labels.png",type="png")
		DimPlot(Data.processed, reduction = "umap", group.by= "Cell_Con")
			savePlot(filename="11c_UMAP_annotated_cellTypes-condition.png",type="png")
		DimPlot(Data.processed, reduction = "umap", group.by= "Cell_Con", label=TRUE)
			savePlot(filename="11d_UMAP_annotated_cellTypes-condition_labels.png",type="png")

	setwd(TSNE_DIR)
		DimPlot(Data.processed, reduction = "tsne", group.by= "CellType")
			savePlot(filename="11a_TSNE_annotated_cellTypes.png",type="png")
		DimPlot(Data.processed, reduction = "tsne", group.by= "CellType", label=TRUE)
			savePlot(filename="11b_TSNE_annotated_cellTypes_labels.png",type="png")
		DimPlot(Data.processed, reduction = "tsne", group.by= "Cell_Con")
			savePlot(filename="11c_TSNE_annotated_cellTypes-condition.png",type="png")
		DimPlot(Data.processed, reduction = "tsne", group.by= "Cell_Con", label=TRUE)
			savePlot(filename="11d_TSNE_annotated_cellTypes-condition_labels.png",type="png")
	setwd(OutputDIR)


Idents(Data.processed) <- Data.processed$CellType

	#######################################################################
	# Saving Seurat object (Data) for later                               #
	saveRDS(Data.processed, file = paste0(Run_ID, "_annotated_Seurat_object_R.rds"))
	####################################################################### 
	# Load the saved Seurat object (Data.processed)
	##Run_ID=Run_ID
	##COMPARE=COMPARE
	##BaseDIR=BaseDIR 
	##OutputDIR=paste0(BaseOutputDIR, "Annotation/")
	##setwd(OutputDIR)
	#Data.processed<-readRDS(paste0(Run_ID, "_annotated_Seurat_object_R.rds"))
	## inspect the imported Data.processed
		## str(Data.processed)
		## Data.processed@meta.data
	#######################################################################






##########################################################
# Get a table with detected genes per Celltype-condition #
##########################################################
OutputDIR=paste0(BaseOutputDIR, "GeneInfo/")
dir.create(OutputDIR)
setwd(OutputDIR)

# add a col to the @meta.data that concatinates the cellType_condition
Idents(Data.processed) <- Data.processed$Cell_Con
Conditions<- unique(Data.processed$Cell_Con)

Genes_Cell_Con<-NULL
for (C in 1:length(Conditions)) {
	print(Conditions[C])
	Select<-NULL
	Select <- subset(Data.processed, subset= Cell_Con == Cell_Con[C])
	SUMexpression<-rowSums(Select@assays$RNA@counts[,,drop=FALSE])
	#ExpressedGenes<-SUMexpression[SUMexpression>0]
	Genes<-SUMexpression

	OUTPUT<-data.frame(matrix(NA, nrow=length(Genes), ncol=10))
	colnames(OUTPUT)<-c("RunID", "Cell_Con", "Celltype", "Dose", "DetectedGene", "AvExpr_posCell", "SDexpr_posCell", "nCells", "nPosCells", "pct.PosCells")
	OUTPUT[,"RunID"]<-rep(Run_ID, length(Genes))
	OUTPUT[,"Cell_Con"]<-rep(Conditions[C], length(Genes))
	OUTPUT[,"Celltype"]<-rep(unlist(strsplit(Conditions[C], "_"))[1], length(Genes))
	OUTPUT[,"Dose"]<-rep(unlist(strsplit(Conditions[C], "_"))[2], length(Genes))

	for (G in 1:length(names(Genes))) {
		GeneInfo<-NULL
		GeneInfo<- as.matrix(Select@assays$RNA@counts[names(Genes)[G],,drop=FALSE])
		OUTPUT[G,"DetectedGene"] <- row.names(GeneInfo)		
		OUTPUT[G,"nCells"] <- ncol(GeneInfo)

		GeneInfo<-GeneInfo[,GeneInfo>0]
			OUTPUT[G,"nPosCells"] <- length(GeneInfo)
			OUTPUT[G,"AvExpr_posCell"] <- mean(GeneInfo)
			OUTPUT[G,"SDexpr_posCell"] <- sd(GeneInfo)
			OUTPUT[G,"pct.PosCells"] <- (OUTPUT[G,"nPosCells"]/OUTPUT[G,"nCells"])*100
	}
	Genes_Cell_Con<-rbind(Genes_Cell_Con, OUTPUT)
}
write.table(Genes_Cell_Con, paste0("14_AllGenes_Cell-Con_", Run_ID, "_fullTable.csv"), sep=",")

GenesDetected_Cell_Con<-Genes_Cell_Con[Genes_Cell_Con$pct.PosCells>10,]
write.table(GenesDetected_Cell_Con, paste0("14_DetectedGenes_Cell-Con_", Run_ID, "_fullTable.csv"), sep=",")





#########################################################################
# Finding differentially expressed features (cluster biomarkers)        #
#########################################################################
#########################################################################
## 2a) Within a single cluster -> what are the DEGs between conditions? #
##             Cluster = CellType per condition                         #
#########################################################################



# Get a list of genes (not DEG)  
testedCelltype<-NULL
CellTypeResponse_All<-NULL

for (COMP in 1:length(Condition1)) {
CellTypeResponse<-data.frame()

	# obtain celltypes 2 test
	C1_types<-Data.processed@meta.data$Cell_Con[grep(paste0("_",Condition1[COMP]), Data.processed@meta.data$Cell_Con)]
	C1<-unique(C1_types[triplicated(C1_types)]) #need at least 3 cells per type
	C1<-sort(gsub(paste0("_",Condition1[COMP]),"", C1))

	C2_types<-Data.processed@meta.data$Cell_Con[grep(paste0("_",Condition2[COMP]), Data.processed@meta.data$Cell_Con)]
	C2<-unique(C2_types[triplicated(C2_types)])
	C2<-sort(gsub(paste0("_",Condition2[COMP]),"", C2))
	merge<-c(C1, C2)
	CellType2test<-merge[duplicated(merge)]


	for (CELLTYPE in CellType2test) {
		testedCelltype<-c(testedCelltype, CELLTYPE)
		genes <-FindMarkers(Data.processed, assay=ASSAY, ident.1=paste0(CELLTYPE, "_", Condition1[COMP]), ident.2=paste0(CELLTYPE, "_", Condition2[COMP]), features=row.names(Data.processed$RNA@data), slot="data", logfc.threshold=0, min.pct=0)
		if (nrow(genes)>0) {
			genes <-cbind(cluster_id = CELLTYPE, Comparison = paste(Condition1[COMP], Condition2[COMP], sep="-"), gene = row.names(genes), genes)
			row.names(genes)<-NULL
		}
		CellTypeResponse<-rbind(CellTypeResponse, genes)
	}
	write.table(CellTypeResponse, paste0("13a_CellTypeResponse_",Condition1[COMP], "-", Condition2[COMP], "_fullTable.csv"), sep=",")
	CellTypeResponse_All<-rbind(CellTypeResponse_All,CellTypeResponse)
}

# Get a list of DEGs  
testedCelltype<-NULL
CellTypeResponse_All<-NULL

for (COMP in 1:length(Condition1)) {
CellTypeResponse<-data.frame()

	# obtain celltypes 2 test
	C1_types<-Data.processed@meta.data$Cell_Con[grep(paste0("_",Condition1[COMP]), Data.processed@meta.data$Cell_Con)]
	C1<-unique(C1_types[triplicated(C1_types)]) #need at least 3 cells per type
	C1<-sort(gsub(paste0("_",Condition1[COMP]),"", C1))

	C2_types<-Data.processed@meta.data$Cell_Con[grep(paste0("_",Condition2[COMP]), Data.processed@meta.data$Cell_Con)]
	C2<-unique(C2_types[triplicated(C2_types)])
	C2<-sort(gsub(paste0("_",Condition2[COMP]),"", C2))
	merge<-c(C1, C2)
	CellType2test<-merge[duplicated(merge)]


	for (CELLTYPE in CellType2test) {
		testedCelltype<-c(testedCelltype, CELLTYPE)
		DEGs <-FindMarkers(Data.processed, assay=ASSAY, ident.1=paste0(CELLTYPE, "_", Condition1[COMP]), ident.2=paste0(CELLTYPE, "_", Condition2[COMP]), slot="data")
		if (nrow(DEGs)>0) {
			DEGs <-cbind(cluster_id = CELLTYPE, Comparison = paste(Condition1[COMP], Condition2[COMP], sep="-"), gene = row.names(DEGs), DEGs)
			row.names(DEGs)<-NULL
		}
		CellTypeResponse<-rbind(CellTypeResponse, DEGs)
	}
	write.table(CellTypeResponse, paste0("13a_CellTypeResponse_",Condition1[COMP], "-", Condition2[COMP], "_fullTable.csv"), sep=",")
	CellTypeResponse_All<-rbind(CellTypeResponse_All,CellTypeResponse)

	DEGs_CellTypeResponse_FDR<-subset(CellTypeResponse,CellTypeResponse$p_val_adj < pAdjValue)
	pct<-cbind(DEGs_CellTypeResponse_FDR$pct.1, DEGs_CellTypeResponse_FDR$pct.2)
	DEGs_CellTypeResponse_FDR<-subset(DEGs_CellTypeResponse_FDR,rowSums(pct) > 0.1)
	write.table(DEGs_CellTypeResponse_FDR, paste0("13b_DEGs_CellTypeResponse_FDR_",Condition1[COMP], "-", Condition2[COMP], ".csv"), sep=",")

	DEGs_CellTypeResponse_FDR<-subset(DEGs_CellTypeResponse_FDR, abs(DEGs_CellTypeResponse_FDR$avg_log2FC)< L2FC)
	write.table(DEGs_CellTypeResponse_FDR, paste0("13b_DEGs_CellTypeResponse_FDR_L2FC_",Condition1[COMP], "-", Condition2[COMP], ".csv"), sep=",")


	print(paste0("Found ", length(unique(DEGs_CellTypeResponse_FDR$gene)), " unique DEGs in the CellTypes when comparing ",Condition1[COMP], "-", Condition2[COMP]))
}


# create overview table of #genes for all celltypes & comparisons
Overview_celltypes<-data.frame(matrix(NA, nrow=length(unique(testedCelltype)), ncol=(length(Condition1)*3)))
row.names(Overview_celltypes)<-unique(testedCelltype)

Cname<-NULL
for (COMP in 1:length(Condition1)) {
	test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
	name<-paste0(test, c("_genes", "_DEGs_FDR", "_DEGs_FDR_L2FC"))
	Cname<-c(Cname, name)
}
colnames(Overview_celltypes)<-Cname

for (COMP in 1:length(Condition1)) {
	test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
	Dose_genes <- CellTypeResponse_All[CellTypeResponse_All$Comparison==test,]

	if (!is.null(Dose_genes)) { 
	for (CELLTYPE in unique(testedCelltype)) {
	C_genes<-Dose_genes[Dose_genes$cluster_id==CELLTYPE,]
	Overview_celltypes[CELLTYPE, paste0(test,"_genes")]<-nrow(C_genes)
	C_genes_p<-subset(C_genes,C_genes$p_val_adj < pAdjValue)
	Overview_celltypes[CELLTYPE, paste0(test,"_DEGs_FDR")]<-nrow(C_genes_p)
	C_genes_pf<-subset(C_genes_p,abs(C_genes_p$avg_log2FC)< L2FC)
	Overview_celltypes[CELLTYPE, paste0(test,"_DEGs_FDR_L2FC")]<-nrow(C_genes_pf)
	}
}
}
write.table(Overview_celltypes, paste0("13c_Overview_celltypes_genes.csv"), sep=",")

# create overview table of DOWNREG #genes for all celltypes & comparisons
Overview_celltypes<-data.frame(matrix(NA, nrow=length(unique(testedCelltype)), ncol=(length(Condition1)*3)))
row.names(Overview_celltypes)<-unique(testedCelltype)

Cname<-NULL
for (COMP in 1:length(Condition1)) {
	test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
	name<-paste0(test, c("_genes", "_DEGs_FDR", "_DEGs_FDR_L2FC"))
	Cname<-c(Cname, name)
}
colnames(Overview_celltypes)<-Cname

for (COMP in 1:length(Condition1)) {
	test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
	Dose_genes <- CellTypeResponse_All[CellTypeResponse_All$Comparison==test,]
	Dose_genes <- Dose_genes[Dose_genes$avg_log2FC<0,]

	if (!is.null(Dose_genes)) { 
	for (CELLTYPE in unique(testedCelltype)) {
	C_genes<-Dose_genes[Dose_genes$cluster_id==CELLTYPE,]
	Overview_celltypes[CELLTYPE, paste0(test,"_genes")]<-nrow(C_genes)
	C_genes_p<-subset(C_genes,C_genes$p_val_adj < pAdjValue)
	Overview_celltypes[CELLTYPE, paste0(test,"_DEGs_FDR")]<-nrow(C_genes_p)
	C_genes_pf<-subset(C_genes_p,abs(C_genes_p$avg_log2FC)< L2FC)
	Overview_celltypes[CELLTYPE, paste0(test,"_DEGs_FDR_L2FC")]<-nrow(C_genes_pf)
	}
}
}
write.table(Overview_celltypes, paste0("13d_Overview_celltypes_DOWNREG_genes.csv"), sep=",")





#########################################################################
## 2b) Within a single cluster -> what are the DEGs between conditions? #
##             Cluster = seurat_cluster number per condition            #
#########################################################################


# add a col to the @meta.data that concatinates the cellType_condition
Data.processed$Clust_Con <-paste0(Data.processed@meta.data[,"seurat_clusters"], "_", Data.processed@meta.data[,COMPARE])
Idents(Data.processed) <- Data.processed$Clust_Con

##################################################
# Get a list of genes (not DEG) for each cluster #
##################################################
	testedCluster<-NULL
	ClusterResponse_All<-NULL

	for (COMP in 1:length(Condition1)) {
		ClusterResponse<-data.frame()

		# obtain clusters 2 test
	C1_clusters<-Data.processed@meta.data$Clust_Con[grep(paste0("_",Condition1[COMP]), Data.processed@meta.data$Clust_Con)]
		C1<-unique(C1_clusters[triplicated(C1_clusters)]) #need at least 3 cells per cluster
		C1<-sort(as.integer(gsub(paste0("_",Condition1[COMP]),"", C1)))

		C2_clusters<-Data.processed@meta.data$Clust_Con[grep(paste0("_",Condition2[COMP]), Data.processed@meta.data$Clust_Con)]
		C2<-unique(C2_clusters[triplicated(C2_clusters)])
		C2<-sort(as.integer(gsub(paste0("_",Condition2[COMP]),"", C2)))
		merge<-c(C1, C2)
		cluster2test<-merge[duplicated(merge)]
		for (CLUSTER in cluster2test) {
			testedCluster<-c(testedCluster, CLUSTER)
			genes <-FindMarkers(Data.processed, assay=ASSAY, ident.1=paste0(CLUSTER, "_", Condition1[COMP]), ident.2=paste0(CLUSTER, "_", Condition2[COMP]), features=row.names(Data.processed$RNA@data), slot="data", logfc.threshold=0, min.pct=0)
			if (nrow(genes)>0) {
			genes <-cbind(cluster_id = CLUSTER, Comparison = paste(Condition1[COMP], Condition2[COMP], sep="-"), gene = row.names(genes), genes)
			row.names(genes)<-NULL
			}
			ClusterResponse<-rbind(ClusterResponse, genes)
		}
		write.table(ClusterResponse, paste0("12a_ClusterResponse_",Condition1[COMP], "-", Condition2[COMP], "_fullTable.csv"), sep=",")
		ClusterResponse_All<-rbind(ClusterResponse_All,ClusterResponse)
	}

##########################################
## Get a list of "DEGs" for each cluster #
##########################################
	testedCluster<-NULL
	ClusterResponse_All<-NULL
	for (COMP in 1:length(Condition1)) {
		ClusterResponse<-data.frame()

		# obtain clusters 2 test
	C1_clusters<-Data.processed@meta.data$Clust_Con[grep(paste0("_",Condition1[COMP]), Data.processed@meta.data$Clust_Con)]
		C1<-unique(C1_clusters[triplicated(C1_clusters)]) #need at least 3 cells per cluster
		C1<-sort(as.integer(gsub(paste0("_",Condition1[COMP]),"", C1)))

		C2_clusters<-Data.processed@meta.data$Clust_Con[grep(paste0("_",Condition2[COMP]), Data.processed@meta.data$Clust_Con)]
		C2<-unique(C2_clusters[triplicated(C2_clusters)])
		C2<-sort(as.integer(gsub(paste0("_",Condition2[COMP]),"", C2)))
		merge<-c(C1, C2)
		cluster2test<-merge[duplicated(merge)]
		
		for (CLUSTER in cluster2test) {
			testedCluster<-c(testedCluster, CLUSTER)
			DEGs <-FindMarkers(Data.processed, assay=ASSAY, ident.1=paste0(CLUSTER, "_", Condition1[COMP]), ident.2=paste0(CLUSTER, "_", Condition2[COMP]), slot="data")
			if (nrow(DEGs)>0) {
			DEGs <-cbind(cluster_id = CLUSTER, Comparison = paste(Condition1[COMP], Condition2[COMP], sep="-"), gene = row.names(DEGs), DEGs)
			row.names(DEGs)<-NULL
			}
			ClusterResponse<-rbind(ClusterResponse, DEGs)
		}
		ClusterResponse_All<-rbind(ClusterResponse_All,ClusterResponse)

		DEGs_ClusterResponse_FDR<-subset(ClusterResponse,ClusterResponse$p_val_adj < pAdjValue)
		pct<-cbind(DEGs_ClusterResponse_FDR$pct.1, DEGs_ClusterResponse_FDR$pct.2)
		DEGs_ClusterResponse_FDR<-subset(DEGs_ClusterResponse_FDR,rowSums(pct) > 0.1)
		write.table(DEGs_ClusterResponse_FDR, paste0("12b_DEGs_ClusterResponse_FDR_",Condition1[COMP], "-", Condition2[COMP], ".csv"), sep=",")

		DEGs_ClusterResponse_FDR_L2FC<-subset(DEGs_ClusterResponse_FDR, abs(DEGs_ClusterResponse_FDR$avg_log2FC)< L2FC)
		write.table(DEGs_ClusterResponse_FDR, paste0("12b_DEGs_ClusterResponse_FDR_L2FC_",Condition1[COMP], "-", Condition2[COMP], ".csv"), sep=",")

		print(paste0("Found ", length(unique(DEGs_ClusterResponse_FDR$gene)), " unique DEGs in the Clusters when comparing ",Condition1[COMP], "-", Condition2[COMP]))
	}

	# create overview table of #genes for all clusters & comparisons
	Overview_clusters<-data.frame(matrix(NA, nrow=length(unique(testedCluster)), ncol=(length(Condition1)*3)))
	testedClusterName<-paste0("C", testedCluster)
	row.names(Overview_clusters)<-unique(testedClusterName)

	Cname<-NULL
	for (COMP in 1:length(Condition1)) {
		test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
		name<-paste0(test, c("_genes", "_DEGs_FDR", "_DEGs_FDR_L2FC"))
		Cname<-c(Cname, name)
	}
	colnames(Overview_clusters)<-Cname

	for (COMP in 1:length(Condition1)) {
		test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
		Dose_genes <- ClusterResponse_All[ClusterResponse_All$Comparison==test,]

		if (!is.null(Dose_genes)) { 
		for (CLUSTER in unique(testedCluster)) {
		C_genes<-Dose_genes[Dose_genes$cluster_id==CLUSTER,]
		Overview_clusters[paste0("C",CLUSTER), paste0(test,"_genes")]<-nrow(C_genes)
		C_genes_p<-subset(C_genes,C_genes$p_val_adj < pAdjValue)
		Overview_clusters[paste0("C",CLUSTER), paste0(test,"_DEGs_FDR")]<-nrow(C_genes_p)
		C_genes_pf<-subset(C_genes_p,abs(C_genes_p$avg_log2FC)< L2FC)
		Overview_clusters[paste0("C",CLUSTER), paste0(test,"_DEGs_FDR_L2FC")]<-nrow(C_genes_pf)
		}
	}
	}
	write.table(Overview_clusters, paste0("12c_Overview_clusters_genes.csv"), sep=",")

	# create overview table of DOWNREG #genes for all clusters & comparisons
	Overview_clusters<-data.frame(matrix(NA, nrow=length(unique(testedCluster)), ncol=(length(Condition1)*3)))
	row.names(Overview_clusters)<-unique(testedClusterName)

	Cname<-NULL
	for (COMP in 1:length(Condition1)) {
		test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
		name<-paste0(test, c("_genes", "_DEGs_FDR", "_DEGs_FDR_L2FC"))
		Cname<-c(Cname, name)
	}
	colnames(Overview_clusters)<-Cname

	for (COMP in 1:length(Condition1)) {
		test<- paste(Condition1[COMP], Condition2[COMP], sep="-")
		Dose_genes <- ClusterResponse_All[ClusterResponse_All$Comparison==test,]
		Dose_genes <- Dose_genes[Dose_genes$avg_log2FC<0,]
		Dose_genes <- Dose_genes[Dose_genes$Comparison==test,]

		if (!is.null(Dose_genes)) { 
		for (CLUSTER in unique(testedCluster)) {
		C_genes<-Dose_genes[Dose_genes$cluster_id==CLUSTER,]
		Overview_clusters[paste0("C",CLUSTER), paste0(test,"_genes")]<-nrow(C_genes)
		C_genes_p<-subset(C_genes,C_genes$p_val_adj < pAdjValue)
		Overview_clusters[paste0("C",CLUSTER), paste0(test,"_DEGs_FDR")]<-nrow(C_genes_p)
		C_genes_pf<-subset(C_genes_p,abs(C_genes_p$avg_log2FC)< L2FC)
		Overview_clusters[paste0("C",CLUSTER), paste0(test,"_DEGs_FDR_L2FC")]<-nrow(C_genes_pf)
		}
	}
	}
	write.table(Overview_clusters, paste0("12d_Overview_clusters_DOWNREG_genes.csv"), sep=",")





#########

#######################################################################
# Saving Seurat object (Data) for later                               #
saveRDS(Data.processed, file = paste0(Run_ID, "_processed_Seurat_object_annotated_R.rds"))
####################################################################### 



print("End of Script -> storing SessionInfo")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

