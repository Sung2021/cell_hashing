### This is a chunk of codes for 
### Cell_hashing analysis
#####################################
setwd('~/Desktop/HMH/')
library(dplyr)
library(ggplot2)
library(reshape)
library(Seurat)

## 1. read tag data 
## data from CITE-seq-Count
hto <- Read10X('rds/2022.cell_hashing/umi_count/', gene.column = 1)
hto %>% dim()

## 2. read transcriptome data
dge <- read.table('rds/2022.cell_hashing/dge/star_gene_exon_tagged.dge.txt.gz',
                  header = T, row.names = 1, check.names = FALSE)
dge %>% dim()

## 3. Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(colnames(hto), colnames(dge))
length(cells.use) 

## 4. Create Seurat Object
obj.srt <- CreateSeuratObject(dge[,cells.use], min.features = 100)
## the number of tag: 2
obj.srt[['HTO']] <- CreateAssayObject(counts = hto[1:2,colnames(obj.srt)])
obj.srt <- HTODemux(obj.srt, assay = 'HTO')
table(obj.srt$HTO_classification.global)

## 5. selecting singlets only
# singlet <- rownames(obj.srt@meta.data[obj.srt@meta.data$HTO_classification.global == "Singlet",])
# singlet <- rownames(obj.srt@meta.data[obj.srt@meta.data$HTO_classification.global %in% c( "Singlet",
                                                                                          "Negative"),])
obj.srt.filtered <- subset(obj.srt, cells= singlet)
dim(obj.srt.filtered) 
## save obj at this point
saveRDS(obj.srt.filtered, 'rds/2022.cell_hashing/cell_hashing.72hr.singlets.negative_included.rds')

## sample prepare for QC
obj.srt <- obj.srt.filtered
hist(obj.srt@meta.data$nCount_RNA, breaks = 100)
summary(obj.srt@meta.data$nCount_RNA)
boxplot(obj.srt@meta.data$nCount_RNA)

## 5.1. QC step
# Create the knee plot 
# Convert to dgCMatrix, which is a compressed, sparse matrix format
mat <- as((obj.srt@assays$RNA@counts), "dgCMatrix")
mat %>% str()
mat <- obj.srt@assays$RNA@counts %>% as.matrix() %>% t()

tot_counts <- rowSums(mat)  ### UMI per cell
Rank <- rank(-tot_counts)
df <- cbind(tot_counts, Rank)
df <- data.frame(df)
df %>% dim()
tag <- substr(obj.srt@meta.data$HTO_classification,1,2)
df$tag <- tag

## knee plot
ggplot(df, aes(tot_counts,Rank, color=tag)) + geom_point() +
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count")

ggplot(df, aes(tot_counts,Rank, color=tag)) + geom_point() +
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count") +facet_grid(.~tag)

ggplot(df, aes(tot_counts,Rank)) + geom_point() +
  scale_x_log10() +  annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count") + 
  geom_vline(xintercept = c(300,400,1000), color='red')


#################################################################
## 6. Analysis
obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst',nfeatures = 3000) 
all.genes <- rownames(obj.srt)
obj.srt <- ScaleData(obj.srt, features = all.genes)
obj.srt <- NormalizeData(obj.srt)
obj.srt <- RunPCA(obj.srt, features = VariableFeatures(object = obj.srt), npcs = 10)
DimPlot(obj.srt, reduction = "pca")

obj.srt <- FindNeighbors(obj.srt, dims = 1:10)
obj.srt <- FindClusters(obj.srt, resolution = c(0.1,0.2,0.3,0.4,0.5,1))
obj.srt <- RunTSNE(obj.srt, dims = 1:10)
obj.srt <- RunUMAP(obj.srt, dims = 1:10)
