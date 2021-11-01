### This is a chunk of codes for 
### 2021 
### Cell_hashing analysis
#####################################
setwd('~/Desktop/HMH/')
library(dplyr)
library(ggplot2)
library(reshape)
library(Seurat)

tmp <- Read10X('~/Desktop/Sung_work/fastq/cell_hashing/19092-44/02/umi_count/', gene.column = 1)
tmp[1:3,]
dge <- read.table('~/Desktop/Sung_work/fastq/cell_hashing/19092-44/temp/star_gene_exon_tagged.dge.txt.gz',
                  header = T, row.names = 1, check.names = FALSE)
dge[1:3,1:3] # 18174 10000

# Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(colnames(tmp), colnames(dge))
length(cells.use) # 8209

obj.srt <- CreateSeuratObject(dge[,cells.use], min.features = 100)
obj.srt[['HTO']] <- CreateAssayObject(counts = tmp[1:3,colnames(obj.srt)])
obj.srt <- HTODemux(obj.srt, assay = 'HTO')
table(obj.srt$HTO_classification.global)

obj.srt@meta.data[1:3,]
###########################################################################
### Filtering 
singlet <- rownames(obj.srt@meta.data[obj.srt@meta.data$HTO_classification.global == "Singlet",])
obj.srt.filtered <- subset(obj.srt, cells= singlet)
dim(obj.srt.filtered)
obj.srt.backup <- obj.srt
saveRDS(obj.srt.backup, 'rds/2021.10.14.cell_hashing.all.before.filter.rds')
obj.srt <- obj.srt.filtered
hist(obj.srt@meta.data$nCount_RNA, breaks = 100)
summary(obj.srt@meta.data$nCount_RNA)
boxplot(obj.srt@meta.data$nCount_RNA)



# Create the knee plot 
# Convert to dgCMatrix, which is a compressed, sparse matrix format
mat<- as(t(obj.srt@assays$RNA@counts), "dgCMatrix")
tot_counts <- rowSums(mat)  ### UMI per cell
Rank <- rank(-tot_counts)
df <- cbind(tot_counts, Rank)
df <- data.frame(df)
df[1:3,]
ggplot(df, aes(tot_counts,Rank)) + geom_point() +
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count")

ggplot(df, aes(tot_counts,Rank)) + geom_point() +
  scale_x_log10() +  annotation_logticks() +
  labs(y = "Barcode rank", x = "Total UMI count") + 
  geom_vline(xintercept = c(300,400,1000), color='red')


obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst',nfeatures = 3000) 
all.genes <- rownames(obj.srt)
obj.srt <- ScaleData(obj.srt, features = all.genes)
obj.srt <- NormalizeData(obj.srt)
obj.srt <- RunPCA(obj.srt, features = VariableFeatures(object = obj.srt), npcs = 10)
DimPlot(obj.srt, reduction = "pca", group.by = 'tag')

obj.srt <- FindNeighbors(obj.srt, dims = 1:10)
obj.srt <- FindClusters(obj.srt, resolution = c(0.1,0.2,0.3,0.4,0.5,1))
obj.srt <- RunTSNE(obj.srt, dims = 1:10)
obj.srt <- RunUMAP(obj.srt, dims = 1:10)

group <- 'RNA_snn_res.0.1'
group <- 'tag'
group <- 'Phase'
reduction = 'tsne'
reduction ='umap'
DimPlot(obj.srt, reduction = reduction, group.by = group)

FeaturePlot(obj.srt, features = c('Tcf7','Ctla4','Id3','Il7r','Bcl2','Nr4a1','Il17a','Il21','Il23a'), reduction = reduction)
DimPlot(obj.srt, group.by = "RNA_snn_res.0.1")
colnames(obj.srt@meta.data)
DimPlot(obj.srt, group.by  = "HTO_classification", 
        split.by = "HTO_classification", ncol = 3)
genes <- c('Bcl6', 'Cxcr5')
group <- "HTO_classification"
FeaturePlot(obj.srt, features = c('Bcl6', 'Cxcr5'))
FeaturePlot(obj.srt, features = genes, split.by = "HTO_classification", ncol = 3, pt.size = 0.3)
VlnPlot(obj.srt, features = genes, group.by =group, ncol = 1)


saveRDS(obj.srt, 'rds/2021.10.13.cell_hashing.rds')




obj.srt.filtered <- subset(obj.srt, nCount_RNA >= 1000)
dim(obj.srt.filtered)
obj.srt <- obj.srt.filtered


obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst') # in the paper
all.genes <- rownames(obj.srt)
obj.srt <- ScaleData(obj.srt, features = all.genes)
obj.srt <- NormalizeData(obj.srt)
obj.srt <- RunPCA(obj.srt, features = VariableFeatures(object = obj.srt))
ElbowPlot(obj.srt)
DimPlot(obj.srt, reduction = "pca")

obj.srt <- FindNeighbors(obj.srt, dims = 1:10)
obj.srt <- FindClusters(obj.srt, resolution = c(0.1,0.2,0.3,0.4,0.5,1,1.5))
obj.srt <- RunTSNE(obj.srt, dims = 1:10)
obj.srt <- RunUMAP(obj.srt, dims = 1:10)

DimPlot(obj.srt, reduction = c('pca','tsne','umap')[3], group.by = 'HTO_maxID')
DimPlot(obj.srt, reduction = 'tsne', group.by = 'RNA_snn_res.0.3')

genes <- c('Tcf7', 'Ctla4','Bcl6','Icos','Pdcd1','Cxcr5','Slamf1')
FeaturePlot(obj.srt, features = genes, reduction = 'tsne', ncol = 4, pt.size = 0.1)
VlnPlot(obj.srt, features = genes, group.by = 'tag',ncol = 4, pt.size = 0.1)


obj.srt[['tag']] <- 'NA'
obj.srt@meta.data[grep('mt',obj.srt@meta.data$HTO_maxID),'tag'] <- 'mt'
obj.srt@meta.data[grep('wt',obj.srt@meta.data$HTO_maxID),'tag'] <- 'wt'
obj.srt@meta.data[grep('naive',obj.srt@meta.data$HTO_maxID),'tag'] <- 'naive'
obj.srt@meta.data$tag <- factor(obj.srt@meta.data$tag, levels = c('naive','wt','mt'))


saveRDS(obj.srt, 'rds/2021.10.14.cell_hashing.singlet_only.umi_filtered.rds')
obj.srt <- readRDS('rds/2021.10.14.cell_hashing.singlet_only.umi_filtered.rds')
Idents(obj.srt) <- 'tag'
all.markers <- FindAllMarkers(obj.srt, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
all.markers[1:3,]
top_genes <- all.markers %>% group_by(cluster) %>% top_n(20, wt=avg_log2FC)

DoHeatmap(obj.srt, features = top_genes$gene) + NoLegend() + 
  theme(text = element_text(size=6)) 

write.csv(top_genes, 'rds/2021.cell_hashing/2021.10.21.cell_hashing.res.0.5.top.gene.20.csv')

Idents(obj.srt) <- 'RNA_snn_res.0.1'
all.markers <- FindAllMarkers(obj.srt, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
all.markers[1:3,]
top_genes <- all.markers %>% group_by(cluster) %>% top_n(30, wt=avg_log2FC)

DoHeatmap(obj.srt, features = top_genes$gene) + NoLegend() + 
  theme(text = element_text(size=8)) 
group='tag'
group='RNA_snn_res.0.5'
DoHeatmap(obj.srt, features = top_genes$gene, group.by = group) + NoLegend() + 
  theme(text = element_text(size=8)) 


genes <- c('Nr4a1','Nr4a2','Nr4a3','Nr3c1','Tsc22d3','Fkbp5')
FeaturePlot(obj.srt, features = genes, reduction = 'tsne', ncol = 4, pt.size = 0.1)
VlnPlot(obj.srt, features = genes, group.by = 'RNA_snn_res.0.1',ncol = 4, pt.size = 0.1)
VlnPlot(obj.srt, features = genes, group.by = 'tag',ncol = 4, pt.size = 0.1)



############ 2021.10.20 #########
obj.srt <- readRDS('rds/2021.cell_hashing/2021.10.14.cell_hashing.singlet_only.umi_filtered.rds')
wt.srt <- readRDS('rds/2021.cell_hashing/2021.10.17.cell_hashing.wt.srt.rds')

reduction='tsne'
DimPlot(wt.srt, reduction = reduction)
group = 'RNA_snn_res.0.1'
DimPlot(wt.srt, reduction = reduction,group.by = group)
gene <- c('Cd4','Cd8a','Icos','Sell','Cd44')
FeaturePlot(wt.srt, features = gene, reduction = 'tsne',pt.size = 0.1)



reduction='tsne'
obj <- obj.srt
group <- 'tag'
DimPlot(obj, reduction = reduction, group= group)

meta <- obj.srt@meta.data
meta[1:3,]
ggplot(meta, aes(tag)) + geom_bar()
group='RNA_snn_res.0.5'
ggplot(meta, aes(tag)) + geom_bar(aes_string(fill=group)) + theme_classic()
ggplot(meta, aes(tag)) + geom_bar(aes_string(fill=group),position = 'fill') + theme_classic()

group <- 'RNA_snn_res.0.5'
group <- 'tag'
DimPlot(obj.srt, reduction = 'tsne', group.by = group) + theme_classic()



############ 2021.10.21 #########
obj.srt <- readRDS('rds/2021.cell_hashing/2021.10.14.cell_hashing.singlet_only.umi_filtered.rds')

obj <- obj.srt
reduction='tsne'
group = 'RNA_snn_res.0.1'

DimPlot(obj, reduction = reduction,group.by = group)
gene <- c('Cd4','Cd8a','Icos','Sell','Cd44')
genes <- read.csv('~/Desktop/HMH/info/gene_info/gene_info.csv', header = T)
gene <- genes$chemotaxis
gene <- genes[,2]
gene <- genes[,3]
FeaturePlot(obj, features = gene, reduction = 'tsne',pt.size = 0.1)



library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)
useMart()
getLDS
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  return(humanx)
}

g2m <- convertHumanGeneList(cc.genes$g2m.genes)
s.genes <- convertHumanGeneList(cc.genes$s.genes)
list(g2m,s.genes)
saveRDS(list(g2m,s.genes), 'rds/2021.cell_hashing/cell_cycle.genes.rds')

cellcycle.genes <- readRDS('rds/2021.cell_hashing/cell_cycle.genes.rds')
obj.srt <- CellCycleScoring(obj.srt, g2m.features = cellcycle.genes[[1]],
                        s.features = cellcycle.genes[[2]])
obj.srt@meta.data %>% head()

group='Phase'
group='RNA_snn_res.0.5'
DimPlot(obj.srt, reduction = reduction,group.by = group, label = T)
gene <- c('Cdkn1a','Ccne2','Atf7')
FeaturePlot(obj, features = gene, reduction = 'tsne',pt.size = 0.1)



group <- 'tag'
Idents(obj.srt) <- group
all.markers <- FindAllMarkers(obj.srt, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
saveRDS(all.markers, 'rds/2021.cell_hashing/2021.10.21.cell_hashing_res.0.5.all.markers.rds')
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='naive') %>% dim()
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='wt') %>% dim()
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='mt') %>% dim()


all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='naive') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='wt') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='mt') %>% dim()

gene1 <- all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='naive') %>% rownames()
gene2 <- all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='wt') %>% rownames()
gene3 <- all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='mt') %>% rownames()


top_genes <- all.markers %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC)

DoHeatmap(obj.srt, features = c(gene1,gene2,gene3)) + NoLegend() + 
  theme(text = element_text(size=8))
DoHeatmap(obj.srt, features = gene2) + NoLegend() + 
  theme(text = element_text(size=8))
DoHeatmap(obj.srt, features = gene3) + NoLegend() + 
  theme(text = element_text(size=8))

top_genes %>% filter(cluster=='6')
top_genes %>% filter(cluster=='10')



library(scater)
library(SingleCellExperiment)
library(monocle)


obj.srt <- obj.srt
data <- as(as.matrix(obj.srt@assays$RNA@data), 'sparseMatrix')

## preparing pd,fData,fd from seurat object
pd <- new('AnnotatedDataFrame', data = obj.srt@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

## monocle set
m.data <- newCellDataSet(data, phenoData = pd, featureData = fd, 
                         lowerDetectionLimit = 0.5, 
                         expressionFamily = VGAM::negbinomial.size())

## size factor
## distribution factor
m.data <- estimateSizeFactors(m.data)
m.data <- estimateDispersions(m.data)

m.data <- detectGenes(m.data, min_expr = 0.1) 


#Run ordering algorithm


var_genes <- obj.srt@assays$RNA@var.features
ordering_genes <- var_genes
HSMM <- setOrderingFilter(m.data, ordering_genes) ## from here, HSMM

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM <- reduceDimension(HSMM,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=4,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

# First decide what you want to color your cells by
print(head(pData(HSMM)))

HSMM <- orderCells(HSMM)
group='tag'
group='Phase'
group='RNA_snn_res.0.5'
group='RNA_snn_res.0.1'
plot_cell_trajectory(HSMM, 
                     color_by = group,
                     theta = -15,
                     show_branch_points = F,
                     show_tree = T,
                     cell_size = 1) +
  theme(legend.position = "right") + facet_wrap(.~tag)
saveRDS(HSMM, 'rds/2021.cell_hashing/2021.10.22.cell_hashing.all.monocle.rds')
temp.s <- data.frame(HSMM@reducedDimS)
temp.s <- t(temp.s)
colnames(temp.s) <- c("C1", "C2", "C3", "C4")
temp.s <- data.frame(temp.s)
gene.exp <- data.frame(obj.srt@assays$RNA@scale.data)
gene <- 'Nr4a1'
ggplot(temp.s, aes(C1,C2, color=t(gene.exp[gene,]))) + geom_point()
obj.monocle <- readRDS('rds/2021.cell_hashing/')

genes
gene <- c('Cxcr5','Il2ra','Prdm1','Id2')
gene <-rownames(obj.srt)[grep('^Dep', rownames(obj.srt))]
gene <- 'Tsc2'
FeaturePlot(obj.srt, features = gene, reduction = reduction)
rownames(obj.srt)
rownames(obj.srt)[grep('^Dep', rownames(obj.srt))]


selected.markers <- FindMarkers(obj.srt, ident.1 = 'wt', ident.2 = 'mt',
                                only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
selected.markers %>% dim()
selected.markers <- FindMarkers(obj.srt, ident.1 = 'mt', ident.2 = 'wt',
                                only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25)
gene <- selected.markers %>% rownames()

DoHeatmap(obj.srt, features = gene) + NoLegend() 

table(obj.srt@meta.data$tag, obj.srt@meta.data$RNA_snn_res.0.5)
df <- table(obj.srt@meta.data$tag, obj.srt@meta.data$RNA_snn_res.0.5) %>% data.frame()
ggplot(df, aes(Var2,Freq, color=Var1)) + geom_col()


############ 2021.10.22 #########
obj.srt <- readRDS('rds/2021.cell_hashing/2021.10.14.cell_hashing.singlet_only.umi_filtered.rds')
obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst',nfeatures = 3000) 
all.genes <- rownames(obj.srt)
obj.srt <- ScaleData(obj.srt, features = all.genes)
obj.srt <- NormalizeData(obj.srt)
obj.srt <- RunPCA(obj.srt, features = VariableFeatures(object = obj.srt), npcs = 10)
DimPlot(obj.srt, reduction = "pca", group.by = 'tag')

obj.srt <- FindNeighbors(obj.srt, dims = 1:10)
obj.srt <- FindClusters(obj.srt, resolution = c(0.1,0.2,0.3,0.4,0.5,1))
obj.srt <- RunTSNE(obj.srt, dims = 1:10)
obj.srt <- RunUMAP(obj.srt, dims = 1:10)
cellcycle.genes <- readRDS('rds/2021.cell_hashing/cell_cycle.genes.rds')
obj.srt <- CellCycleScoring(obj.srt, g2m.features = cellcycle.genes[[1]],
                            s.features = cellcycle.genes[[2]])

saveRDS(obj.srt, 'rds/2021.10.21.cell_hashing.vst3000.rds')

obj.srt <- readRDS('rds/2021.10.21.cell_hashing.vst3000.rds')
group <- 'RNA_snn_res.0.1'
group <- 'tag'
group <- 'Phase'
reduction = 'tsne'
reduction ='umap'
reduction ='pca'
DimPlot(obj.srt, reduction = reduction, group.by = group)
FeaturePlot(obj.srt, features = c('Tcf7','Ctla4','Id3','Il7r','Bcl2',
                                  'Nr4a1','Nr4a2','Nr4a3','Il17a','Il21','Il23a','Nr3c1',
                                  'Ccne2','Ifng','Cxcr5','Ccr5'), 
            reduction = reduction, pt.size = 0.1)

genes <- read.csv('info/gene_info/gene_info.csv', header = T)
gene <- c(genes$Tfh, 'Cxcr5')
gene <- c(genes$Th17, 'Cxcr5')
gene <- c(genes$Th1, 'Cxcr5')
gene <- c(genes$Th2, 'Cxcr5')
gene <- c(genes$iTreg, 'Cxcr5')
gene <- c(genes$Tr, 'Cxcr5')
FeaturePlot(obj.srt, features = gene, 
            reduction = reduction, pt.size = 0.2)

FeaturePlot(obj.srt, features = c('Tcf7','Ctla4','Id3','Il7r','Bcl2','Nr4a1','Il17a','Il21','Il23a'), reduction = reduction)
DimPlot(obj.srt, group.by = "RNA_snn_res.0.1")
colnames(obj.srt@meta.data)
DimPlot(obj.srt, group.by  = "HTO_classification", 
        split.by = "HTO_classification", ncol = 3)
genes <- c('Bcl6', 'Cxcr5')
group <- "HTO_classification"
FeaturePlot(obj.srt, features = c('Bcl6', 'Cxcr5'))
FeaturePlot(obj.srt, features = genes, split.by = "HTO_classification", ncol = 3, pt.size = 0.3)
VlnPlot(obj.srt, features = genes, group.by =group, ncol = 1)


saveRDS(obj.srt, 'rds/2021.10.13.cell_hashing.rds')

VariableFeatures(obj.srt)
VariableFeaturePlot(obj.srt)


genes <- c('Nr4a1','Nr4a2','Nr4a3','Nr3c1','Tsc22d3','Fkbp5')
FeaturePlot(obj.srt, features = genes, reduction = 'tsne', ncol = 4, pt.size = 0.1)
VlnPlot(obj.srt, features = genes, group.by = 'RNA_snn_res.0.1',ncol = 4, pt.size = 0.1)
VlnPlot(obj.srt, features = genes, group.by = 'tag',ncol = 4, pt.size = 0.1)

var.features <- VariableFeatures(obj.srt)
data.data <- obj.srt@assays$RNA@data[var.features,] %>% data.frame()
data.data <- t(data.data) %>% data.frame()
data.data %>% dim() 
data.data %>% head()
input.data <- data.data
## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}

zscore_result <- zscore(input.data  = input.data)

ggplot(data.data, aes(Cdkn1a, Ccne2)) + geom_point()

pheatmap::pheatmap(zscore_result, cluster_rows = T, cluster_cols = T, show_rownames = F, show_colnames = F)

cor.data <- cor(data.data[1:100,1:100], method = 'pearson')
heatmap(cor.data)



tmp <- obj.srt@assays$RNA@data

var.features <- VariableFeatures(obj.srt)
library(Rmagic)
library(ggplot2)

MAGIC_data <- magic(tmp, knn = 5, k = 5)
magic.data <- t(MAGIC_data$result)

obj.reductions <- obj.srt@reductions$tsne@cell.embeddings
library(viridis)
gene <-'Cxcr5'
gene <- 'Il6'
gene <- 'Il21'
gene <- 'Bcl6'
gene <- 'Ccr7'
gene <- 'Il12a'
gene <- 'Tbx21'

ggplot(data.frame(obj.reductions), aes(tSNE_1, tSNE_2, color=data.frame(magic.data)[,gene])) + geom_point(size=0.5) + 
  scale_color_viridis(option="B") + theme_classic() +ggtitle(gene)
ggplot(data.frame(magic.data), aes(Tcf7, Lef1)) + geom_point()
ggplot(data.frame(magic.data), aes(Cxcr5, Ccr5,
                                   color=obj.srt@meta.data$tag)) + 
  geom_point(size=0.5) + theme_classic() + facet_grid(.~obj.srt@meta.data$tag)

saveRDS(magic.data, 'rds/2021.cell_hashing/2021.10.22.naive.wt.mt.all.magic.rds')



Idents(obj.srt) <- 'RNA_snn_res.0.1'
all.markers <- FindAllMarkers(obj.srt, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
all.markers[1:3,]
all.markers %>% dim()
saveRDS(all.markers, 'rds/2021.cell_hashing/2021.10.22.all.markers.rds')

all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='0') %>% dim()
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='1') %>% dim()
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='2') %>% dim()
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='3') %>% dim()


all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='0') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='1') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='2') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='3') %>% dim()

top_genes <- all.markers %>% group_by(cluster) %>% top_n(15, wt=avg_log2FC)

DoHeatmap(obj.srt, features = rownames(obj.srt)[grep('Il',rownames(obj.srt))]) + NoLegend() + 
  theme(text = element_text(size=8)) 

obj.srt@assays$RNA@data
temp <- obj.srt@assays$RNA@data[rownames(obj.srt)[grep('Il',rownames(obj.srt))],]
pheatmap::pheatmap(temp)



rownames(obj.srt)[grep('Il',rownames(obj.srt))]

rownames(obj.srt)[grep('Stat',rownames(obj.srt))]

reduction='tsne'
gene <- rownames(obj.srt)[grep('Il6',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Stat',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Il7',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Gata',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Cxcr',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Ccr',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Cx',rownames(obj.srt))]
gene <- rownames(obj.srt)[grep('Bcl',rownames(obj.srt))]

# <- cellcycle.genes[[2]]
gene <-rownames(obj.srt)[grep('Cd40',rownames(obj.srt))]
gene <- c('Il2ra','Stat5a','Bcl2l1')
gene <-rownames(obj.srt)[c(grep('Il12',rownames(obj.srt)),
                           grep('Stat4',rownames(obj.srt)),
                           grep('Tbx21',rownames(obj.srt)))]
gene <-rownames(obj.srt)[c(grep('Il4',rownames(obj.srt)),
                           grep('Stat6',rownames(obj.srt)),
                           grep('Gata3',rownames(obj.srt)))]
gene <-rownames(obj.srt)[c(grep('Il6',rownames(obj.srt)),
                           grep('Stat3',rownames(obj.srt)),
                           grep('Ror',rownames(obj.srt)))]
FeaturePlot(obj.srt, features = gene, reduction = reduction, pt.size = 0.2, ncol = 3)



gene <-rownames(obj.srt)[c(grep('Inf',rownames(obj.srt)),
                           grep('Gzmb',rownames(obj.srt)),
                           grep('Il12r',rownames(obj.srt)))]
gene <-rownames(obj.srt)[c(grep('Il4',rownames(obj.srt)),
                           grep('Il13',rownames(obj.srt)),
                           grep('Il5',rownames(obj.srt)),
                           grep('Ccl1',rownames(obj.srt)))]
gene <-rownames(obj.srt)[c(grep('Il17r',rownames(obj.srt)),
                           grep('Il23',rownames(obj.srt)),
                           grep('Ifn',rownames(obj.srt)),
                           grep('Il12r',rownames(obj.srt)))]
FeaturePlot(obj.srt, features = gene, reduction = reduction, pt.size = 0.1, ncol = 5)

reduction ='tsne'
group='RNA_snn_res.0.5'
DimPlot(obj.srt, reduction = reduction, group.by = group, label = T)


Idents(obj.srt) <- group
all.markers <- FindAllMarkers(obj.srt, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
top_genes <- all.markers %>% group_by(cluster) %>% top_n(20, wt=avg_log2FC)

DoHeatmap(obj.srt, features = top_genes$gene) + NoLegend() + 
  theme(text = element_text(size=6)) 


all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='0') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='5') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='1') %>% dim()
all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='3') %>% dim()



all.markers %>% filter(avg_log2FC >= 0.5849625 | avg_log2FC <= -0.5849625) %>%  filter(cluster=='5') %>% rownames()
all.markers %>% filter(avg_log2FC >= 1 | avg_log2FC <= -1) %>%  filter(cluster=='5') %>% rownames()
cluster <- '6'
all.markers %>% filter(cluster == cluster) %>% filter(avg_log2FC >=1) %>% hist()
all.markers %>% filter(cluster == cluster) %>% filter(avg_log2FC >=1) %>% rownames()
gene <- all.markers %>% filter(cluster == cluster) %>% filter(avg_log2FC >= 1) %>% rownames()
DoHeatmap(obj.srt, features = gene) + NoLegend() + 
  theme(text = element_text(size=6)) 
gene <- c('Klrg1','Itgae','Id2','Fos')
FeaturePlot(obj.srt, features = gene, reduction = reduction, ncol = 5)

saveRDS(obj.srt, 'rds/2021.cell_hashing/2021.10.22.all.rds')




############ 2021.10.25 #########
obj.srt <- readRDS('rds/2021.cell_hashing/2021.10.22.all.rds')
reduction='tsne'
group = 'RNA_snn_res.0.5'
group='Phase'

DimPlot(obj.srt, reduction = reduction,group.by = group, label = T)
gene <- c('Cd4','Cd8a','Icos','Sell','Cd44','Cd69','Stat3','Lamtor1','Rheb',rownames(obj.srt)[grep('mtor', rownames(obj.srt))])
genes <- read.csv('~/Desktop/HMH/info/gene_info/gene_info.csv', header = T)
gene <- genes$chemotaxis
gene <- genes[,2]
gene <- genes[,3]
gene <- c('Cd4','Il21','Icos','Sell','Cd44','Cd69','Stat3','Cxcr5','Ctla4')

FeaturePlot(obj.srt, features = gene, reduction = 'tsne',pt.size = 0.1)


rownames(obj.srt)[grep('mtor', rownames(obj.srt))]
rownames(obj.srt)[grep('Rheb', rownames(obj.srt))]

magic.data <- readRDS('rds/2021.cell_hashing/2021.10.22.naive.wt.mt.all.magic.rds')
HSMM <-  readRDS('rds/2021.cell_hashing/2021.10.22.cell_hashing.all.monocle.rds')
temp.s <- data.frame(HSMM@reducedDimS)
temp.s <- t(temp.s)
colnames(temp.s) <- c("C1", "C2", "C3", "C4")
temp.s <- data.frame(temp.s)
gene.exp <- data.frame(obj.srt@assays$RNA@scale.data)
meta <- data.frame(obj.srt@meta.data)
gene <- 'Il21'
gene <- 'Tnfrsf4'
gene <- 'Fos'
ggplot(temp.s, aes(C1,C2, color=t(gene.exp[gene,]))) + 
  geom_point(size=0.5, alpha=0.5) + theme_classic() + 
  scale_color_gradient(low = 'grey', high = 'dark blue') + ggtitle(gene)

gene <- 'Cxcr5'
ggplot(temp.s, aes(C1,C2, color=magic.data[,gene])) + 
  geom_point(size=0.5, alpha=0.5) + theme_classic() + 
  scale_color_gradient(low = 'grey', high = 'dark blue') + ggtitle(gene)

group <- 'Phase'
ggplot(temp.s, aes(C1,C2, color=meta[,group])) + 
  geom_point(size=0.2) + theme_classic()



ggplot(data.frame(magic.data), aes(Nr4a1)) + geom_density(aes(color=meta$tag))



meta <- data.frame(meta)

p1<- ggplot(meta, aes(tag, fill=RNA_snn_res.0.1)) + geom_bar(position = 'fill')
p2 <-ggplot(meta, aes(tag, fill=RNA_snn_res.0.2)) + geom_bar(position = 'fill')
p3 <-ggplot(meta, aes(tag, fill=RNA_snn_res.0.3)) + geom_bar(position = 'fill')
p4 <-ggplot(meta, aes(tag, fill=RNA_snn_res.0.4)) + geom_bar(position = 'fill')
p5 <-ggplot(meta, aes(tag, fill=RNA_snn_res.0.5)) + geom_bar(position = 'fill')
p6 <-ggplot(meta, aes(tag, fill=RNA_snn_res.1)) + geom_bar(position = 'fill')
cowplot::plot_grid(p1,p2,p3,p4,p5,p6, ncol = 3)
ggplot(meta, aes(RNA_snn_res.0.4,fill=tag)) + geom_bar(position = 'fill')



gene <- c('Cd4','Il21','Icos','Sell','Cd44','Cd69','Stat3','Cxcr5','Ctla4')
rownames(obj.srt)[(grep('Cdk', rownames(obj.srt)))]
rownames(obj.srt)[(grep('Pik', rownames(obj.srt)))]
rownames(obj.srt)[(grep('Nf', rownames(obj.srt)))]
rownames(obj.srt)[(grep('Il2', rownames(obj.srt)))]
rownames(obj.srt)[(grep('Il6', rownames(obj.srt)))]
gene <- c('Cd28', "Pik3ca","Nfkbia"  , "Nfkbib" ,  "Nfkbid"  , "Nfkbie" )

gene <- c('Cd28', "Pik3ca",'Akt1','Cdk1','Cdk2','Cdk4','Cdk9', 'Cdkn1b', 'Ccne2')
gene <- c(rownames(obj.srt)[(grep('Il2', rownames(obj.srt)))], rownames(obj.srt)[(grep('Il6', rownames(obj.srt)))])
FeaturePlot(obj.srt, features = gene, reduction = 'tsne',pt.size = 0.1, ncol = 5)
FeaturePlot(obj.srt, features = rownames(obj.srt)[(grep('Il', rownames(obj.srt)))], 
            reduction = 'tsne',pt.size = 0.1, ncol = 8)

group = 'RNA_snn_res.0.4'
group2 = 'tag'
gene <- c('Cd28', "Pik3ca",'Akt1','Cdk1','Cdk2','Cdk4','Cdk9', 'Cdkn1b')
VlnPlot(obj.srt, group.by = group, features = gene, pt.size = 0.1, ncol = 4, split.by = group2)


ggplot(data.frame(magic.data), aes(Il2ra, Il6ra, color=meta$tag)) + 
  geom_point(size=0.2) +facet_grid(meta$RNA_snn_res.0.1~meta$tag) +ggpubr::stat_cor(method = 'pearson')


ggplot(meta, aes(tag,fill=Phase)) + geom_bar(position = 'fill')




#########################################################################
##################  subset naive and wt ###############



obj <- readRDS('rds/2021.cell_hashing/2021.10.22.all.rds')

obj.tmp <- subset(obj, tag %in% c('naive','wt'))
obj.tmp %>% dim()


obj.tmp <- FindVariableFeatures(obj.tmp, selection.method = 'vst', nfeatures = 3000)
all.genes <- rownames(obj.tmp)
obj.tmp <- ScaleData(obj.tmp, features = all.genes)
obj.tmp <- NormalizeData(obj.tmp)
obj.tmp <- RunPCA(obj.tmp, features = VariableFeatures(object = obj.tmp))
ElbowPlot(obj.tmp)
DimPlot(obj.tmp, reduction = "pca", group.by = 'tag')

obj.tmp <- FindNeighbors(obj.tmp, dims = 1:10)
obj.tmp <- FindClusters(obj.tmp, resolution = c(0.1,0.2,0.3,0.4,0.5,1,1.5))
obj.tmp <- RunTSNE(obj.tmp, dims = 1:10)
obj.tmp <- RunUMAP(obj.tmp, dims = 1:10)
saveRDS(obj.tmp, 'rds/2021.cell_hashing/2021.11.01.naive_wt.rds')

p1 <- DimPlot(obj.tmp, reduction = c('pca','tsne','umap')[2], group.by = 'tag')
DimPlot(obj.tmp, reduction = 'tsne', group.by = 'RNA_snn_res.0.3')
DimPlot(obj.tmp, reduction = 'tsne', group.by = 'Phase')


clusters <- colnames(obj.tmp@meta.data[c(12,17,18,19,16,13,14)])
f <- lapply(clusters, function(cluster){
  p <- DimPlot(obj.tmp, reduction = 'tsne', group.by = cluster, pt.size = 0.2)
  print(p)
})
cowplot::plot_grid(p1, f[[1]],f[[2]],f[[3]],f[[4]],f[[5]],f[[6]],f[[7]], ncol = 3)



genes <- c('Tcf7', 'Ctla4','Bcl6','Icos','Pdcd1','Cxcr5','Slamf1')
FeaturePlot(obj.tmp, features = genes, reduction = 'tsne', ncol = 4, pt.size = 0.1)
VlnPlot(obj.tmp, features = genes, group.by = 'tag',ncol = 4, pt.size = 0.1)


clusters <- colnames(obj.tmp@meta.data[c(12,17,18,19,16,13,14)])
clusters
i <- 1
Idents(obj.tmp) <- clusters[i]
markers <- FindAllMarkers(obj.tmp, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)
write.csv(markers, paste0('rds/2021.cell_hashing/res/2021.11.01.naive_wt_markers.res.', clusters[i], '.csv'))

for(i in seq_along(clusters)){
  Idents(obj.tmp) <- clusters[i]
  markers <- FindAllMarkers(obj.tmp, only.pos = TRUE, min.pct = 0.25, 
                            logfc.threshold = 0.25)
  write.csv(markers, paste0('rds/2021.cell_hashing/res/2021.11.01.naive_wt_markers.res.', clusters[i], '.csv'))
}

top_genes <- markers %>% group_by(cluster) %>% top_n(20, wt=avg_log2FC)
DoHeatmap(obj.tmp, features = top_genes$gene) + NoLegend() + 
  theme(text = element_text(size=6)) 


files <- list.files('rds/2021.cell_hashing/res', full.names = T)
files
i <- 1
substr(files,55,nchar(files)-4)


for(i in 1:length(files)){
  res <- read.csv(files[i], header = T)
  assign(substr(files[i],55,nchar(files[i])-4), res)
}
ls(pattern = 'res.RNA')


for(i in 1:length(files)){
  markers <- read.csv(files[i], header = T)
  p <- vector()
  p0 <- print(substr(files[i],55,nchar(files[i])-4))
  p1 <- markers %>% filter(p_val_adj <=0.05) %>% dim()
  p2 <- markers %>% filter(p_val_adj <=0.05) %>% filter(avg_log2FC >=1 | avg_log2FC <=-1) %>% dim()
  p3 <- markers %>% filter(p_val_adj <=0.05) %>% filter(avg_log2FC >=1.5 | avg_log2FC <=-1.5) %>% dim()
  p4 <-markers %>% filter(p_val_adj <=0.05) %>% filter(avg_log2FC >=2 | avg_log2FC <=-2) %>% dim()
  assign(substr(files[i],55,nchar(files[i])-4), markers)
  p <- cbind(p,p0,p1,p2,p3,p4)
  print(p)
}

ls(pattern = 'res.RNA')

meta <- obj.tmp@meta.data %>% data.frame()
clusters <- colnames(obj.tmp@meta.data[c(12,17,18,19,16,13,14)])

f <- list()
f <- lapply(clusters, function(cluster){
  p <- ggplot(meta, aes_string(cluster)) + geom_bar(aes(fill=meta$tag)) + RotatedAxis()
  print(p)
})
cowplot::plot_grid(f[[1]],f[[2]],f[[3]],f[[4]],
                   f[[5]],f[[6]],f[[7]], ncol = 3)


files <- list.files('rds/2021.cell_hashing/res', full.names = T)
files
for(i in 1:length(files)){
  markers <- read.csv(files[i], header = T)
  markers.fdr <- markers %>% filter(p_val_adj <=0.05) 
  markers.res.1 <- markers %>% filter(p_val_adj <=0.05) %>% filter(avg_log2FC >=1 | avg_log2FC <=-1)
  markers.res.1.5 <- markers %>% filter(p_val_adj <=0.05) %>% filter(avg_log2FC >=1.5 | avg_log2FC <=-1.5)
  markers.res.2 <- markers %>% filter(p_val_adj <=0.05) %>% filter(avg_log2FC >=2 | avg_log2FC <=-2)
  write.csv(markers.fdr, paste0('rds/2021.cell_hashing/res/',substr(files[i],55,nchar(files)-4), '.fdr.csv'))
  write.csv(markers.res.1, paste0('rds/2021.cell_hashing/res/',substr(files[i],55,nchar(files)-4), '.res.1.csv'))
  write.csv(markers.res.1.5, paste0('rds/2021.cell_hashing/res/',substr(files[i],55,nchar(files)-4), '.res.1.5.csv'))
  write.csv(markers.res.2, paste0('rds/2021.cell_hashing/res/',substr(files[i],55,nchar(files)-4), '.res.2.csv'))
}



cellcycle.genes <- readRDS('rds/2021.cell_hashing/cell_cycle.genes.rds')
obj.tmp <- CellCycleScoring(obj.tmp, g2m.features = cellcycle.genes[[1]],
                            s.features = cellcycle.genes[[2]])
