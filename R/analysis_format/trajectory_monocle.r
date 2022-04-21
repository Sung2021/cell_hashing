## trajectory pipeline using monocle
## 1. load library
library(scater)
library(SingleCellExperiment)
library(monocle)


############################################################
## 2. prepare input data
obj.srt <- input.data # seurat object
data <- as(as.matrix(obj.srt@assays$RNA@data), 'sparseMatrix')

## preparing pd,fData,fd from seurat object
pd <- new('AnnotatedDataFrame', data = obj.srt@meta.data[,c(1:6)])
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

## 3. run monocle algorithm
# Run ordering algorithm
var_genes <- obj.srt@assays$RNA@var.features
ordering_genes <- var_genes

## create monocle object
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

## order cells
HSMM <- orderCells(HSMM)

group='time'
group='RNA_snn_res.1'
group='time_cluster'

## color by obj.srt meta data

plot_cell_trajectory(HSMM, 
                     color_by = obj.srt@meta.data[,group],
                     theta = -15,
                     show_branch_points = F,
                     show_tree = T,
                     cell_size = 1) +
  theme(legend.position = "right") 

## save trajectory dimension to data frame

ddtree <- data.frame(t(HSMM@reducedDimS))
ddtree[1:3,]
dim(ddtree)

## plot it
ggplot(ddtree, aes(X1,X2, color=obj.srt@meta.data[,group])) + geom_point(size=0.5, alpha=0.5) +
  facet_wrap(.~obj.srt@meta.data[,group],ncol=7) + theme_classic()

pData(HSMM)

## save monocle output to seurat meta data
obj.srt[['monocle_Pseudotime']] <- pData(HSMM)$Pseudotime
obj.srt[['monocle_State']] <- pData(HSMM)$State



######################################################################
###### 2022.04.21 #####
library(scater)
library(SingleCellExperiment)
library(monocle)


############################################################
## 2. prepare input data
# seurat object
data <- as(as.matrix(obj.srt@assays$RNA@data), 'sparseMatrix')

## preparing pd,fData,fd from seurat object
pd <- new('AnnotatedDataFrame', data = obj.srt@meta.data[,c(1:6)])
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

## 3. run monocle algorithm
# Run ordering algorithm
var_genes <- obj.srt@assays$RNA@var.features
ordering_genes <- var_genes

## create monocle object
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

## order cells
HSMM <- orderCells(HSMM)

group='time'
group='RNA_snn_res.1'
group='time_cluster'

## color by obj.srt meta data

plot_cell_trajectory(HSMM, 
                     color_by = obj.srt@meta.data[,group],
                     theta = -15,
                     show_branch_points = F,
                     show_tree = T,
                     cell_size = 1) +
  theme(legend.position = "right") 

## save trajectory dimension to data frame

ddtree <- data.frame(t(HSMM@reducedDimS))
ddtree[1:3,]
dim(ddtree)

HSMM %>% saveRDS('2022.cell_hashing/imm_timecourse/imm_timecourse_monocle.22.04.21.rds')
