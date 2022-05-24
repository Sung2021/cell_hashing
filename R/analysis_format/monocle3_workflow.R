library(monocle3)
library(SingleCellExperiment)
## 2. prepare input data
obj.srt # seurat object
data <- as(as.matrix(obj.srt@assays$RNA@data), 'sparseMatrix')

## preparing pd,fData,fd from seurat object
cell_metadata <- new('AnnotatedDataFrame', 
                 data = obj.srt@meta.data[,c(2:7,20)])
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
gene_metadata <- new('AnnotatedDataFrame', data = fData)
## convert the AnnotatedDataFrame to data.frame
gene_metadata <- as(gene_metadata, "data.frame")
cell_metadata <- as(cell_metadata, "data.frame")
## monocle set
cds <- new_cell_data_set(data, 
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_metadata)
