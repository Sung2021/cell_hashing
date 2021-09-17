#!/bin/bash
## for a single paired end (R1, R2) files
path=~/Desktop//Sung_work/fastq/cell_hashing/19092-44
## fastqc
## trimming 

CITE-seq-Count -R1 $path/19092FL-44-01-01-01_S64_L007_R1_001.fastq.gz -R2 $path/19092FL-44-01-01-01_S64_L007_R2_001.fastq.gz \
-t ~/Desktop//Sung_work/fastq/cell_hashing/tag_list.csv \
-cbf 1 -cbl 16 -umif 17 -umil 26 -cells 2000 \
-o ~/Desktop//Sung_work/fastq/cell_hashing/output/
