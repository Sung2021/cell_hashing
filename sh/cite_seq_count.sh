#!/bin/bash

path=~/Desktop/Sung_work/fastq/cell_hashing/19092-44/01

## all.R1.fastq, all.R2.fastq : the combinded 4 fastq files  
## expected cell number : 10,000 
## when checking fastq, found that tag sequence starting after 10 bp : --start trim 10 applied

CITE-seq-Count --start-trim 10 \
-R1 $path/all.tag.R1.fastq.gz \
-R2 $path/all.tag.R2.fastq.gz \
-t ~/Desktop/Sung_work/fastq/cell_hashing/tag_list.csv \
-cbf 1 -cbl 16 -umif 17 -umil 26 -cells 10000 \
-o $path/
