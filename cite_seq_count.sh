#!/bin/bash

path=~/Desktop/Sung_work/fastq/cell_hashing/19092-44/01

CITE-seq-Count -R1 $path/all.R1.fastq.gz \
-R2 $path/all.R2.fastq.gz \
-t ~/Desktop/Sung_work/fastq/cell_hashing/tag_list.csv \
-cbf 1 -cbl 16 -umif 17 -umil 26 -cells 10000 \
-o $path/
