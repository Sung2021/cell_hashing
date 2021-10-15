#!/bin/bash

path=~/Desktop/Sung_work/fastq/cell_hashing/19092-44
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar FastqToSam \
F1= $path/01/all.R1.fastq.gz \ 
F2= $path/01/all.R2.fastq.gz \
O= $path/temp/all.unsorted.bam SM=CH_transcriptome

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagBamWithReadSequenceExtended \
INPUT= $path/temp/all.unsorted.bam \
OUTPUT= $path/temp/all.celltag.bam \
SUMMARY= $path/temp/all.celltag.summary.txt \
BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagBamWithReadSequenceExtended \
INPUT= $path/temp/all.celltag.bam \
OUTPUT= $path/temp/all.cellmoleculartag.bam \
SUMMARY= $path/temp/all.cellmoleculartag.summary.txt \
BASE_RANGE=17-28 BASE_QUALITY=10 \
BARCODED_READ=1 DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1

### FilterBam
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/FilterBam \
TAG_REJECT=XQ \
INPUT= $path/temp/all.cellmoleculartag.bam \
OUTPUT= $path/temp/all.cellmoleculartag.filtered.bam 

### TrimStartingSequence ## error on 10X 
#/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TrimStartingSequence \
#INPUT= $path/temp/all.cellmoleculartag.filtered.bam \
#OUTPUT= $path/temp/all.cellmoleculartag.trimmed.bam \ 
#OUTPUT_SUMMARY= $path/temp/all.cellmoleculartag.trimming_report.txt \
#SEQUENCE= AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 

### PolyATrimmer 
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/PolyATrimmer INPUT= $path/temp/all.cellmoleculartag.filtered.bam OUTPUT= $path/temp/all.cellmoleculartag.polyA.trimmed.bam OUTPUT_SUMMARY= $path/temp/all.cellmoleculartag.polyA.trimming_report.txt MISMATCHES=0 NUM_BASES=6 NEW=true


### SamToFastq
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar SamToFastq INPUT=$path/temp/all.cellmoleculartag.polyA.trimmed.bam FASTQ=$path/temp/all.fastq


### STAR alignment
### STAR only support uncompressed fastq
### gzip -d fastq.gz # if it is compressed
~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 8 \
--readFilesIn $path/temp/all.fastq \
--genomeDir ~/Desktop/ref/mm10/star_index \
--outFileNamePrefix $path/temp/all.star\
--outSAMtype BAM SortedByCoordinate

### SortSam
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar SortSam I= $path/temp/all.star_--outSAMtypeAligned.out.sam O=$path/temp/all.sorted.bam SO=queryname

### MergeBamAlignment
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=~/Desktop/ref/mm10/mm10.fa UNMAPPED_BAM= $path/temp/all.cellmoleculartag.polyA.trimmed.bam ALIGNED_BAM= $path/temp/all.sorted.bam OUTPUT=$path/temp/merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false


### TagReadWithGeneExon
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagReadWithGeneFunction I=$path/temp/merged.bam O=$path/temp/star_gene_exon_tagged.bam ANNOTATIONS_FILE=~/Desktop/ref/mm10/mm10.refFlat

### DigitalExpression
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/DigitalExpression I=$path/temp/star_gene_exon_tagged.bam O=$path/temp/star_gene_exon_tagged.dge.txt.gz SUMMARY=$path/temp/dge.summary.txt NUM_CORE_BARCODES=10000 LOCUS_FUNCTION_LIST=null STRAND_STRATEGY=null
