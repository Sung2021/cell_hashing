#!/bin/bash

## 1. concatenate all R1 and R2 to one file each 
cat *_R1_001.fastq.gz > all.R1.fastq.gz
cat *_R2_001.fastq.gz > all.R2.fastq.gz

## 2. fastqtosam
path=~/Desktop/admera/19092-58/19092-58-02/02_01
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar FastqToSam F1=all.R1.fastq.gz F2=all.R2.fastq.gz O= temp/all.bam SM=ch_transcriptome

## 3. adding cell tag
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagBamWithReadSequenceExtended INPUT=all.bam OUTPUT=all.cell.bam SUMMARY=all.cell.summary.txt BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 

## 4. adding molecular tag
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagBamWithReadSequenceExtended INPUT=all.cell.bam OUTPUT=all.cell.molecule.bam SUMMARY=all.cell.molecule.summary.txt BASE_RANGE=17-28 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

## 5. FilterBam
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/FilterBam TAG_REJECT=XQ INPUT= all.cell.molecule.bam OUTPUT= all.cell.mol.filtered.bam 

## 6. TrimStartingSequence ## error on 10X 
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TrimStartingSequence INPUT=all.cell.mol.filtered.bam OUTPUT=all.trimmed.bam OUTPUT_SUMMARY=all.trimmed_report.txt SEQUENCE= AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5

## 7. PolyATrimmer 
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/PolyATrimmer INPUT= all.trimmed.bam OUTPUT=all.polytrimmed.bam OUTPUT_SUMMARY= all.polytrimmed.summary.txt MISMATCHES=0 NUM_BASES=6 


## 8. SamToFastq
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar SamToFastq INPUT=all.polytrimmed.bam FASTQ=all.fastq

## 9. STAR alignment
ulimit -n 1000
 ~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 10 --readFilesIn all.fastq --genomeDir ~/Desktop/ref/mm10/star_index --outFileNamePrefix all.star_ --outSAMtype BAM SortedByCoordinate

## 10. SortSam
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar SortSam I= all.star_Aligned.sortedByCoord.out.bam O= all.star.sorted.bam SO=queryname 

## 11. MergeBamAlignment
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=~/Desktop/ref/mm10/mm10.fa UNMAPPED_BAM= all.polytrimmed.bam ALIGNED_BAM= all.star.sorted.bam OUTPUT=merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

## 12. TagReadWithGeneExon
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagReadWithGeneFunction I=merged.bam O=star_gene_exon_tagged.bam ANNOTATIONS_FILE=~/Desktop/ref/mm10/mm10.refFlat

## 13. DigitalExpression
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/DigitalExpression I=star_gene_exon_tagged.bam O=star_gene_exon_tagged.dge.txt.gz SUMMARY=dge.summary.txt NUM_CORE_BARCODES=10000 LOCUS_FUNCTION_LIST=null STRAND_STRATEGY=null


## CITE-seq-Count

#!/bin/bash
path=~/Desktop/admera/19092-58/19092-58-02/02_02

cat *_R1_001.fastq.gz > all.tag.R1.fastq.gz
cat *_R2_001.fastq.gz > all.tag.R2.fastq.gz

CITE-seq-Count --start-trim 10 \
-R1 $path/all.tag.R1.fastq.gz \
-R2 $path/all.tag.R2.fastq.gz \
-t $path/tag_list.csv \
-cbf 1 -cbl 16 -umif 17 -umil 26 -cells 10000 -o .
