#!/bin/bash

cat *_R1_001.fastq.gz > all.R1.fastq.gz
cat *_R2_001.fastq.gz > all.R2.fastq.gz

java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar FastqToSam F1=all.R1.fastq.gz F2=all.R2.fastq.gz O= temp/all.bam SM=ch_transcriptome

cd temp/
/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagBamWithReadSequenceExtended INPUT=all.bam OUTPUT=all.cell.bam SUMMARY=all.cell.summary.txt BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagBamWithReadSequenceExtended INPUT=all.cell.bam OUTPUT=all.cell.molecule.bam SUMMARY=all.cell.molecule.summary.txt BASE_RANGE=17-28 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/FilterBam TAG_REJECT=XQ INPUT= all.cell.molecule.bam OUTPUT= all.cell.mol.filtered.bam 

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TrimStartingSequence INPUT=all.cell.mol.filtered.bam OUTPUT=all.trimmed.bam OUTPUT_SUMMARY=all.trimmed_report.txt SEQUENCE= AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/PolyATrimmer INPUT= all.trimmed.bam OUTPUT=all.polytrimmed.bam OUTPUT_SUMMARY= all.polytrimmed.summary.txt MISMATCHES=0 NUM_BASES=6 

java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar SamToFastq INPUT=all.polytrimmed.bam FASTQ=all.fastq

ulimit -n 1000
 ~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 10 --readFilesIn all.fastq --genomeDir ~/Desktop/ref/mm10/star_index --outFileNamePrefix all.star_ --outSAMtype BAM SortedByCoordinate

ulimit -n 1000
java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar SortSam I= all.star_Aligned.sortedByCoord.out.bam O= all.star.sorted.bam SO=queryname 

java -jar ~/Desktop/software/Drop-seq_tools-2.4.1/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=~/Desktop/ref/mm10/mm10.fa UNMAPPED_BAM= all.polytrimmed.bam ALIGNED_BAM= all.star.sorted.bam OUTPUT=merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/TagReadWithGeneFunction I=merged.bam O=star_gene_exon_tagged.bam ANNOTATIONS_FILE=~/Desktop/ref/mm10/mm10.refFlat

/Users/sung/Desktop/software/Drop-seq_tools-2.4.1/DigitalExpression I=star_gene_exon_tagged.bam O=star_gene_exon_tagged.dge.txt.gz SUMMARY=dge.summary.txt NUM_CORE_BARCODES=10000 LOCUS_FUNCTION_LIST=null STRAND_STRATEGY=null


