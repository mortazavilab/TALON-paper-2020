#!/bin/bash
#$ -q sam128
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N kallisto_GM12878
#$ -pe one-node-mpi 16
#$ -R y


module load kallisto/0.43.1
mkdir -p Kallisto/Rep2

illumina1=GM12878_rep2_illumina_1.fastq.gz
illumina2=GM12878_rep2_illumina_2.fastq.gz
transcriptomeFasta=../../refs/gencode.v29.transcripts.fa.gz
indexName=../../refs/Kallisto/gencode_v29.idx

kallisto quant -i $indexName -t 16 -o Kallisto/Rep2 -b 100 $illumina1 $illumina2
