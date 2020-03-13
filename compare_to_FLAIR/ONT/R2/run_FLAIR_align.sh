#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N flair_align_ONT_GM12878_2
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools

genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa
reads=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/ONT_RNA02_GM12878_R2/unmapped_reads/ONT47.fastq

python ~/flair/flair.py align -g $genome \
                              -r $reads \
                              -t 16 \
                              -m ~/minimap2-2.17 \
                              -v1.3 


