#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N GM12878_flair_collapse
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

source activate FLAIR
module load samtools
module load bedtools

reads=R1-R2-concat.fastq
gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
query=R1-R2_flair_all_corrected.psl
genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

python ~/flair/flair.py collapse -f $gtf \
                              -q $query \
                              -t 16 \
                              -g $genome \
                              -m ~/minimap2-2.17 \
                              -r $reads


