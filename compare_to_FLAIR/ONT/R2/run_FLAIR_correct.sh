#!/bin/bash
#$ -q sam
#$ -pe one-node-mpi 16
#$ -R y
#$ -N ONT_GM12878_2_flair_correct
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

source activate FLAIR
module load bedtools

gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
query=flair.aligned.bed
chromSizes=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.chromSizes
genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

python ../../flair-1.4/flair.py correct -f $gtf \
                              -q $query \
                              -t 16 \
                              -c $chromSizes \
                              -g $genome

