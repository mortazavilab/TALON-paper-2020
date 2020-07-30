#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N flair_r2
#$ -M freese@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools
module load minimap2/2.17

flair_dir=/dfs3/pub/freese/mortazavi_lab/bin/flair/
genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa
reads=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2/rep2.fasta

python ${flair_dir}flair.py align -g $genome \
                              -r $reads \
                              -t 16 \
                              -v1.3 \
                              -o r2_aligned 
