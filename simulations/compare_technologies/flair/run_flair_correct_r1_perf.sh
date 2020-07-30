#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N flair_correct_r1_p
#$ -M freese@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load bedtools

flair_dir=/dfs3/pub/freese/mortazavi_lab/bin/flair/

gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
query=r1_perf_aligned.bed
chromSizes=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.chromSizes
genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

python ${flair_dir}flair.py correct -f $gtf \
                              -q $query \
                              -t 16 \
                              -c $chromSizes \
                              -g $genome \
                              -o r1_perf