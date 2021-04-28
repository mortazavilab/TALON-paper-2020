#!/bin/bash
#$ -q som,bio,free64
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -ckpt restart
#$ -N sirv_1_correct

module load bedtools
module load samtools 
cd /data/users/freese/mortazavi_lab/bin/TALON-paper-2020/compare_to_FLAIR/SIRV/r1/

gtf=/data/users/freese/mortazavi_lab/bin/TALON-paper-2020/compare_to_FLAIR/SIRV/sirv.gtf
query=flair.aligned.bed
chromSizes=/data/users/freese/mortazavi_lab/bin/TALON-paper-2020/compare_to_FLAIR/SIRV/sirv_chrom_sizes.tsv
genome=/data/users/freese/mortazavi_lab/bin/TALON-paper-2020/compare_to_FLAIR/SIRV/sirv.fa

python ~/mortazavi_lab/bin/flair/flair.py correct -f $gtf \
                              -q $query \
                              -t 4 \
                              -c $chromSizes \
                              -g $genome