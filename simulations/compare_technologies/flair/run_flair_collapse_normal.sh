#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N norm_flair_collapse
#$ -M freese@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools
module load minimap2/2.17


flair_dir=/dfs3/pub/freese/mortazavi_lab/bin/flair/
reads=normal.fasta

gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
query=normal_corrected.psl
genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

python ${flair_dir}flair.py collapse -f $gtf \
                              -q $query \
                              -t 16 \
                              -g $genome \
                              -r $reads \
                              -o norm_collapse