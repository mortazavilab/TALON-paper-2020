#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -N flair_q_perf
#$ -M freese@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools
module load minimap2/2.17

mkdir -p tmp_perf
config=perf_config.tsv
isoforms=perf_collapse.isoforms.fa
flair_dir=/dfs3/pub/freese/mortazavi_lab/bin/flair/

python ${flair_dir}flair.py quantify -i $isoforms \
                              -t 4 \
                              -r $config \
                              --tpm \
                              --temp_dir tmp_perf \
                              -o perf_counts.tsv