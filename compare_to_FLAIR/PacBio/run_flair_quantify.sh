#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N PB_flair_q
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

source activate FLAIR
module load samtools
module load bedtools

mkdir -p tmp
config=config.tsv
isoforms=flair.collapse.isoforms.fa

python ~/flair/flair.py quantify -i $isoforms \
                              -t 4 \
                              -m ~/minimap2-2.17 \
                              -r $config \
                              --tpm \
                              --temp_dir tmp
