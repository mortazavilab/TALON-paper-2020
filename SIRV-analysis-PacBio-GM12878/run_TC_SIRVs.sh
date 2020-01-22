#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N SIRV-TC-dry
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load dwyman/anaconda/3
source activate mypython3.7.2

mkdir -p TC_SIRV/Rep1
mkdir -p TC_SIRV/Rep2
mkdir -p TC_SIRV/both_reps

cd TC_SIRV/Rep1
time python /data/users/dwyman/TranscriptClean-2.0.2/TranscriptClean.py \
    --sam ../../PacBio_Sequel2_GM12878_R1_SIRV_reads.sam \
    --genome ../refs/SIRV_ERCC/SIRV.fa \
    --spliceJns ../refs/SIRV_ERCC/SIRV_SJs.tsv  \
    -t 16 \
    --canonOnly \
    --deleteTmp \
    --outprefix TC
time Rscript /data/users/dwyman/TranscriptClean-2.0.2/generate_report.R TC

cd ../Rep2
time python /data/users/dwyman/TranscriptClean-2.0.2/TranscriptClean.py \
    --sam ../../PacBio_Sequel2_GM12878_R2_SIRV_reads.sam \
    --genome ../refs/SIRV_ERCC/SIRV.fa \
    --spliceJns ../refs/SIRV_ERCC/SIRV_SJs.tsv  \
    -t 16 \
    --canonOnly \
    --deleteTmp \
    --outprefix TC
time Rscript /data/users/dwyman/TranscriptClean-2.0.2/generate_report.R TC
