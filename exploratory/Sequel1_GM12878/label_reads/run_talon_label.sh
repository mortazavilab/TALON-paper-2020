#!/bin/bash
#$ -q sam
#$ -pe one-node-mpi 16
#$ -R y
#$ -N TALON-label
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load dwyman/anaconda/3
source activate python3.6

talon_label_reads \
    --f /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/data/PacBio_GM12878_1/TC_v2.0.2/TC_clean.sam \
    --g /share/crsp/lab/seyedam/dwyman/TALON-paper-2019/refs/hg38/hg38.fa \
    --t 16 \
    --o GM12878_R1_Sequel1

talon_label_reads \
    --f /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/data/PacBio_GM12878_2/TC_v2.0.2/TC_clean.sam \
    --g /share/crsp/lab/seyedam/dwyman/TALON-paper-2019/refs/hg38/hg38.fa \
    --t 16 \
    --o GM12878_R2_Sequel1

source deactivate
qstat -j $JOB_ID
