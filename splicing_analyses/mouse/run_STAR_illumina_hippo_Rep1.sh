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
#$ -N hipp_1

DPATH=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/
curpath=~/mortazavi_lab/bin/TALON-paper-2020/splicing_analyses/mouse/
reads="${DPATH}illumina_Hippocampus/Rep1/389_S9_R1.fastq.gz ${DPATH}illumina_Hippocampus/Rep1/389_S9_R2.fastq.gz"

module load STAR/2.5.2a
STAR \
	--runThreadN 8 \
	--genomeDir /data/users/freese/mortazavi_lab/ref/mm10/STAR_mm10_ENCODE/ \
	--readFilesIn $reads \
	--sjdbGTFfile /data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outFileNamePrefix ${curpath}hippo_1_aligned \
	--outSAMattributes NH HI NM MD jM jI \
	--outSAMtype SAM \
	--readFilesCommand zcat
