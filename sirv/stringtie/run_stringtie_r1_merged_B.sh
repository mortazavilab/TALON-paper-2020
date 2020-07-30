#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -M freese@uci.edu
#$ -N st_r1_m

sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/
r1=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/label_reads/PacBio_Rep1_labeled.bam
r2=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/label_reads/PacBio_Rep2_labeled.bam
ref_annot=merged_stringtie.gtf

# rep 1
prefix='r1_merged_e_B'
${sdir}./stringtie \
	-L \
	-p 16 \
	-e \
	-b ${prefix}_ballgown \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv \
	$r1