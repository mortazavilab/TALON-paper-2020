#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -M freese@uci.edu
#$ -N st_m

printf "r1_stringtie.gtf\nr2_stringtie.gtf" > stringtie_gtfs.tsv

gtf_list=stringtie_gtfs.tsv
sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/

${sdir}./stringtie \
	--merge \
	-p 16 \
	-c 1 \
	-o merged_stringtie_no_ref.gtf \
	-l st_novel_ \
	$gtf_list