#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -M freese@uci.edu
#$ -N st_m

gtf_list=stringtie_gtfs.tsv
ref_annot=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/

${sdir}./stringtie \
	--merge \
	-p 16 \
	-c 1 \
	-G ${ref_annot} \
	-o merged_stringtie.gtf \
	-l st_novel_ \
	$gtf_list