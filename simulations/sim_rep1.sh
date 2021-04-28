#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -N sim_1

conda activate nanosim
ns_path=/data/users/freese/mortazavi_lab/bin/NanoSim/src/
model=mortazavi_model/ont_training
ref_genome=~/mortazavi_lab/ref/hg38/hg38.fa
ref_transcriptome=gencode.v29.transcripts.fa 

python ${ns_path}simulator.py transcriptome \
	-rt $ref_transcriptome \
	-e expression_abundance.tsv \
	-c $model \
	-o rep1/ont \
	-n 2000000 \
	-b guppy \
	-r dRNA \
	--no_model_ir \
	--uracil \
	-t 16
