#!/bin/sh
#SBATCH -n 16
#SBATCH -A seyedam_lab
#SBATCH --output=210411_flair.out
#SBATCH --error=210411_flair.err
#SBATCH --time=7-00:00:00
#SBATCH -J flair
#SBATCH --mem=64G
#SBATCH --mail-type=START,END
#SBATCH --mail-user=freese@uci.edu

# source activate FLAIR
module load samtools
# module load bedtools
module load minimap2


# Don't time this step because it's just the mapping
# # flair align
# genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa
# reads=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/unmapped_reads/flnc.fastq

f_dir=/dfs6/pub/freese/mortazavi_lab/bin/flair/

# python ${f_dir}flair.py align -g $genome \
#                               -r $reads \
#                               -t 16 \
#                               -o r1_flair.aligned \
#                               -v1.3 
                              
# reads=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/unmapped_reads/flnc.fastq

# python ${f_dir}flair.py align -g $genome \
#                               -r $reads \
#                               -t 16 \
#                               -o r2_flair.aligned \
#                               -v1.3 

# # Don't time this step because it's equivalent to TranscriptClean
# # flair correct
# gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
# query=r1_flair.aligned.bed
# chromSizes=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.chromSizes
# genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

# python ${f_dir}flair.py correct -f $gtf \
#                               -q $query \
#                               -t 16 \
#                               -o r1_flair \
#                               -c $chromSizes \
#                               -g $genome
                              
# gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
# query=r2_flair.aligned.bed
# chromSizes=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.chromSizes
# genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

# python ${f_dir}flair.py correct -f $gtf \
#                               -q $query \
#                               -t 16 \
#                               -o r2_flair \
#                               -c $chromSizes \
#                               -g $genome

# flair collapse
cat r1_flair_all_corrected.psl r2_flair_all_corrected.psl > r1r2_flair_all_corrected.psl
cat /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/unmapped_reads/flnc.fastq \
    /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/unmapped_reads/flnc.fastq \
    > r1r2_concat.fastq
    
reads=r1r2_concat.fastq
gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
query=r1r2_flair_all_corrected.psl
genome=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/hg38_SIRV/hg38_SIRV.fa

python ${f_dir}flair.py collapse -f $gtf \
                              -q $query \
                              -t 16 \
                              -g $genome \
                              -o r1r2_flair.collapse \
                              -r $reads
                              
# make config file
# printf "GM12878_Rep1\tGM12878\tbatch1    /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/unmapped_reads/flnc.fastq
# GM12878_Rep2GM12878\tbatch1\t/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/unmapped_reads/flnc.fastq" >> 210411_config.tsv

# run flair quantify
mkdir -p 210411_tmp
config=210411_config.tsv
isoforms=r1r2_flair.collapse.isoforms.fa

python ${f_dir}flair.py quantify -i $isoforms \
                              -t 16 \
                              -r $config \
                              --tpm \
                              --temp_dir 210411_tmp \
                              -o 210411_counts.tsv