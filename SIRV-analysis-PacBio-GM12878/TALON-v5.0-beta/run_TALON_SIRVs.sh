#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N TALON-SIRV
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

set -e
module load dwyman/anaconda/3
source activate python3.6

DATE=$(date +%F\(%a\) | awk -F"(" '{print $1}')
time talon_initialize_database \
    --f ../../refs/SIRV_ERCC/SIRV-spikein_gene-sep_reformatted.gtf \
    --g SIRV \
    --a SIRV \
    --l 0 \
    --5p 500 \
    --3p 300 \
    --o SIRV_${DATE}

time talon --f config.csv \
           --db SIRV_${DATE}.db \
           --t 16 \
           --build SIRV --cov 0.9 --identity 0.8 \
           --o SIRV

# Unfiltered abundance file
time talon_abundance \
    --db SIRV_${DATE}.db \
    -a SIRV \
    -b SIRV \
    --o SIRV

# Make whitelist
talon_filter_transcripts \
    --db SIRV_${DATE}.db \
    -a SIRV \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    --maxFracA 0.5 \
    --minCount 5 \
    --minDatasets 2 \
    --o SIRV_whitelist.csv

# Filter abundance
time talon_abundance \
    --db SIRV_${DATE}.db \
    --w SIRV_whitelist.csv \
    -a SIRV \
    -b SIRV \
    --o SIRV

source deactivate
