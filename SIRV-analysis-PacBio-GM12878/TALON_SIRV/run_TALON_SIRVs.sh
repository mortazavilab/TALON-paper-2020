#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N TALON-SIRV
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

sewt -e
module load dwyman/anaconda/3
source activate TALON-v4.4.2

DATE=$(date +%F\(%a\) | awk -F"(" '{print $1}')
time talon_initialize_database \
    --f /share/crsp/lab/seyedam/share/PACBIO/genomes/SIRV/SIRV-spikein.gtf \
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

time talon_abundance \
    --db SIRV_${DATE}.db \
    -a SIRV \
    -b SIRV \
    --o SIRV

source deactivate
