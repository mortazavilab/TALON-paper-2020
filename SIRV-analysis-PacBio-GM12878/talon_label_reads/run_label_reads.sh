#$ -q sam
#$ -pe one-node-mpi 4
#$ -R y
#$ -N LR-SIRV
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

set -e
module load dwyman/anaconda/3
source activate python3.6

talon_label_reads --f /share/crsp/lab/seyedam/dwyman/TALON-paper-2020/SIRV-analysis-PacBio-GM12878/TC_SIRV_general_minimap/Rep1/TC_clean.sam \
    --g ../../refs/hg38_SIRV/hg38_SIRV.fa  \
    --t 8 \
    --ar 10 \
    --deleteTmp \
    --o Rep1

talon_label_reads --f /share/crsp/lab/seyedam/dwyman/TALON-paper-2020/SIRV-analysis-PacBio-GM12878/TC_SIRV_general_minimap/Rep2/TC_clean.sam \
    --g ../../refs/hg38_SIRV/hg38_SIRV.fa  \
    --t 8 \
    --ar 10 \
    --deleteTmp \
    --o Rep2


source deactivate
