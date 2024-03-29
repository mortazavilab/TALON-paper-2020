# Liftover the RNA-PET data from hg19 to hg38
module load jje/kent
set -e

cd FANTOM5
liftOver hg19_CAGE.bed ../hg19ToHg38.over.chain hg38_CAGE.bed unlifted.bed
rm unlifted.bed
