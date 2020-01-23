# Time consuming- run in screen

set -e
mkdir -p TranscriptClean
cd TranscriptClean

# Paths specific to Dana's system.
TC_dir=~/TranscriptClean-2.0.2
module load dwyman/anaconda/3
source activate mypython3.7.2

# For human
python $TC_dir/accessory_scripts/get_SJs_from_gtf.py \
      --f ../GENCODE_v29/gencode.v29.annotation.gtf \
      --g ../hg38/hg38.fa \
      --o gencode_v29_SJs.tsv
