# Download RNA-PET data from the ENCODE portal. Mapped to hg19
# Modify file to remove the score field because it crashes liftOver

mkdir -p data
cd data

# GM12878: clone-free
wget https://www.encodeproject.org/files/ENCFF001TIL/@@download/ENCFF001TIL.bed.gz
gunzip ENCFF001TIL.bed.gz
awk -v OFS="\t" '{print $1,$2,$3,$4,0,$6}' ENCFF001TIL.bed > GM12878.bed
rm ENCFF001TIL.bed

