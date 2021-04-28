infile=$1
outfile=$2

module load samtools
samtools view -H $infile > $outfile
grep -v “@” ont_nanosim_mapped.sam | awk ‘{if ($2 ==0 || $2 == 1 || $2 == 16) print $0}' >> $outfile