# Count known genes and transcripts from a TALON read annot file
# Genomic transcripts are not included
set -e 

read_annot=$1
dataset=$2

# Check header first
if [[  $(head -1 $read_annot | awk '{print $4}') != "chrom" ]]; then
    echo "Error: 4th column of file must be 'chrom'"
    exit  
fi
if [[  $(head -1 $read_annot | awk '{print $10}') != "gene_ID" ]]; then
    echo "Error: 10th column of file must be 'gene_ID'"
    exit
fi
if [[  $(head -1 $read_annot | awk '{print $11}') != "transcript_ID" ]]; then
    echo "Error: 11th column of file must be 'transcript_ID'"
    exit
fi
if [[  $(head -1 $read_annot | awk '{print $16}') != "gene_novelty" ]]; then
    echo "Error: 16th column of file must be 'gene_novelty'"
    exit
fi
if [[  $(head -1 $read_annot | awk '{print $17}') != "transcript_novelty" ]]; then
    echo "Error: 17th column of file must be 'transcript_novelty'"
    exit
fi


known_genes=$(grep -v "Genomic" $read_annot | \
              grep $dataset | \
              grep -v "gSpikein_ERCC" | \
              awk '{if($4 !~ "SIRV" && $16 == "Known") print $10}' | sort -u | wc -l)
known_transcripts=$(grep $dataset $read_annot | \
                    grep -v "gSpikein_ERCC" | \
                    awk '{if($4 !~ "SIRV" && $17 == "Known") print $11}' | sort -u | wc -l)


echo "n Known genes in $dataset: $known_genes"
echo "n Known transcripts in $dataset: $known_transcripts"
