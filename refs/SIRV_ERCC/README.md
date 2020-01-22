# SIRVs
## Download the SIRV set that our lab uses from Lexogen
```
wget https://www.lexogen.com/wp-content/uploads/2018/08/SIRV_Set3_Sequences_170612a-ZIP.zip
unzip SIRV_Set3_Sequences_170612a-ZIP.zip
rm SIRV_Set3_Sequences_170612a-ZIP.zip
mv SIRV_Set3_Sequences_170612a\ \(ZIP\)/ SIRV_Set3_Sequences_170612a
cp SIRV_Set3_Sequences_170612a/SIRV_isoforms_multi-fasta_170612a.fasta SIRV.fa
cp SIRV_Set3_Sequences_170612a/SIRV_ERCCs_multi-fasta_170612a.fasta ERCC.fa
```
These are the SIRV/ERCC files we will be using:
* SIRV_ERCCs_multi-fasta_170612a.fasta, multi-fasta files with individual lines for each of the 92 ERCC gene sequences
* SIRV_isoforms_multi-fasta_170612a.fasta, multi-fasta files with individual lines for each of the 7 SIRV isoform gene sequences
* SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf, gtf file with correct isoform annotations for the SIRVs only

## GTF file processing
First, create an ERCC GTF by running Diane's script on the ERCC fasta (merge_encode_annotations.py found [here](https://github.com/detrout/long-rna-seq-condor/blob/master/woldrnaseq/merge_encode_annotations.py), last commit: c990527).
(requires xopen from pypi to be installed)
```
python3 merge_encode_annotations.py \
  -o ERCC-spikein.gtf \
  ERCC.fa
```

It came to light that the SIRV GTF allows transcripts from different strands to belong to the same gene, which TALON does not support. Therefore, we need to manually reconfigure the GTF file so that genes belong to only one strand.
```
cp SIRV_Set3_Sequences_170612a/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf SIRV-spikein_raw.gtf
python separate_multistrand_genes.py --f SIRV-spikein_raw.gtf --o SIRV-spikein_gene-sep.gtf
```

For TALON to work, we need the GTF file to have gene and transcript entries. So we need to run the TALON utility on SIRV-spikein_gene-sep.gtf and on ERCC-spikein_raw.gtf to get this.
```
source activate TALON-v4.4.2
talon_reformat_gtf -gtf SIRV-spikein_gene-sep.gtf
talon_reformat_gtf -gtf ERCC-spikein.gtf
```
This gives us SIRV-spikein_gene-sep_reformatted.gtf and ERCC-spikein_reformatted.gtf

## Make SIRV splice junction file
```
source activate mypython3.7.2
TC_dir=~/TranscriptClean-2.0.2
python $TC_dir/accessory_scripts/get_SJs_from_gtf.py \
      --f SIRV-spikein_gene-sep_reformatted.gtf \
      --g SIRV.fa \
      --o SIRV_SJs.tsv
```
