# Supplement

Files/paths:
```bash
PLOTPATH=../plotting_scripts
DATA=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON

db=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/full_gencode_v29_SIRV_2020-02-29.db

abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/ont_talon_abundance.tsv
filt_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/ont_talon_abundance_filtered.tsv

kallisto1=../Illumina/GM12878/Kallisto/Rep1/abundance.tsv
kallisto2=../Illumina/GM12878/Kallisto/Rep2/abundance.tsv
```

## Figure S3 and Figure S7a-b: TALON read length distributions for PacBio and ONT GM12878 datasets

Remove SIRV and ERCC reads
```bash
awk '{if($4 !~ "ERCC" && $4 !~ "SIRV") print $0}' \
    $DATA/GM12878_talon_read_annot.tsv > read_lengths/GM12878_main_read_annot.tsv
```
Plot
```bash
python ${PLOTPATH}/plot_read_length_distributions.py \
    --r read_lengths/GM12878_main_read_annot.tsv \
     --datasets PB_GM12878_R1,PB_GM12878_R2,ONT_GM12878_R1,ONT_GM12878_R2 \
     --map read_lengths/read_length_name_mapping.csv \
     --o read_lengths/
```
See resulting plots [here](https://github.com/dewyman/TALON-paper-2020/tree/master/Supplement/read_lengths).


# Figure S5B: Number of exons per transcript model detected in PacBio GM12878 transcriptomes. Transcripts are grouped by novelty type assignment.
```bash
pb_filt_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_talon_abundance_filtered.tsv
Rscript ${PLOTPATH}/plot_n_exons_by_novelty.R \
    --f ${pb_filt_abundance} \
    -o figures/
```
<img align="center" width="600" src="figures/transcript_exonCount_by_novelty_type.png">

# Figure S7: Transcript and gene quantification by PacBio/ONT and TALON 
```bash
abundance=${DATA}/pb_ont_talon_abundance.tsv
filt_abundance=${DATA}/pb_ont_talon_abundance_filtered.tsv
```

## Panel C: Expression level of known and ISM transcript models in PacBio/ONT in GM12878
```bash
Rscript ${PLOTPATH}/plot_longread_transcript_expression_corr.R \
         --f ${filt_abundance} \
         --d1 PB_GM12878_R1 \
         --d1_type 'PacBio' \
         --d2 ONT_GM12878_R1 \
         --d2_type 'ONT' \
         --celltype GM12878 \
         --ISM \
         -o figures/
```
<img align="center" width="400" src="figures/PB_GM12878_R1-ONT_GM12878_R1_Known-ISM_transcript_correlationPlot.png">
Correlations are in PB_GM12878_R1-ONT_GM12878_R1_Known-ISM_transcript_correlations.txt. 

## Figure S9:Known gene expression correlation between long read platforms and Illumina
### PacBio
```bash
mkdir -p lr_sr_corr
abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_talon_abundance.tsv
Rscript ${PLOTPATH}/plot_longread_illumina_gene_correlation.R \
    --f ${abundance} \
    --ik1 ${kallisto1} \
    --ik2 ${kallisto2} \
    --color blue \
    --r1 PB_GM12878_R1 \
    --r2 PB_GM12878_R2 \
    --srtype Illumina \
    --lrtype PacBio \
    -o lr_sr_corr
```
### ONT
```bash
mkdir -p lr_sr_corr
abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/ont_talon_abundance.tsv
Rscript ${PLOTPATH}/plot_longread_illumina_gene_correlation.R \
    --f ${abundance} \
    --ik1 ${kallisto1} \
    --ik2 ${kallisto2} \
    --color blue \
    --r1 ONT_GM12878_R1 \
    --r2 ONT_GM12878_R2 \
    --srtype Illumina \
    --lrtype ONT \
    -o lr_sr_corr
```
