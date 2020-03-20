# Running FLAIR on GM12878 PacBio data

FLAIR was cloned from https://github.com/BrooksLabUCSC/flair on 8/5/2019.

1. Run align and correct steps separately on replicates
```bash
qsub R1/./run_FLAIR_align.sh
qsub R2/./run_FLAIR_align.sh
```
```bash
qsub R1/./run_FLAIR_correct.sh
qsub R2/./run_FLAIR_correct.sh
```
2. Then, run collapse step on concatenated files from both reps.
```bash
cat R1/flair_all_corrected.psl R2/flair_all_corrected.psl > R1-R2_flair_all_corrected.psl
cat /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/unmapped_reads/flnc.fastq \
    /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/unmapped_reads/flnc.fastq \
    > R1-R2-concat.fastq
qsub ./run_flair_collapse.sh
```
3. Finally, run quantify step. To do this, you need to create a tab-delimited config file with fields dataset name, condition, batch, and fastq reads file. This is what the GM12878 file looks like:
```bash
GM12878_Rep1    GM12878    batch1    /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/unmapped_reads/flnc.fastq
GM12878_Rep2    GM12878    batch1    /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/unmapped_reads/flnc.fastq
```
```bash
qsub ./run_flair_quantify.sh
```

4. In order to determine how well PacBio + FLAIR detects genes known to be expressed in short-read data, we converted the FLAIR output to a TALON-like format, and then ran a custom R script:
```bash
# Remove SIRVs and ERCCs
grep -v "gSpikein_ERCC" counts_matrix.tsv | \
        grep -v "SIRV" > counts_matrix_nospikes.tsv

python ../format_flair_matrix_like_talon.py counts_matrix_nospikes.tsv counts_matrix_talon_abd.tsv
```

5. We also want to see how reproducible our datasets are as characterized by FLAIR. Run the following gene and transcript correlations to see:
```bash
PLOTPATH=../../plotting_scripts
Rscript $PLOTPATH/plot_longread_gene_expression_corr.R \
    --f counts_matrix_talon_abd.tsv \
    --color blue \
    --d1 GM12878_Rep1_GM12878_batch1 \
    --d2 GM12878_Rep2_GM12878_batch1 \
    --celltype GM12878 \
    --d1_type "PacBio Rep1" \
    --d2_type "PacBio Rep2" \
    -o plots

Rscript $PLOTPATH/plot_longread_transcript_expression_corr.R \
   --f counts_matrix_talon_abd.tsv \
   --color blue \
    --d1 GM12878_Rep1_GM12878_batch1 \
    --d2 GM12878_Rep2_GM12878_batch1 \
    --celltype GM12878 \
    --d1_type "PacBio Rep1" \
    --d2_type "PacBio Rep2" \
    --outdir plots
```

6. Finally, let's see how the gene detection sensitivity compares between FLAIR and TALON on the same datasets. We will need the PacBio GM12878 unfiltered TALON abundance file
```bash
abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_talon_abundance.tsv

mkdir -p plots
Rscript ../compare_TALON_FLAIR_detection_to_Illumina.R \
    --talon ${abundance} \
    --flair counts_matrix_talon_abd.tsv \
    --talonD PB_GM12878_R1,PB_GM12878_R2 \
    --flairD GM12878_Rep1_GM12878_batch1,GM12878_Rep2_GM12878_batch1 \
    --ik1 ../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    -o plots
```

7. Get the number of known transcripts detected in TALON, FLAIR, and both
```
python ../compare_known_transcripts.py \
   --f counts_matrix_talon_abd.tsv \
   --t /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_talon_abundance_filtered.tsv
```
