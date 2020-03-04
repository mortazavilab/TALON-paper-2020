# Figure 3
A step-by-step description of how we generated the panels of Figure 3 in the TALON manuscript

Firstbash, some filepaths:
```bash
PLOTPATH=../plotting_scripts
abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/ont_talon_abundance.tsv
filt_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/ont_talon_abundance_filtered.tsv
pb_ont_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_ont_talon_abundance.tsv
pb_ont_filt_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_ont_talon_abundance_filtered.tsv
kallisto1=../Illumina/GM12878/Kallisto/Rep1/abundance.tsv
kallisto2=../Illumina/GM12878/Kallisto/Rep2/abundance.tsv
mkdir plots
```
Software versions:
* R v3.5.1  

## Panel A: Expression level of known genes (GENCODE v29) in each biological replicate of GM12878
```R
Rscript ${PLOTPATH}/plot_longread_gene_expression_corr.R \
          --f ${abundance} \
          --color blue \
          --d1 ONT_GM12878_R1 \
          --d2 ONT_GM12878_R2 \
          --celltype GM12878 \
          --d1_type 'Rep1 ONT' \
          --d2_type 'Rep2 ONT' \
          -o plots/
```

## Panel B
```R
Rscript ${PLOTPATH}/plot_longread_transcript_expression_corr.R  \
          --f ${filt_abundance} \
          --d1 ONT_GM12878_R1 \
          --d2 ONT_GM12878_R2 \
          --d1_type 'Rep1 ONT' \
          --d2_type 'Rep2 ONT' \
          -o plots/ 
```

## Panel C
```R
Rscript ${PLOTPATH}/plot_longread_gene_expression_corr.R \
          --f ${pb_ont_abundance} \
          --color blue \
          --d1 PB_GM12878_R1 \
          --d2 ONT_GM12878_R1 \
          --celltype GM12878 \
          --d1_type 'GM12878 PacBio' \
          --d2_type 'GM12878 ONT' \
          -o plots/
```

<!-- ```R
Rscript ${PLOTPATH}/plot_longread_gene_expression_corr.R \
          --f ${pb_ont_abundance} \
          --color blue \
          --d1 PB_GM12878_R2 \
          --d2 ONT_GM12878_R1 \
          --celltype GM12878 \
          --d1_type 'Rep2 PB' \
          --d2_type 'Rep1 ONT' \
          -o plots/
``` -->

<!-- ```R
Rscript ${PLOTPATH}/plot_longread_gene_expression_corr.R \
          --f ${pb_ont_abundance} \
          --color blue \
          --d1 PB_GM12878_R1 \
          --d2 ONT_GM12878_R2 \
          --celltype GM12878 \
          --d1_type 'Rep1 PB' \
          --d2_type 'Rep2 ONT' \
          -o plots/
```

```R
Rscript ${PLOTPATH}/plot_longread_gene_expression_corr.R \
          --f ${pb_ont_abundance} \
          --color blue \
          --d1 PB_GM12878_R2 \
          --d2 ONT_GM12878_R2 \
          --celltype GM12878 \
          --d1_type 'Rep2 PB' \
          --d2_type 'Rep2 ONT' \
          -o plots/
``` -->


## Panel D
```R
Rscript ${PLOTPATH}/plot_longread_transcript_expression_corr.R  \
          --f ${pb_ont_filt_abundance} \
          --d1 PB_GM12878_R1 \
          --d2 ONT_GM12878_R1 \
          --d1_type 'GM12878 PacBio' \
          --d2_type 'GM12878 ONT' \
          -o plots/ \
          --antisense 
```

<!-- ```R
Rscript ${PLOTPATH}/plot_longread_transcript_expression_corr.R  \
          --f ${pb_ont_filt_abundance} \
          --d1 PB_GM12878_R2 \
          --d2 ONT_GM12878_R1 \
          --d1_type 'Rep2 PB' \
          --d2_type 'Rep1 ONT' \
          -o plots/ 
```

```R
Rscript ${PLOTPATH}/plot_longread_transcript_expression_corr.R  \
          --f ${pb_ont_filt_abundance} \
          --d1 PB_GM12878_R1 \
          --d2 ONT_GM12878_R2 \
          --d1_type 'Rep1 PB' \
          --d2_type 'Rep2 ONT' \
          -o plots/ 
```

```R
Rscript ${PLOTPATH}/plot_longread_transcript_expression_corr.R  \
          --f ${pb_ont_filt_abundance} \
          --d1 PB_GM12878_R2 \
          --d2 ONT_GM12878_R2 \
          --d1_type 'Rep2 PB' \
          --d2_type 'Rep2 ONT' \
          -o plots/ 
``` -->

