# Figure 2
A step-by-step description of how we generated the panels of Figure 2 in the TALON manuscript

First, some filepaths:
```
PLOTPATH=../plotting_scripts
abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_talon_abundance.tsv
filt_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_talon_abundance_filtered.tsv
kallisto1=../Illumina/GM12878/Kallisto/Rep1/abundance.tsv
kallisto2=../Illumina/GM12878/Kallisto/Rep2/abundance.tsv
```
Software versions:
* R v3.5.1  

## Panel A: Expression level of known genes (GENCODE v29) in each biological replicate of GM12878
```
Rscript ${PLOTPATH}/plot_longread_gene_expression_corr.R \
          --f ${abundance} \
          --color blue \
          --d1 PB_GM12878_R1 \
          --d2 PB_GM12878_R2 \
          --celltype GM12878 \
          --d1_type 'Rep1 PacBio' \
          --d2_type 'Rep2 PacBio' \
          -o plots/
```

## Panel B

## Panel C

## Panel D
