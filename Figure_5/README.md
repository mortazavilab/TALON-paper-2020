# Figure 5
A step-by-step description of how we generated the panels of Figure 5 in the TALON manuscript

First, some filepaths:
```
PLOTPATH=../plotting_scripts
abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/mouse_TALON/Brain_talon_abundance.tsv
filt_abundance=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/mouse_TALON/Brain_talon_abundance_filtered.tsv
```
Software versions:
* R v3.6.0
* edgeR v3.28.1 

## Panel A: Total number of PacBio reads assigned to each novelty category after transcript filtering in Cortex
```
Rscript ${PLOTPATH}/plot_novelty_category_read_counts.R \
         --f ${filt_abundance}  \
         --dataset PacBio_Cortex_Rep1,PacBio_Cortex_Rep2 \
         --o plots/

```
<img align="center" width="400" src="plots/PacBio_Cortex_Rep1-PacBio_Cortex_Rep2_reads_by_isoform_category.png">

## Panel B: Total number of PacBio reads assigned to each novelty category after transcript filtering in Hippocampus
```
Rscript ${PLOTPATH}/plot_novelty_category_read_counts.R \
         --f ${filt_abundance}  \
         --dataset PacBio_Hippocampus_Rep1,PacBio_Hippocampus_Rep2 \
         --o plots/
```
<img align="center" width="400" src="plots/PacBio_Hippocampus_Rep1-PacBio_Hippocampus_Rep2_reads_by_isoform_category.png">	

## Panel C: Differential Transcript Expression
Check the Rscript 

<img align="center" src="plots/volcanoPlot_novelty_labels.png">

## Panel D: 
```
cp S30.py and TALONClass to the directory where ${filt_abundance} is then run:

python S30.py
 
```
<img align="center" src="plots/cxhcNoveltyVenn.png">

## Panel E: Psnir UCSC Genome Browser screenshot and expression
Colored tracks for UCSC Genome browser where generated using the script located in
```
../analysis_scripts/gen_novelty_tracks_gtf.py
```

<img align="center" src="plots/Psnir.pdf">

To plot the expression of each isoform see the Rscript Barplot_pnisr.R

<img align="center" src="plots/Pnisr_expression.png">
