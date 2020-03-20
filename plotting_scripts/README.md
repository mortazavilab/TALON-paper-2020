# Plotting scripts
This directory contains scripts used for data analysis and visualization in the TALON manuscript. For R scripts, you can view usage instructions by running `Rscript scriptname.R --h`. For Python scripts, this can be achieved by running `python scriptname.py --h`. Here is a brief guide to the scripts:

| Script        | Description  |
| ------------- |-------------|
| [longread_v_illumina_gene_counts_edgeR.R](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/longread_v_illumina_gene_counts_edgeR.R)| Runs differential gene expression analysis on long-read vs short-read data. Generates an MA plot and DE output file. | 
| [longread_v_illumina_transcript_counts_edgeR.R](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/longread_v_illumina_transcript_counts_edgeR.R) | Runs differential transcript expression analysis on long-read vs short-read data. Generates an MA plot and DE output file. |
| [plot_GC_content_by_DE.py](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/plot_GC_content_by_DE.py) | Uses the output of a differential expression analysis to plot the GC content of genes by their DE category |
| [plot_detection_by_TPM_for_datasets.R](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/plot_detection_by_TPM_for_datasets.R) | Bins Illumina-expressed genes by TPM and then plots the fraction in each bin that was detected in one or two long-read replicates |
| [plot_frac_As_by_read_annot_category.py](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/plot_frac_As_by_read_annot_category.py) | Plot the distribution of post-alignment fraction As per TALON novelty category (exploratory plot) |
| [plot_fraction_novel_vs_total_gene_exp.py](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/plot_fraction_novel_vs_total_gene_exp.py) | Plot total read count for gene vs the fraction of reads annotated as novel for that gene (exploratory plot) |
| [plot_gene_or_transcript_length_by_DE.py](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/plot_gene_or_transcript_length_by_DE.py) | Uses the output of a differential expression analysis to plot the length of genes/transcripts by their DE category |
| [plot_longread_gene_expression_corr.R](https://github.com/dewyman/TALON-paper-2020/blob/master/plotting_scripts/plot_longread_gene_expression_corr.R) | Plot gene expression correlation for two long-read datasets |
