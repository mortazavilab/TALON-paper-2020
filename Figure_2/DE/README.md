# Differential Expression

Here, we investigate whether factors like gene length or GC content play a role in the different gene/transcript expression measurements we get from PacBio and Illumina.

The starting point is an output file that was generated from the code that made the MA plots.
```
PLOTPATH=../../plotting_scripts
mkdir -p plots
python $PLOTPATH/plot_gene_length_by_DE.py \
    --f ../plots/edgeR_PacBio_illumina_gene_counts.tsv \
    --col median_length \
    --platform PacBio \
    --o plots/
```
