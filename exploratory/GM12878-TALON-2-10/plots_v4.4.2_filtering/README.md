Run filtering settings comparable to v4.4.2
```
source activate python3.6

# Make whitelist
talon_filter_transcripts \
    --db ../full_gencode_v29_SIRV_2020-02-11.db \
    -a gencode_v29_SIRV \
    --datasets GM12878_R1,GM12878_R2 \
    --maxFracA 1 \
    --minCount 1 \
    --minDatasets 2 \
    --o whitelist_v4.4.2.csv

# Filter abundance
time talon_abundance \
    --db  ../full_gencode_v29_SIRV_2020-02-11.db \
    --w whitelist_v4.4.2.csv \
    -a gencode_v29_SIRV  \
    -b hg38_SIRV \
    --o GM12878
```

# Plot number of reads by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f GM12878_talon_abundance_filtered.tsv \
    --ymax 10000000 \
    --datasets GM12878_R1,GM12878_R2 \
    -o .

# Plot number of reads by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f GM12878_talon_abundance_filtered.tsv \
    --ymax 10000000 \
    --datasets GM12878_R1 \
    -o .

# Plot number of reads by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f GM12878_talon_abundance_filtered.tsv \
    --ymax 10000000 \
    --datasets GM12878_R2 \
    -o .
```
