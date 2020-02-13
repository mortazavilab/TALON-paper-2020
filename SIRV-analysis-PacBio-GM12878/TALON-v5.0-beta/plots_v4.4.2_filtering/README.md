
Run filtering settings comparable to v4.4.2
```
source activate python3.6

# Make whitelist
talon_filter_transcripts \
    --db ../SIRV_2020-02-10.db \
    -a SIRV \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    --maxFracA 1 \
    --minCount 1 \
    --minDatasets 2 \
    --o SIRV_whitelist_v4.4.2.csv

# Filter abundance
time talon_abundance \
    --db ../SIRV_2020-02-10.db \
    --w SIRV_whitelist_v4.4.2.csv \
    -a SIRV \
    -b SIRV \
    --o SIRV
```
## Isoforms detected
```
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f SIRV_talon_abundance_filtered.tsv \
        --datasets SIRV_Rep1,SIRV_Rep2 \
        --ymax 350 \
        -o .
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f SIRV_talon_abundance_filtered.tsv \
        --datasets SIRV_Rep1 \
        --ymax 350 \
        -o .
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f SIRV_talon_abundance_filtered.tsv \
        --datasets SIRV_Rep2 \
        --ymax 350 \
        -o .
```

## Reads by category
```
# Plot number of reads by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f SIRV_talon_abundance_filtered.tsv \
    --ymax 10000 \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    -o .

# Plot number of reads by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f SIRV_talon_abundance_filtered.tsv \
    --ymax 10000 \
    --datasets SIRV_Rep1 \
    -o .

# Plot number of reads by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f SIRV_talon_abundance_filtered.tsv \
    --ymax 10000 \
    --datasets SIRV_Rep2 \
    -o .
```
