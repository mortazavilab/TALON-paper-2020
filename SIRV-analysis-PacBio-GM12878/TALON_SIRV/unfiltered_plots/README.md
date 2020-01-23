# Visualization of TALON results

## Distinct isoforms by category
Plot number of UNFILTERED distinct isoforms by category (both reps together)
```
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f ../SIRV_talon_abundance.tsv \
        --datasets SIRV_Rep1,SIRV_Rep2 \
        -o .
```
Plot number of UNFILTERED distinct isoforms by category (Rep1)
```
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f ../SIRV_talon_abundance.tsv \
        --datasets SIRV_Rep1 \
        -o .
```
Plot number of UNFILTERED distinct isoforms by category (Rep2)
```
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f ../SIRV_talon_abundance.tsv \
        --datasets SIRV_Rep2 \
        -o .
```

## Reads by category
Plot number of UNFILTERED reads by category (both reps together)
```
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../SIRV_talon_abundance.tsv \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    -o .
```
Plot number of UNFILTERED reads by category (Rep1)
```
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts_one_dataset.R \
    --f ../SIRV_talon_abundance.tsv \
    --dataset SIRV_Rep1 \
    -o .
```
Plot number of UNFILTERED reads by category (Rep2)
```
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts_one_dataset.R \
    --f ../SIRV_talon_abundance.tsv \
    --dataset SIRV_Rep2 \
    -o .
```

