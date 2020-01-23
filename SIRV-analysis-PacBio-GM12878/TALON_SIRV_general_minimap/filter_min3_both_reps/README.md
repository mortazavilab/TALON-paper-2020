# Visualization of TALON results: new filtering strategy

## Filter abundance file
Require transcript to be 1) Known, or 2) Appear in each replicate at least 3x.  
```
filtered=SIRV_talon_abundance_custom_filter.tsv
head -1 ../SIRV_talon_abundance.tsv > $filtered
awk '{ if ($10 == 'Known') print $0 }' ../SIRV_talon_abundance.tsv >> $filtered
tail -n +2 ../SIRV_talon_abundance.tsv | \
      awk '{ if ($10 != 'Known' && $12 > 2 && $13 > 2) print $0 }' >> $filtered
```

## Distinct isoforms by category
```
# Plot number of distinct isoforms by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f $filtered \
        --datasets SIRV_Rep1,SIRV_Rep2 \
        -o .

# Plot number of distinct isoforms by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f $filtered \
        --datasets SIRV_Rep1 \
        -o .

# Plot number of distinct isoforms by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f $filtered \
        --datasets SIRV_Rep2 \
        -o .
```

## Reads by category
```
# Plot number of UNFILTERED reads by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f $filtered \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    -o .

# Plot number of UNFILTERED reads by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts_one_dataset.R \
    --f $filtered \
    --dataset SIRV_Rep1 \
    -o .

# Plot number of UNFILTERED reads by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts_one_dataset.R \
    --f $filtered \
    --dataset SIRV_Rep2 \
    -o .
```

