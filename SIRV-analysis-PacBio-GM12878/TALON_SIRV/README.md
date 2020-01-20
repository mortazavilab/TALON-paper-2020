# Visualization of TALON results

Plot number of UNFILTERED distinct isoforms by category (both reps together)
```
Rscript ../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f SIRV_talon_abundance.tsv \
        --datasets SIRV_Rep1,SIRV_Rep2 \
        -o plots
```
Plot number of UNFILTERED distinct isoforms by category (Rep1)
```
Rscript ../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f SIRV_talon_abundance.tsv \
        --datasets SIRV_Rep1 \
        -o plots
```
Plot number of UNFILTERED distinct isoforms by category (Rep1)
```
Rscript ../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f SIRV_talon_abundance.tsv \
        --datasets SIRV_Rep2 \
        -o plots
```
