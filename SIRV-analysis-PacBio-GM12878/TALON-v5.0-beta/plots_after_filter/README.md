## Reads by category
```
# Plot number of UNFILTERED reads by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../SIRV_talon_abundance_filtered.tsv \
    --ymax 10000 \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    -o .

# Plot number of UNFILTERED reads by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../SIRV_talon_abundance_filtered.tsv \
    --ymax 10000 \
    --datasets SIRV_Rep1 \
    -o .

# Plot number of UNFILTERED reads by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../SIRV_talon_abundance_filtered.tsv \
    --ymax 10000 \
    --datasets SIRV_Rep2 \
    -o .
```
