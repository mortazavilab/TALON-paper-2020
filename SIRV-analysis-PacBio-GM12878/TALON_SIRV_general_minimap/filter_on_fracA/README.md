```
python filter_on_fracA.py --f ../SIRV_talon_read_annot.tsv --a ../unfiltered_plots/SIRV_Rep1_fraction_As_10bp_after_transcript.tsv --o SIRV

# Filter abundance
source activate TALON-v4.4.2
filtered=SIRV_talon_abundance_filtered.tsv
time talon_abundance \
    --db ../SIRV_2020-01-23.db \
    --w SIRV_frac_A_whitelist.csv \
    -a SIRV \
    -b SIRV \
    --o SIRV
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
# Plot number of  reads by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f $filtered \
    --datasets SIRV_Rep1,SIRV_Rep2 \
    -o .

# Plot number of  reads by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts_one_dataset.R \
    --f $filtered \
    --dataset SIRV_Rep1 \
    -o .

# Plot number of  reads by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts_one_dataset.R \
    --f $filtered \
    --dataset SIRV_Rep2 \
    -o .
```
