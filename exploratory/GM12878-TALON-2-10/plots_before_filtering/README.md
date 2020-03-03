## Fraction As at end of transcript analysis
```
python ../../../plotting_scripts/plot_frac_As_by_read_annot_category.py \
    --f ../GM12878_talon_read_annot.tsv \
    --outprefix GM12878

python ../../../plotting_scripts/plot_percent_IP_by_read_annot_category.py \
    --f ../GM12878_talon_read_annot.tsv \
    --outprefix GM12878
```

## Isoforms detected
```
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f ../GM12878_talon_abundance.tsv \
        --datasets GM12878_R1,GM12878_R2 \
        --ymax 400000 \
        -o .
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f ../GM12878_talon_abundance.tsv \
        --datasets GM12878_R1 \
        --ymax 400000 \
        -o .
Rscript ../../../plotting_scripts/plot_novelty_categories_distinct_isoforms.R \
        --f ../GM12878_talon_abundance.tsv \
        --datasets GM12878_R2 \
        --ymax 400000 \
        -o .
```

# Plot number of reads by category (both reps together)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../GM12878_talon_abundance.tsv \
    --ymax 10000000 \
    --datasets GM12878_R1,GM12878_R2 \
    -o .

# Plot number of reads by category (Rep1)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../GM12878_talon_abundance.tsv \
    --ymax 10000000 \
    --datasets GM12878_R1 \
    -o .

# Plot number of reads by category (Rep2)
Rscript ../../../plotting_scripts/plot_novelty_category_read_counts.R \
    --f ../GM12878_talon_abundance.tsv \
    --ymax 10000000 \
    --datasets GM12878_R2 \
    -o .
```
