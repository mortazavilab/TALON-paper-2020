# Testing TALON filter on SIRV reads

## Panel A: Reads by category (unfiltered)
```
PLOTPATH=../../plotting_scripts
abd=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_sirv_talon_abundance.tsv
mkdir -p plots

Rscript $PLOTPATH/plot_novelty_category_read_counts.R \
    --f $abd \
    --ymax 15000 \
    --datasets PB_GM12878_R1,PB_GM12878_R2 \
    -o plots/
```

## Panel B: Distinct transcript models, unfiltered
```
PLOTPATH=../../plotting_scripts
abd=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_sirv_talon_abundance.tsv
mkdir -p plots

Rscript $PLOTPATH/plot_novelty_categories_distinct_isoforms.R \
        --f $abd \
        --datasets PB_GM12878_R1,PB_GM12878_R2 \
        --ymax 350 \
        -o plots/
```

## Panel C: Fraction internal priming by category (unfiltered)
```
PLOTPATH=../../plotting_scripts
read_annot=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/GM12878_talon_read_annot.tsv
mkdir -p plots

# Get read_annot file of only SIRVs from GM12878
head -1 ${read_annot} > SIRV_talon_read_annot.tsv
awk '{if($4 ~ "SIRV" && $2 ~ "PB") print $0}' ${read_annot} >> SIRV_talon_read_annot.tsv

source activate python3.6
python $PLOTPATH/plot_percent_IP_by_read_annot_category.py \
    --f SIRV_talon_read_annot.tsv \
    --omitGenomic \
    --outprefix plots/SIRV
```

## Panel D: Reads by category (filtered)
```
PLOTPATH=../../plotting_scripts
filt_abd=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_sirv_talon_abundance_filtered.tsv
mkdir -p post-filter

Rscript $PLOTPATH/plot_novelty_category_read_counts.R \
    --f $filt_abd \
    --ymax 15000 \
    --datasets PB_GM12878_R1,PB_GM12878_R2 \
    -o post-filter/
```

## Panel E: Distinct transcript models, filtered
```
PLOTPATH=../../plotting_scripts
filt_abd=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_sirv_talon_abundance_filtered.tsv
mkdir -p post-filter

Rscript $PLOTPATH/plot_novelty_categories_distinct_isoforms.R \
        --f $filt_abd \
        --datasets PB_GM12878_R1,PB_GM12878_R2 \
        --ymax 350 \
        -o post-filter/
```
