# Supplement

Files/paths:
```
PLOTPATH=../plotting_scripts
DATA=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON
```

## Figure S? and Figure S?: TALON read length distributions for PacBio and ONT GM12878 datasets

Remove SIRV and ERCC reads
```
awk '{if($4 !~ "ERCC" && $4 !~ "SIRV") print $0}' \
    $DATA/GM12878_talon_read_annot.tsv > read_lengths/GM12878_main_read_annot.tsv
```
Plot
```
python ${PLOTPATH}/plot_read_length_distributions.py \
    --r read_lengths/GM12878_main_read_annot.tsv \
     --datasets PB_GM12878_R1,PB_GM12878_R2,ONT_GM12878_R1,ONT_GM12878_R2 \
     --map read_lengths/read_length_name_mapping.csv \
     --o read_lengths/
```
See resulting plots [here](https://github.com/dewyman/TALON-paper-2020/tree/master/Supplement/read_lengths).
