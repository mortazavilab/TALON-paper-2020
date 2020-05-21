

1. Create a file listing the TSS and TES position for each annotated transcript:
```
python get_annotated_TSSs_and_TESs.py \
    --gtf ../refs/GENCODE_v29/gencode.v29.SIRV.ERCC.annotation.gtf \
    --o gencode.v29.SIRV.ERCC.TSS_and_TES.tsv
```

2. Plot distance of actual read start/emd from the annotated start/end associated with the annotated intron chain
```
python plot_TSS_and_TES_annot_dists.py \
    --f /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/GM12878_talon_read_annot.tsv \
    --ref gencode.v29.SIRV.ERCC.TSS_and_TES.tsv \
    --datasets PB_GM12878_R1,PB_GM12878_R2 \
    --xmax 800 --ymax 750000 --o GM12878_PacBio
```
