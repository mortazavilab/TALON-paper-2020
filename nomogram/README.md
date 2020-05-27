We want to investigate what PacBio sequencing depth is necessary to reach plateuing expression levels. We subsample the reads and evaluate what proportion of genes or transcript expression levels are quantified within 5% of the expression levels we recorded with all of our reads.

1. Grab the read annots for PacBio and ONT. Remove SIRVs and ERCCs
```bash
annot=/data/users/freese/TALON_data/revisions_1-20/human_TALON/GM12878_talon_read_annot.tsv
head -1 $annot > pb_annots.tsv
head -1 $annot > ont_annots.tsv
grep -v ERCC $annot | grep -v -E "SIRV[0-7]" | grep 'PB_GM12878_R1\|PB_GM12878_R2' >> pb_annots.tsv
grep -v ERCC $annot | grep -v -E "SIRV[0-7]" | grep 'ONT_GM12878_R1\|ONT_GM12878_R2' >> ont_annots.tsv
cat pb_annots.tsv | cut -f2,10,11,16,17 > pb_mini_annots.tsv
cat ont_annots.tsv | cut -f2,10,11,16,17 > ont_mini_annots.tsv
```

2. Subsample reads at various depths
```bash 
pb_reads=pb_mini_annots.tsv
ont_reads=ont_mini_annots.tsv

python subsample_and_plot.py \
	--f $pb_reads \
	--prefix pb

python subsample_and_plot.py \
	--f $ont_reads \
	--intervals 500000,750000,1000000,1250000 \
	--whitelist $ont_whitelist \
	--prefix ont

```