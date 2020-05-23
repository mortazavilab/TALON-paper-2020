# Running FLAIR on PacBio SIRV data

FLAIR was cloned from https://github.com/BrooksLabUCSC/flair on 8/5/2019.

1. Pull only the SIRV reads from the PacBio FLAIR run and format it like talon.
```bash
# Remove SIRVs and ERCCs
head -1 ../PacBio/counts_matrix.tsv > sirv_counts_matrix.tsv
grep -v "gSpikein_ERCC" ../PacBio/counts_matrix.tsv | \
        grep "SIRV" >> sirv_counts_matrix.tsv

python format_flair_matrix_like_talon.py sirv_counts_matrix.tsv sirv_counts_matrix_talon.tsv
```

2. Get the number of known transcripts detected in TALON, FLAIR, and both
```bash
python compare_transcripts.py \
   --f sirv_counts_matrix_talon.tsv \
   --t /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_sirv_talon_abundance_filtered.tsv
```
