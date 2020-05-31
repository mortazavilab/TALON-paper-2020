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

## Comparing counts per SIRV isoform in TALON and FLAIR
1. Merge together the two abundance files for known SIRV transcripts only
```python
import pandas as pd
flair_abd = pd.read_csv("sirv_counts_matrix_talon.tsv", sep = '\t', header = 0)
talon_abd = pd.read_csv("/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/pb_sirv_talon_abundance_filtered.tsv", sep = '\t', header = 0)

# Filter to known transcripts only and remove extra columns
flair_abd = flair_abd.loc[flair_abd.transcript_novelty == "Known"][["transcript_ID","GM12878_Rep1_GM12878_batch1", "GM12878_Rep2_GM12878_batch1"]]
talon_abd = talon_abd.loc[talon_abd.transcript_novelty == "Known"][["annot_transcript_id","PB_GM12878_R1", "PB_GM12878_R2"]]

# Rename cols
flair_abd = flair_abd.rename(columns={"GM12878_Rep1_GM12878_batch1": "FLAIR R1", 
                                      "GM12878_Rep2_GM12878_batch1": "FLAIR R2",
                                      "transcript_ID": "annot_transcript_id"})
talon_abd = talon_abd.rename(columns={"PB_GM12878_R1": "TALON R1",
                                      "PB_GM12878_R2": "TALON R2"})

# Merge
merged = pd.merge(talon_abd, flair_abd, on = "annot_transcript_id", how = "outer").fillna(0)
merged.to_csv("merged_talon_flair_pacbio_GM12878_isoform_abd.tsv", sep = "\t", index = None, header = None)
```

