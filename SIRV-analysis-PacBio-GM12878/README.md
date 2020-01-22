# SIRV-analysis-PacBio-GM12878

The goal of this analysis is to sanity-check the pipeline on these known isoforms. This is also a test of how Minimap behaves on the SIRVs to determine whether we need to map them with different parameter settings (SIRVs have noncanonical SJs that may be mishandled by the new Minimap2 parameters). I'm going to look at the reads pre-TranscriptClean.

1. Get the reads that mapped to a SIRV
```
cat <(samtools view -H /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/Minimap2/mapped_FLNC.sam) <(awk '{if ($3 ~ "SIRV") print $0}' /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/Minimap2/mapped_FLNC.sam) > PacBio_Sequel2_GM12878_R1_SIRV_reads.sam

cat <(samtools view -H /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/Minimap2/mapped_FLNC.sam) <(awk '{if ($3 ~ "SIRV") print $0}' /share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/Minimap2/mapped_FLNC.sam) > PacBio_Sequel2_GM12878_R2_SIRV_reads.sam
```

2. Is there much multimapping going on, or was that minimal?
```
samtools view PacBio_Sequel2_GM12878_R1_SIRV_reads.sam | awk '{if($2 != '0' && $2 != '16') print $0}' | wc -l

samtools view PacBio_Sequel2_GM12878_R2_SIRV_reads.sam | awk '{if($2 != '0' && $2 != '16') print $0}' | wc -l
```
For Rep1, only 4 out of 6,924 alignments in the file are non-primary. In Rep2, 9 out of 13,821 alignments were non-primary. So overall, this does not seem like a big deal at all.

3. Profile the sequence errors in the SIRVs using Dry Run mode in TranscriptClean
```
qsub ./run_TC_dry_mode.sh
```

4. Do a SIRV-only TALON run directly on the mapped reads
```


``` 