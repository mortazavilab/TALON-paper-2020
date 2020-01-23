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

3. Profile the sequence errors in the SIRVs using TranscriptClean. This command includes report generation.
```
qsub ./run_TC_SIRVs.sh
```
Observing the reports, I notice that the NCSJ rate is around 5% in both replicates, and TC correction brings that down to about 1.5%. I think that the 5% number is similar to what we see in Sequel datasets, although I would need to check. The 1.5% number is much lower than we were typically able to get in the past, which would potentially suggest that the SequelII has a slightly different error profile (with more NCSJs coming from small indel errors and fewer coming from RT priming errors). That being said, I noticed that there is a read with a 1270 bp insertion, which obviously is a terrible mapping. 

I loaded the SIRV reads (before/after TC) into the IGV browser and inspected them. The corrections are appearing as expected. It will be interesting to quantify how much fake novelty we see with TALON. Some of the SIRVs seem very well-behaved in Rep2, for instance, but SIRV3 has a fair amount of junk.

4. Do a SIRV-only TALON run on the corrected reads
```


``` 
