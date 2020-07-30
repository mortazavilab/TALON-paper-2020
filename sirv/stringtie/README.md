# Running StringTie2 on our data

According to this workflow http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de

StringTie was downloaded from git on 7/29/2020

1. Run StringTie2 on each mapped sample
```bash 
ref_annot=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/
r1=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/label_reads/PacBio_Rep1_labeled.bam
r2=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/label_reads/PacBio_Rep2_labeled.bam

# rep 1
prefix='r1'
${sdir}./stringtie \
	-L \
	-p 16 \
	-c 1 \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv
	$r1

# rep 2
prefix='r2'
${sdir}./stringtie \
	-L \
	-p 16 \
	-c 1 \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv
	$r2
```

Neither r1_stringtie.gtf nor r2_stringtie.gtf have any transcripts with coverage == 0. Only one transcript has coverage < 1, which is 0.994865  STRG.11292.3  SIRV4.1-SIRV404 from r1_stringtie.gtf.

I also ran my mini script to see how many known and novel SIRVs were detected in each replicate
```bash
prefix='r1'
python get_stringtie_sirvs_onefile.py \
	-gtf ${prefix}_stringtie.gtf \
	-ref_id_field reference_id

prefix='r2'
python get_stringtie_sirvs_onefile.py \
	-gtf ${prefix}_stringtie.gtf \
	-ref_id_field reference_id
```

```
# Replicate 1 
StringTie2 found 60 known SIRV isoforms
StringTie2 found 3 novel SIRV isoforms

# Replicate 2 
StringTie2 found 63 known SIRV isoforms
StringTie2 found 2 novel SIRV isoforms
```

2. Merge the StringTie output gtfs to get the union of the two annotations
```bash
printf "r1_stringtie.gtf\nr2_stringtie.gtf" > stringtie_gtfs.tsv

gtf_list=stringtie_gtfs.tsv
ref_annot=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf
sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/

${sdir}./stringtie \
	--merge \
	-p 16 \
	-c 1 \
	-G ${ref_annot} \
	-o merged_stringtie.gtf \
	-l st_novel_ \
	$gtf_list
```

I wanted to check what the number of novel and known isoforms was from the resulting merged GTF. The reference IDs are now inexplicably in the "transcript_id" field as opposed to the "reference_id" field.

```bash
python get_stringtie_sirvs_onefile.py \
	-gtf merged_stringtie.gtf \
	-ref_id_field transcript_id
```

```
StringTie2 found 69 known SIRV isoforms
StringTie2 found 5 novel SIRV isoforms
```

Good job, StringTie!

3. Rerun each of the replicates using the merged GTF as the reference annotation. This will ensure concordance of novel transcript ids across the datasets as well as allow for us to quantify on each replicate. In this case, we will use the `-e` option, as StringTie has already found all the novel transcripts we needed to find, and therefore we only want to include models that have already been annotated by StringTie.
```bash
sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/
r1=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/label_reads/PacBio_Rep1_labeled.bam
r2=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/label_reads/PacBio_Rep2_labeled.bam
ref_annot=merged_stringtie.gtf

# rep 1
prefix='r1_merged'
${sdir}./stringtie \
	-L \
	-p 16 \
	-c 1 \
	-e \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv \
	$r1

# rep 2
prefix='r2_merged'
${sdir}./stringtie \
	-L \
	-p 16 \
	-c 1 \
	-e \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv \
	$r2
```

Finally, we'll examine the output GTFs for novelty. Which again, inexplicably, now use the "transcript_id" field to list the id from the reference.

```bash
prefix='r1_merged'
python get_stringtie_sirvs_onefile.py \
	-gtf ${prefix}_stringtie.gtf \
	-ref_id_field transcript_id

prefix='r2_merged'
python get_stringtie_sirvs_onefile.py \
	-gtf ${prefix}_stringtie.gtf \
	-ref_id_field transcript_id
```

And what do we have here... a bunch of transcripts that have 0 coverage, what even? 
```
# rep 1
Removed 15 0-coverage transcripts
StringTie2 found 56 known SIRV isoforms
StringTie2 found 3 novel SIRV isoforms

# rep 2
Removed 12 0-coverage transcripts
StringTie2 found 60 known SIRV isoforms
StringTie2 found 2 novel SIRV isoforms
```

Finally, I'll merge the gtfs on transcript id alone to see how many unique known and novel isoforms were found in the combined replicates

```bash
python get_stringtie_sirvs.py \
	-gtfs 'r1_merged_stringtie.gtf,r2_merged_stringtie.gtf' \
	-ref_id_field transcript_id
```

```
Removed 15 0-coverage transcripts
Removed 12 0-coverage transcripts

StringTie2 found 60 known SIRV isoforms
StringTie2 found 3 novel SIRV isoforms
```

Interestingly it did worse when being assembled with the merged GTF?


4. Just kidding, maybe this is the last step. Now, I'll try removing the `-e` option to see if maybe this is causing the 0-coverage isoforms to make it through the pipeline. 

```bash
sdir=/dfs3/pub/freese/mortazavi_lab/bin/stringtie-2.1.4.Linux_x86_64/
r1=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/label_reads/PacBio_Rep1_labeled.bam
r2=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/label_reads/PacBio_Rep2_labeled.bam
ref_annot=merged_stringtie.gtf

# rep 1
prefix='r1_merged_no_e'
${sdir}./stringtie \
	-L \
	-p 16 \
	-c 1 \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv \
	$r1

# rep 2
prefix='r2_merged_no_e'
${sdir}./stringtie \
	-L \
	-p 16 \
	-c 1 \
	-G ${ref_annot} \
	-o ${prefix}_stringtie.gtf \
	-A ${prefix}_abundance.tsv \
	$r2
```

Again, examine each replicate for number of known/novel isoforms detected. This time with "reference_id"

```bash
prefix='r1_merged_no_e'
python get_stringtie_sirvs_onefile.py \
	-gtf ${prefix}_stringtie.gtf \
	-ref_id_field reference_id

prefix='r2_merged_no_e'
python get_stringtie_sirvs_onefile.py \
	-gtf ${prefix}_stringtie.gtf \
	-ref_id_field reference_id
```

Nope, no 0-coverage transcripts.

```
# rep 1
StringTie2 found 49 known SIRV isoforms
StringTie2 found 13 novel SIRV isoforms

# rep 2
StringTie2 found 53 known SIRV isoforms
StringTie2 found 11 novel SIRV isoforms
```

Finally, get the detected merged isoforms. 

```bash
python get_stringtie_sirvs.py \
	-gtfs 'r1_merged_no_e_stringtie.gtf,r2_merged_no_e_stringtie.gtf' \
	-ref_id_field reference_id
```

```
StringTie2 found 54 known SIRV isoforms
StringTie2 found 103 novel SIRV isoforms
```