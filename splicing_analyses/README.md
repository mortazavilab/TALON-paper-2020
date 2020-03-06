# Splicing Analysis

We did a lot of extra analysis of our data in terms of looking at splice junctions. The first step for running any of the sub-analyses in this directory is dependent on first obtaining the splice junctions from each dataset, which will be detailed here. 

1. Get the tables from the supplemental tables file that we'll be using, and set other paths that we'll be using (to the hg38 reference genome).
```bash
mkdir figures

# download the supplementary tables and change this path!
DATA=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/

gm_pb_gtf=${DATA}pb_talon.gtf
gm_ont_gtf=${DATA}ont_talon.gtf

REFPATH=~/mortazavi_lab/ref/hg38/
```

2. Extract splice junctions from GM12878 PacBio and ONT gtfs using TranscriptClean
```bash
python get_SJs_from_gtf.py \
	--f ${gm_pb_gtf} \
	--g ${REFPATH}hg38.fa \
	--o pb_talon_GM12878_sjs.tab

python get_SJs_from_gtf.py \
	--f ${gm_ont_gtf} \
	--g ${REFPATH}hg38.fa \
	--o ont_talon_GM12878_sjs.tab
```

3. Now, let's get the splice junctions present in the Illumina data by mapping with STAR. 
```bash
qsub run_STAR_illumina_GM12878_Rep1.sh
qsub run_STAR_illumina_GM12878_Rep2.sh
```

4. Filter out novel Illumina SJs that don't have support from both reps and known Illumina SJs that have no read support.
```bash
python filter_illumina_sjs.py \
	-sj_1 GM12878_Rep1_alignedSJ.out.tab \
	-sj_2 GM12878_Rep2_alignedSJ.out.tab 
```