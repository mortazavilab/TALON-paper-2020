# Compare different technologies' performance on simulated data

## Running FLAIR on the simulated data

1. Run align and correct steps separately on replicates.
```bash
# normal data
qsub run_flair_align_r1.sh
qsub run_flair_align_r2.sh

qsub run_flair_correct_r1.sh
qsub run_flair_correct_r2.sh

# perfect data
qsub run_flair_align_r1_perf.sh
qsub run_flair_align_r2_perf.sh

qsub run_flair_correct_r1_perf.sh
qsub run_flair_correct_r2_perf.sh
```	

2. Concatenate files from each replicate and run collapse
```bash
# normal data
r1_fasta=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1/rep1.fasta
r2_fasta=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2/rep2.fasta
cat r1_all_corrected.psl r2_all_corrected.psl > normal_corrected.psl
cat $r1_fasta $r2_fasta > normal.fasta
qsub run_flair_collapse_normal.sh

# perfect data
r1_p_fasta=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_perf/rep1_perf.fasta
r2_p_fasta=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_perf/rep2_perf.fasta
cat r1_perf_all_corrected.psl r2_perf_all_corrected.psl > perf_corrected.psl
cat $r1_p_fasta $r2_p_fasta > perf.fasta
qsub run_flair_collapse_perf.sh
```

3. Quantify
Config files for this step look like this:

Normal
```
norm_r1	norm	batch1	/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1/rep1.fasta
norm_r2	norm	batch1 /dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2/rep2.fasta
```

Perfect 
```
perf_r1	perf	batch1	/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_perf/rep1_perf.fasta
perf_r2	perf	batch1 /dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_perf/rep2_perf.fasta
```

Then run the quantify scripts
```bash
qsub run_flair_quantify_norm.sh
qsub run_flair_quantify_perf.sh
```

4. Remove spike ins (these shouldn't be in here but we'll do so to be sure) and format abundance matrices to look like TALON abundance matrices
```bash
flair_dir=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/compare_to_FLAIR/

# normal
grep -v "gSpikein_ERCC" norm_counts.tsv | \
        grep -v "SIRV" > norm_counts_nospikes.tsv
python ${flair_dir}format_flair_matrix_like_talon.py norm_counts_nospikes.tsv norm_talon_abundance.tsv

# perfect
grep -v "gSpikein_ERCC" perf_counts.tsv | \
        grep -v "SIRV" > perf_counts_nospikes.tsv
python ${flair_dir}format_flair_matrix_like_talon.py perf_counts_nospikes.tsv perf_talon_abundance.tsv
```

5. Performance of FLAIR on isoform detection
```bash
adir=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/compare_technologies/

r1=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1/rep1_headers.fasta
r2=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2/rep2_headers.fasta
python ${adir}compare_isoform_detection.py \
	-sim_files ${r1},${r2} \
	-sim_name normal \
	-ab norm_talon_abundance.tsv \
	-tech_name FLAIR \
	-datasets norm_r1_norm_batch1,norm_r2_norm_batch1

r1_perf=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_perf/rep1_perf_headers.fasta
r2_perf=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_perf/rep2_perf_headers.fasta
python ${adir}compare_isoform_detection.py \
	-sim_files ${r1},${r2} \
	-sim_name perfect \
	-ab perf_talon_abundance.tsv \
	-tech_name FLAIR \
	-datasets perf_r1_perf_batch1,perf_r2_perf_batch1
```

6. Performance of FLAIR on transcript and gene level quantification
```bash
adir=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/compare_technologies/
gtf=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/gencode.v29.SIRV.ERCC.annotation.gtf

r1=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1/rep1_headers.fasta
r2=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2/rep2_headers.fasta
python ${adir}compare_quantification.py \
	-sim_files ${r1},${r2} \
	-sim_name normal \
	-t_ab norm_talon_abundance.tsv \
	-g_ab norm_talon_abundance.tsv \
	-tech_name FLAIR \
	-ref_gtf $gtf \
	-datasets 'norm_r1_norm_batch1,norm_r2_norm_batch1'

r1_perf=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_perf/rep1_perf_headers.fasta
r2_perf=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_perf/rep2_perf_headers.fasta
python ${adir}compare_quantification.py \
	-sim_files ${r1_perf},${r2_perf} \
	-sim_name perfect \
	-t_ab perf_talon_abundance.tsv \
	-g_ab perf_talon_abundance.tsv \
	-tech_name FLAIR \
	-ref_gtf $gtf \
	-datasets 'perf_r1_perf_batch1,perf_r2_perf_batch1'
```

```
# Normal
Simulated dataset: normal, norm_r1_norm_batch1, Technology: FLAIR
Gene correlation: 0.9536902976522196, gene pval: 0.0
Transcript correlation: 0.5936932792166355, transcript pval: 0.0
Simulated dataset: normal, norm_r2_norm_batch1, Technology: FLAIR
Gene correlation: 0.984062616549627, gene pval: 0.0
Transcript correlation: 0.684465595334091, transcript pval: 0.0

# Perfect
Gene correlation: 0.9060328826821611, gene pval: 0.0
Transcript correlation: 0.48897560416057523, transcript pval: 0.0
Simulated dataset: perfect, perf_r2_perf_batch1, Technology: FLAIR
Gene correlation: 0.9080603176920531, gene pval: 0.0
Transcript correlation: 0.49528866804526006, transcript pval: 0.0
```