## Illumina/Gencode splice junction support for whole transcripts

We want to determine what the percentage of PacBio TALON transcripts with all of their splice junctions supported by short-read data is. We decided to subset on novelty to provide some more insight. We'll show this step-by-step example for GM12878 first.

1. Use the GM12878 PacBio GTF and Illumina SJs to see how many transcripts of each category have all of their splice junctions supported by Illumina data.
```bash
DATA=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/
gm_gtf=${DATA}pb_talon.gtf 

python ../../get_isoform_sj_support.py \
	-gtf ${gm_gtf} \
	-ref_sj_1 GM12878_alignedSJ.out.tab \
	-ref_sj_2 gencode_v29_sjs.tab \
	-sample pb_GM12878
```

2. Generate a summary table counting how many transcripts in each novelty category are supported by GENCODE or Illumina
```bash
python ../../gen_isoform_support_table.py \
	-csv pb_GM12878_isoform_sj_support.csv \
	-sample pb_GM12878
```

3. Visualize Illumina support vs. isoform novelty type. Keep in mind that a transcript towards having Illumina support must have ALL of its splice junctions supported by the Illumina data.
```bash
python ../../plot_isoform_sj_support_by_novelty.py \
	-c pb_GM12878_isoform_sj_support_summary.csv \
	-s "PacBio GM12878"
```

<img align="center" width="500" src="figures/PacBio_GM12878_sj_novelty_illumina_isoform_novelty_support.png">

4. Also run this on the ONT data
```bash
# ONT GM12878
ont_gm_gtf=${DATA}ont_talon.gtf
python ../../get_isoform_sj_support.py \
	-gtf ${ont_gm_gtf} \
	-ref_sj_1 GM12878_alignedSJ.out.tab \
	-ref_sj_2 gencode_v29_sjs.tab \
	-sample ont_GM12878

python ../../gen_isoform_support_table.py \
	-csv ont_GM12878_isoform_sj_support.csv \
	-sample ont_GM12878

python ../../plot_isoform_sj_support_by_novelty.py \
	-c ont_GM12878_isoform_sj_support_summary.csv \
	-s "ONT GM12878"
```

<img align="center" width="500" src="figures/ONT_GM12878_sj_novelty_illumina_isoform_novelty_support.png">