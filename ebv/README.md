# EBV

We wanted to see if TALON and long-read sequencing could be used to detect, characterize, and quantify transcripts from the EBV chromosome used to immortalize the GM12878 cell line. 


## 1. Initialize some vars
```bash
rep1_sam=/data/users/freese/TALON_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R1/label_reads/PacBio_Rep1_labeled.sam
rep2_sam=/data/users/freese/TALON_data/revisions_1-20/data/PacBio_Sequel2_GM12878_R2/label_reads/PacBio_Rep2_labeled.sam
```

## 2. Grab all reads that mapped to chrEBV, along with the suyam headers
```bash
grep ^@ $rep1_sam > ebv_rep1.sam
grep chrEBV $rep1_sam >> ebv_rep1.sam
grep ^@ $rep2_sam > ebv_rep2.sam
grep chrEBV $rep2_sam >> ebv_rep2.sam

sed -i 's/chrEBV/chr1/' ebv_rep1.sam
sed -i 's/chrEBV/chr1/' ebv_rep2.sam
```

## 3. Initialize TALON db with custom EBV gtf (using chr1 because sam files are fake mapped to chr1)
```bash
sed 's/chrEBV/chr1/' ebv.gtf > ebv_chr1.gtf
talon_initialize_database \
    --f ebv_chr1.gtf \
    --g HHV4 \
    --a HHV4 \
    --o ebv
```

### Run TALON on EBV reads extracted from mapped GM12878 sam files
```bash
talon \
    --f ebv_talon_config.csv \
    --db ebv.db \
    --build HHV4 \
    --o ebv
```

## Post-TALON files

### Get a whitelist of transcripts via TALON filtering
```bash
talon_filter_transcripts \
    --db ebv.db \
    -a HHV4 \
    --maxFracA 0.5 \
    --o ebv_whitelist 
```

### Get filtered GTF file
```bash
talon_create_GTF \
      --db ebv.db \
      --b HHV4 \
      --a HHV4 \
      --o ebv \
      --whitelist ebv_whitelist
```

### Get unfiltered abundance file
```bash
talon_abundance \
    --db ebv.db \
    --a HHV4 \
    --b HHV4 \
    --o ebv
```

### Get filtered abundance file
```bash
talon_abundance \
    --db ebv.db \
    --a HHV4 \
    --b HHV4 \
    --whitelist ebv_whitelist \
    --o ebv
```

### GM12878 files to compare with 
```bash
gm_gtf=/data/users/freese/TALON_data/revisions_1-20/human_TALON/pb_talon.gtf
gm_ab=/data/users/freese/TALON_data/revisions_1-20/human_TALON/pb_talon_abundance.tsv
gm_filt_ab=/data/users/freese/TALON_data/revisions_1-20/human_TALON/pb_talon_abundance_filtered.tsv
```

## Plotting
### Generate EBV abundance violin plots
```bash
python ebv_compute_tpms.py \
    --human_gtf $gm_gtf \
    --human_filt_ab $gm_filt_ab \
    --human_ab $gm_ab \
    --ebv_filt_ab ebv_talon_abundance_filtered.tsv \
    --ebv_ab ebv_talon_abundance.tsv \
    --datasets PB_GM12878_R1,PB_GM12878_R2 \
    --o pb
```
```bash
Rscript plot_ebv_v_human_abundances.R \
          --gene_csv pb_ebv_human_gene_abundance.csv \
          --transcript_csv pb_ebv_human_transcript_abundance.csv \
          --datasets combined \
          --o pb
```

<img align="center" width="300" src="combined_genes_ebv_human.png "><img align="center" width="300" src="combined_transcripts_ebv_human.png ">

## Genome browser trackline
### Generate tracklines using above GTF
```bash
python ../analysis_scripts/gen_novelty_tracks_gtf.py \
          --c ebv_gtf_track_config.csv
url=`cut -d, -f5 ebv_gtf_track_config.csv`
n=`cut -d, -f2 ebv_gtf_track_config.csv`
cp ebv_chr1.gtf ebv_talon_tracks/
printf 'track name="EBV Reference" visibility=pack color=0,0,128\n%s/ebv_chr1.gtf' "$url" >> ebv_talon_tracks/ebv_talon_${n}_tracks
```

From here, you can open the genome browser up and display your tracklines, and use the genome browser's PDF functionality or take a screenshot to get the genome browser shot. 

<img align="center" width="700" src="ebv_tracks.png ">



<!-- 
1. Download the ENCODE CAGE GM12878 data in bed format, and extract only EBV TSSs, and cat them together
```bash
wget https://www.encodeproject.org/files/ENCFF383NVJ/@@download/ENCFF383NVJ.bed.gz
wget https://www.encodeproject.org/files/ENCFF016XXM/@@download/ENCFF016XXM.bed.gz
gunzip ENCFF383NVJ.bed.gz
gunzip ENCFF016XXM.bed.gz

grep chrEBV ENCFF383NVJ.bed > ENCFF383NVJ_ebv_only.bed
grep chrEBV ENCFF016XXM.bed > ENCFF016XXM_ebv_only.bed

cat ENCFF383NVJ_ebv_only.bed > cage_ebv_TSSs.bed
cat ENCFF016XXM_ebv_only.bed >> cage_ebv_TSSs.bed
```

2. Replace chr1 with chrEBV in output TALON gtf so that it matches up with the EBV CAGE peaks.
```bash
sed 's/chr1/chrEBV/g' ebv_talon_observedOnly.gtf > ebv_talon_observedOnly_ebv.gtf
```

3. Intersect EBV CAGE peaks with TSSs found in our GM12878 data and plot!
```bash
mkdir figures

python run_CAGE_analysis.py \
  --gtf ebv_talon_observedOnly_ebv.gtf \
  --cage cage_ebv_TSSs.bed \
  --maxdist 100 \
  --o ./ebv

Rscript plot_support_by_novelty_type.R \
  --f ebv_CAGE_results.csv \
  --t CAGE \
  --novelty transcript_beds/ebv_novelty.csv \
  --splitISM \
  --ymax 15 \
  -o figures/ebv
```

<img align="center" width="700" src="figures/ebv_CAGE_support.png">

## Computational PAS analysis

See if TALON EBV transcript 3' ends are supported by poly-A recognition motifs.

```bash
python run_computational_PAS_analysis.py \
  --gtf ebv_talon_observedOnly_ebv.gtf \
  --genome ebv.fasta \
  --maxdist 35 \
  --o ./ebv

Rscript plot_support_by_novelty_type.R \
  --f ebv_polyA_motif.csv \
  --t PAS-comp \
  --novelty transcript_beds/ebv_novelty.csv \
  --ymax 15 \
  --splitISM \
  -o figures/ebv
```

<img align="center" width="700" src="figures/ebv_PAS-comp_support.png">



 -->