# Figure 4: 5' and 3' completeness by novelty category

Files/paths used to generate the panels of this figure:
```bash
PLOTPATH=../plotting_scripts
OUTPLOTS=plots
mkdir -p ${OUTPLOTS}

data_dir=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON
PB_GTF=$data_dir/pb_talon.gtf
ONT_GTF=$data_dir/ont_talon.gtf
CAGE=../CAGE_data/FANTOM5/hg38_CAGE.bed
RNAPET=../RNA-PET_data/data/GM12878_hg38.bed
GENOME=../refs/hg38_SIRV/hg38_SIRV.fa
```
GTF files are available as supplementary tables of the TALON paper.  
To obtain FANTOM5 CAGE data, please see instructions [here](https://github.com/dewyman/TALON-paper-2020/blob/master/CAGE_data).  
To obtain ENCODE RNA-PET data, please see instructions [here](https://github.com/dewyman/TALON-paper-2020/tree/master/RNA-PET_data).  

Software versions:  
* Python 3.7.2
* Pyfasta 0.5.2 
* Bedtools v2.27.1
* R v3.5.1

## Panel A: Percentage of TALON transcript models with CAGE support for their 5' end by novelty category (GM12878 PacBio)
```bash
source activate mypython3.7.2
OUT=CAGE/PacBio_GM12878
mkdir -p ${OUT}
python run_CAGE_analysis.py \
        --gtf ${PB_GTF} \
        --cage ${CAGE} \
        --maxdist 100 \
        --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_CAGE_results.csv \
    --t CAGE \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --splitISM \
    --ymax 31000 \ 
    -o ${OUTPLOTS}/GM12878_PacBio
```
<img align="center" width="400" src="plots/GM12878_PacBio_CAGE_support.png">

## Panel B: Percentage of TALON transcript models with CAGE support for their 5' end by novelty category (GM12878 ONT)
```bash
OUT=CAGE/ONT_GM12878
mkdir -p ${OUT}
python run_CAGE_analysis.py \
        --gtf ${ONT_GTF} \
        --cage ${CAGE} \
        --maxdist 100 \
        --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_CAGE_results.csv \
    --t CAGE \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --splitISM \
    --ymax 31000 \
    -o ${OUTPLOTS}/GM12878_ONT
```
<img align="center" width="400" src="plots/GM12878_ONT_CAGE_support.png">

## Panel C: Percentage of TALON transcript models with a poly(A) motif identified at their 3' end (GM12878 PacBio)
```bash
OUT=PAS-comp/PacBio_GM12878
mkdir -p ${OUT}
python run_computational_PAS_analysis.py \
        --gtf ${PB_GTF} \
        --genome ${GENOME} \
        --maxdist 35 \
        --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_polyA_motif.csv \
    --t PAS-comp \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --ymax 31000 \
    --splitISM \
    -o ${OUTPLOTS}/GM12878_PacBio

```
<img align="center" width="400" src="plots/GM12878_PacBio_PAS-comp_support.png">

## Panel D: Percentage of TALON transcript models with a poly(A) motif identified at their 3' end (GM12878 ONT)
```bash
OUT=PAS-comp/ONT_GM12878
mkdir -p ${OUT}
python run_computational_PAS_analysis.py \
        --gtf ${ONT_GTF} \
        --genome ${GENOME} \
        --maxdist 35 \
        --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_polyA_motif.csv \
    --t PAS-comp \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --ymax 31000 \
    --splitISM \
    -o ${OUTPLOTS}/GM12878_ONT

```
<img align="center" width="400" src="plots/GM12878_ONT_PAS-comp_support.png">

## Panel E: Percentage of TALON transcript models with RNA-PET support for their 5'-3' end pair (GM12878 PacBio)

```bash
OUT=RNA-PET/PacBio_GM12878
mkdir -p ${OUT}
python run_RNA-PET_analysis.py \
    --gtf ${PB_GTF} \
    --rnapet ${RNAPET} \
    --maxdist 100 \
    --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_RNA-PET_results.csv \
    --t RNA-PET \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --ymax 31000 \
    --splitISM \
    -o ${OUTPLOTS}/GM12878_PacBio
```
<img align="center" width="400" src="plots/GM12878_PacBio_RNA-PET_support.png">


## Panel F: Percentage of TALON transcript models with RNA-PET support for their 5'-3' end pair (GM12878 ONT)
```bash
OUT=RNA-PET/ONT_GM12878
mkdir -p ${OUT}
python run_RNA-PET_analysis.py \
    --gtf ${ONT_GTF} \
    --rnapet ${RNAPET} \
    --maxdist 100 \
    --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_RNA-PET_results.csv \
    --t RNA-PET \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --ymax 31000 \
    --splitISM \
    -o ${OUTPLOTS}/GM12878_ONT
```
<img align="center" width="400" src="plots/GM12878_ONT_RNA-PET_support.png">
