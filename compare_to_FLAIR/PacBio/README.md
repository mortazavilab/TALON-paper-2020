# Running FLAIR on GM12878 PacBio data

FLAIR was cloned from https://github.com/BrooksLabUCSC/flair on 8/5/2019.

# Running FLAIR on GM12878 PacBio data
FLAIR version: 1.4

1. Run align and correct steps separately on replicates
```bash
qsub PacBio_GM12878_1/./run_FLAIR_align.sh
qsub PacBio_GM12878_2/./run_FLAIR_align.sh
```
```bash
qsub PacBio_GM12878_1/./run_FLAIR_correct.sh
qsub PacBio_GM12878_2/./run_FLAIR_correct.sh
```
2. Then, run collapse step on concatenated files from both reps.
```bash
cat PacBio_GM12878_1/flair_all_corrected.psl PacBio_GM12878_2/flair_all_corrected.psl > PacBio_GM12878_1-PacBio_GM12878_2_flair_all_corrected.psl
cat PacBio_GM12878_1/ENCFF281TNJ.fastq PacBio_GM12878_2/ENCFF475ORL.fastq > PacBio_GM12878_1-PacBio_GM12878_2-concat.fastq
qsub ./run_flair_collapse.sh
```
3. Finally, run quantify step. To do this, you need to create a tab-delimited config file with fields dataset name, condition, batch, and fastq reads file. This is what the GM12878 file looks like:
```bash
GM12878_Rep1    GM12878    batch1    /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/PacBio_GM12878_1/ENCFF281TNJ.fastq
GM12878_Rep2    GM12878    batch1    /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/PacBio_GM12878_2/ENCFF475ORL.fastq
```
```bash
qsub ./run_flair_quantify.sh
```
