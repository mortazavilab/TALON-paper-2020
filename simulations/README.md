# Simulations

We downloaded release 2.6 from https://github.com/bcgsc/NanoSim.git.

Format our Illumina expression data like the input for NanoSim
```bash
r1=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/Illumina/GM12878/Kallisto/Rep1/abundance.tsv
r2=/dfs3/pub/freese/mortazavi_lab/bin/TALON-paper-2020/Illumina/GM12878/Kallisto/Rep2/abundance.tsv
python create_nanosim_expression.py \
	-r1 $r1 \
	-r2 $r2 
```

Copy read data to current dir bc the hpc dumb af
```bash
ont_r1_fastq=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/ONT_RNA02_GM12878_R1/unmapped_reads/ONT45_6.fastq
ont_r2_fastq=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/data/ONT_RNA02_GM12878_R2/unmapped_reads/ONT47.fastq
cp $ont_r1_fastq ont_r1_fastq
cp $ont_r2_fastq ont_r2_fastq
```

Aggregate ONT reps (both fastqs and sams)
```bash
module load samtools
ont_r1_fastq=ont_r1.fastq
ont_r2_fastq=ont_r2.fastq
ont_r1_sam=/data/users/freese/TALON_data/revisions_1-20/data/ONT_RNA02_GM12878_R1/label_reads/ONT_Rep1_labeled.sam
ont_r2_sam=/data/users/freese/TALON_data/revisions_1-20/data/ONT_RNA02_GM12878_R2/label_reads/ONT_Rep2_labeled.sam
both_fastq=combined.fastq
cat $ont_r1_fastq > $both_fastq
cat $ont_r2_fastq >> $both_fastq

both_sam=combined.sam
samtools view -H $ont_r1_sam > $both_sam
samtools view $ont_r1_sam >> $both_sam
samtools view $ont_r2_sam >> $both_sam
```

For whatever reason, do this
```bash 
samtools view -H $both_sam > both_sam_back
samtools view $both_sam | grep -v = >> both_sam_back 
mv both_sam_back $both
```

For whatever other reason, do this
```bash 
ref_transcriptome=~/mortazavi_lab/ref/gencode.v29/gencode.v29.transcripts.fa 
python format_gencode_ref_transcriptome.py \
	-rt $ref_transcriptome \
	-ofile gencode.v29.transcripts.fa
```

Build the read profiles for the aggregate ONT data
```bash
conda activate nanosim
module load minimap2/2.17
both_fastq=combined.fastq
both_sam=combined.sam
ref_genome=~/mortazavi_lab/ref/hg38/hg38.fa
ref_transcriptome=gencode.v29.transcripts.fa
ns_path=/data/users/freese/mortazavi_lab/bin/NanoSim/src/

python ${ns_path}read_analysis.py transcriptome \
	-i $both_fastq \
	-rg $ref_genome \
	-rt $ref_transcriptome \
	-ga $both_sam \
	--no_intron_retention \
	-o ont_training \
	-t 16
```

## Simulate the reads!

Perfect reads
```bash
#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -N p_rep1

conda activate nanosim
ns_path=/data/users/freese/mortazavi_lab/bin/NanoSim2.6/NanoSim/src/
model=~/mortazavi_lab/bin/TALON-paper-2020/simulations/mortazavi_model/ont_training
ref_genome=~/mortazavi_lab/ref/hg38/hg38.fa
ref_transcriptome=~/mortazavi_lab/ref/gencode.v29/gencode.v29.transcripts.fa 
exp=~/mortazavi_lab/bin/TALON-paper-2020/simulations/expression_abundance.tsv

python ${ns_path}simulator.py transcriptome \
	-rt $ref_transcriptome \
	-e $exp \
	-c $model \
	-o rep1_perf/sim \
	-n 2000000 \
	--perfect \
	-b guppy \
	-r dRNA \
	-s 0.996 \
	--no_model_ir \
	--seed 1 \
	-t 16

#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -N p_rep2

conda activate nanosim
ns_path=/data/users/freese/mortazavi_lab/bin/NanoSim2.6/NanoSim/src/
model=~/mortazavi_lab/bin/TALON-paper-2020/simulations/mortazavi_model/ont_training
ref_genome=~/mortazavi_lab/ref/hg38/hg38.fa
ref_transcriptome=~/mortazavi_lab/ref/gencode.v29/gencode.v29.transcripts.fa 
exp=~/mortazavi_lab/bin/TALON-paper-2020/simulations/expression_abundance.tsv

python ${ns_path}simulator.py transcriptome \
	-rt $ref_transcriptome \
	-e $exp \
	-c $model \
	-o rep2_perf/sim \
	-n 2000000 \
	--perfect \
	-b guppy \
	-r dRNA \
	-s 0.996 \
	--no_model_ir \
	--seed 2 \
	-t 16
```

Not perfect reads
```bash
#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -N rep1

conda activate nanosim
ns_path=/data/users/freese/mortazavi_lab/bin/NanoSim2.6/NanoSim/src/
model=~/mortazavi_lab/bin/TALON-paper-2020/simulations/mortazavi_model/ont_training
ref_genome=~/mortazavi_lab/ref/hg38/hg38.fa
ref_transcriptome=~/mortazavi_lab/ref/gencode.v29/gencode.v29.transcripts.fa 
exp=~/mortazavi_lab/bin/TALON-paper-2020/simulations/expression_abundance.tsv

python ${ns_path}simulator.py transcriptome \
	-rt $ref_transcriptome \
	-e $exp \
	-c $model \
	-o rep1/sim \
	-n 2000000 \
	-b guppy \
	-r dRNA \
	-s 0.996 \
	--no_model_ir \
	--seed 1 \
	-t 16

#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -N rep2

conda activate nanosim
ns_path=/data/users/freese/mortazavi_lab/bin/NanoSim2.6/NanoSim/src/
model=~/mortazavi_lab/bin/TALON-paper-2020/simulations/mortazavi_model/ont_training
ref_genome=~/mortazavi_lab/ref/hg38/hg38.fa
ref_transcriptome=~/mortazavi_lab/ref/gencode.v29/gencode.v29.transcripts.fa 
exp=~/mortazavi_lab/bin/TALON-paper-2020/simulations/expression_abundance.tsv

python ${ns_path}simulator.py transcriptome \
	-rt $ref_transcriptome \
	-e $exp \
	-c $model \
	-o rep2/sim \
	-n 2000000 \
	-b guppy \
	-r dRNA \
	-s 0.996 \
	--no_model_ir \
	--seed 2 \
	-t 16
```

Concatenate unaligned and aligned reads into one fasta (non-perfect run only!)
```bash
r1_aligned=rep1/sim_aligned_reads.fasta
r1_unaligned=rep1/sim_unaligned_reads.fasta
r1_combined=rep1/rep1.fasta

cat $r1_aligned > $r1_combined
cat $r1_unaligned >> $r1_combined

r2_aligned=rep2/sim_aligned_reads.fasta
r2_unaligned=rep2/sim_unaligned_reads.fasta
r2_combined=rep2/rep2.fasta

cat $r2_aligned > $r2_combined
cat $r2_unaligned >> $r2_combined
```

Just rename the perfect files
```bash
r1=rep1_perf/sim_aligned_reads.fasta
mv $r1 rep1_perf/rep1_perf.fasta

r2=rep2_perf/sim_aligned_reads.fasta
mv $r2 rep2_perf/rep2_perf.fasta
```

<!-- Convert U to T before mapping
```bash
awk '{ gsub("U","T",$10); print $10 }' rep1/ont_reads.fasta
awk '{ gsub("U","T",$10); print $10 }' rep2/ont_reads.fasta
``` -->

Run minimap2 on the reads
```bash 
ref='~/mortazavi_lab/ref/hg38/hg38.fa'

# perfect
minimap2 \
	-ax splice \
	-uf \
	-k14 \
	--MD \
	--secondary=no \
	${ref} \
	rep1_perf/rep1_perf.fasta > rep1/rep1_perf.sam
minimap2 \
	-ax splice \
	-uf \
	-k14 \
	--MD \
	--secondary=no \
	${ref} \
	rep2_perf/rep2_perf.fasta > rep2_perf/rep2_perf.sam

# non perfect
minimap2 \
	-ax splice \
	-uf \
	-k14 \
	--MD \
	--secondary=no \
	${ref} \
	rep1/rep1.fasta > rep1/rep1.sam
minimap2 \
	-ax splice \
	-uf \
	-k14 \
	--MD \
	--secondary=no \
	${ref} \
	rep2/rep2.fasta > rep2/rep2.sam
```

Sort sams
```bash
module load samtools
# r1_dir=rep1/
# r1=rep1/rep1.sam
# r1_perf_dir=rep1_perf/
# r1_perf=rep1_perf/rep1_perf.sam
# r2_dir=rep2/
# r2=rep2/rep2.sam
# r2_perf_dir=rep2_perf/
# r2_perf=rep2_perf/rep2_perf.sam

r1_pref=rep1/rep1
r1_perf_pref=rep1_perf/rep1_perf/
r2_pref=rep2/rep2
r2_perf_pref=rep2_perf/rep2_perf

samtools view -h -b ${r1_pref}.sam > ${r1_pref}.bam
samtools sort ${r1_pref}.bam > ${r1_pref}_sorted.bam
samtools view ${r1_pref}_sorted.bam > ${r1_pref}_sorted.sam

samtools view -h -b ${r2_pref}.sam > ${r2_pref}.bam
samtools sort ${r2_pref}.bam > ${r2_pref}_sorted.bam
samtools view ${r2_pref}_sorted.bam > ${r2_pref}_sorted.sam

samtools view -h -b ${r1_perf_pref}.sam > ${r1_perf_pref}.bam
samtools sort ${r1_perf_pref}.bam > ${r1_perf_pref}_sorted.bam
samtools view ${r1_perf_pref}_sorted.bam > ${r1_perf_pref}_sorted.sam

samtools view -h -b ${r2_perf_pref}.sam > ${r2_perf_pref}.bam
samtools sort ${r2_perf_pref}.bam > ${r2_perf_pref}_sorted.bam
samtools view ${r2_perf_pref}_sorted.bam > ${r2_perf_pref}_sorted.sam
```

Run TranscriptClean on the reads
```bash
ref='~/mortazavi_lab/ref/hg38/hg38.fa'
sjs='/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/TC/gencode_v29_SIRV_SJs.tsv'
var='/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/refs/NA12878.vcf'
tc_path='~/mortazavi_lab/bin/TranscriptClean/'

python ${tc_path}TranscriptClean.py \
	-t 16 \
	--sam rep1/rep1_sorted.sam \
	--genome $ref \
	--spliceJns $sjs \
	--variants $var \
	--outprefix rep1/rep1 \
	--primaryOnly \
	--canonOnly

python ${tc_path}TranscriptClean.py \
	-t 16 \
	--sam rep2/rep2_sorted.sam \
	--genome $ref \
	--spliceJns $sjs \
	--variants $var \
	--outprefix rep2/rep2 \
	--primaryOnly \
	--canonOnly

python ${tc_path}TranscriptClean.py \
	-t 16 \
	--sam rep1_perf/rep1_perf_sorted.sam \
	--genome $ref \
	--spliceJns $sjs \
	--variants $var \
	--outprefix rep1_perf/rep1_perf \
	--primaryOnly \
	--canonOnly

python ${tc_path}TranscriptClean.py \
	-t 16 \
	--sam rep2_perf/rep2_perf_sorted.sam \
	--genome $ref \
	--spliceJns $sjs \
	--variants $var \
	--outprefix rep2_perf/rep2_perf \
	--primaryOnly \
	--canonOnly
```

Run talon_label_reads
```bash
talon_label_reads \
	--f rep1/nanosim_rep1_clean.sam \
	--g ~/mortazavi_lab/ref/hg38/hg38.fa \
	--t 16 \
	--tmpDir rep1/tmp_label_reads/ \
	--deleteTmp \
	--o rep1/nanosim

talon_label_reads \
	--f rep2/nanosim_rep2_clean.sam \
	--g ~/mortazavi_lab/ref/hg38/hg38.fa \
	--t 16 \
	--tmpDir rep2/tmp_label_reads/ \
	--deleteTmp \
	--o rep2/nanosim
```

## perf strand
```bash
rep1_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep1_perf_strand/Minimap2/"
rep2_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep2_perf_strand/Minimap2/"
r1_sam=${rep1_dir}ont_nanosim_sorted.sam
r2_sam=${rep2_dir}ont_nanosim_sorted.sam
r1_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_perf_strand/labeled/"
r2_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_perf_strand/labeled/"

mkdir $r1_dir
mkdir $r2_dir

talon_label_reads \
	--f ${r1_sam} \
	--g ~/mortazavi_lab/ref/hg38/hg38.fa \
	--t 16 \
	--tmpDir ${r1_dir}tmp_label_reads \
	--deleteTmp \
	--o ${r1_dir}ont_nanosim.sam

rep1_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep1_perf_strand/Minimap2/"
rep2_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep2_perf_strand/Minimap2/"
r1_sam=${rep1_dir}ont_nanosim_sorted.sam
r2_sam=${rep2_dir}ont_nanosim_sorted.sam
r1_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_perf_strand/labeled/"
r2_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_perf_strand/labeled/"

talon_label_reads \
	--f ${r2_sam} \
	--g ~/mortazavi_lab/ref/hg38/hg38.fa \
	--t 16 \
	--tmpDir ${r2_dir}tmp_label_reads \
	--deleteTmp \
	--o ${r2_dir}ont_nanosim.sam
```

## strand
```bash
rep1_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep1_strand/Minimap2/"
rep2_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep2_strand/Minimap2/"
r1_sam=${rep1_dir}ont_nanosim_sorted.sam
r2_sam=${rep2_dir}ont_nanosim_sorted.sam
r1_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_strand/labeled/"
r2_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_strand/labeled/"
 
mkdir $r1_dir
mkdir $r2_dir

talon_label_reads \
	--f ${r1_sam} \
	--g ~/mortazavi_lab/ref/hg38/hg38.fa \
	--t 16 \
	--tmpDir ${r1_dir}tmp_label_reads \
	--deleteTmp \
	--o ${r1_dir}ont_nanosim.sam

rep1_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep1_strand/Minimap2/"
rep2_dir="/data/users/freese/TALON_data/revisions_1-20/nanosim/rep2_strand/Minimap2/"
r1_sam=${rep1_dir}ont_nanosim_sorted.sam
r2_sam=${rep2_dir}ont_nanosim_sorted.sam
r1_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep1_strand/labeled/"
r2_dir="~/mortazavi_lab/bin/TALON-paper-2020/simulations/rep2_strand/labeled/"

talon_label_reads \
	--f ${r2_sam} \
	--g ~/mortazavi_lab/ref/hg38/hg38.fa \
	--t 16 \
	--tmpDir ${r2_dir}tmp_label_reads \
	--deleteTmp \
	--o ${r2_dir}ont_nanosim.sam
```


Create talon config file

## perf strand
```bash
printf "ont_1,Simulated ONT,ONT,rep1_perf_strand/nanosim_labeled.sam\n" > ps_config.csv
printf "ont_2,Simulated ONT,ONT,rep2_perf_strand/nanosim_labeled.sam" >> ps_config.csv
```

## strand
```bash
printf "ont_1,Simulated ONT,ONT,rep1_strand/nanosim_labeled.sam\n" > s_config.csv
printf "ont_2,Simulated ONT,ONT,rep2_strand/nanosim_labeled.sam" >> s_config.csv
```


```bash
printf "ont_1,Simulated ONT,ONT,rep1/nanosim_labeled.sam\n" > s_config.csv
printf "ont_2,Simulated ONT,ONT,rep2/nanosim_labeled.sam" >> s_config.csv
```

Init talon db from gencode v29
```bash
ref='/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'
talon_initialize_database \
	--f $ref \
	--g hg38 \
	--a gencode_v29 \
	--o talon
```

## perfect strand
```bash
ref='/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'
mkdir perfect_strand
talon_initialize_database \
	--f $ref \
	--g hg38 \
	--a gencode_v29 \
	--o perfect_strand/talon
```

## strand
```bash
ref='/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'
mkdir strand
talon_initialize_database \
	--f $ref \
	--g hg38 \
	--a gencode_v29 \
	--o strand/talon
```

Run talon
```bash
talon \
	--f config.csv \
	--db talon.db \
	--build hg38 \
	--t 16 \
	--o sim
```

## perf strand
```bash
talon \
	--f ps_config.csv \
	--db perfect_strand/talon.db \
	--build hg38 \
	--t 16 \
	--o perfect_strand/sim
```

## strand
```bash
talon \
	--f s_config.csv \
	--db strand/talon.db \
	--build hg38 \
	--t 16 \
	--o sim
```

Perform filtering using TALON

```bash
db=talon/simulated.db
talon_filter_transcripts \
	--db $db \
	-a v29 \
	--datasets rep1_perf,rep2_perf \
	--o talon/perf_whitelist.csv
talon_filter_transcripts \
	--db $db \
	-a v29 \
	--datasets rep1,rep2 \
	--o talon/whitelist.csv
```

Create an unfiltered abundance file for each set of datasets
```bash
printf "rep1_perf\nrep2_perf" > perf_datasets
talon_abundance \
	--db $db \
	-a v29 \
	-b hg38 \
	-d perf_datasets \
	--o talon/perf
```

```bash
printf "rep1\nrep2" > datasets
talon_abundance \
	--db $db \
	-a v29 \
	-b hg38 \
	-d datasets \
	--o talon/
```

Create a filtered abundance file for each set of datasets
```bash
talon_abundance \
	--db $db \
	-a v29 \
	-b hg38 \
	-d perf_datasets \
	--o talon/perf \
	--whitelist talon/perf_whitelist.csv
```

```bash
talon_abundance \
	--db $db \
	-a v29 \
	-b hg38 \
	-d datasets \
	--o talon/ \
	--whitelist talon/whitelist.csv
```

Isolate the fasta headers for quantification accuracy
```bash
r1=rep1/rep1.fasta
r2=rep2/rep2.fasta
r1_perf=rep1_perf/rep1_perf.fasta
r2_perf=rep2_perf/rep2_perf.fasta

grep \> $r1 > rep1/rep1_headers.fasta
grep \> $r2 > rep2/rep2_headers.fasta
grep \> $r1_perf > rep1_perf/rep1_perf_headers.fasta
grep \> $r2_perf > rep2_perf/rep2_perf_headers.fasta
```

Get quantification correlations
```bash
python compare_quantification.py 
```

```
rep1
Gene correlation: 0.9752182655068521, gene pval: 0.0
Transcript correlation: 0.839616947752433, transcript pval: 0.0

rep2
Gene correlation: 0.9793314342913847, gene pval: 0.0
Transcript correlation: 0.7571022861073621, transcript pval: 0.0

rep1_perf
Gene correlation: 0.9942027438596236, gene pval: 0.0
Transcript correlation: 0.9475228329665103, transcript pval: 0.0

rep2_perf
Gene correlation: 0.9929908876334156, gene pval: 0.0
Transcript correlation: 0.9367575707601811, transcript pval: 0.0
```

