
odir=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/

cd $odir

db=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_1-20/human_TALON/full_gencode_v29_SIRV_2020-02-29.db

## PACBIO ##

printf "PB_GM12878_R1\nPB_GM12878_R2" > pb_datasets

# unfiltered PB abundance
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d pb_datasets \
	--o  pb
cp pb_talon_abundance.tsv backup/
grep -v 'SIRV\|gSpikein_ERCC' pb_talon_abundance.tsv > pb_talon_abundance_no_sirv_no_ercc.tsv
head -1 pb_talon_abundance.tsv > pb_sirv_talon_abundance.tsv
grep SIRV pb_talon_abundance.tsv >> pb_sirv_talon_abundance.tsv
head -1 pb_talon_abundance.tsv > pb_ercc_talon_abundance.tsv
grep gSpikein_ERCC pb_talon_abundance.tsv >> pb_ercc_talon_abundance.tsv
mv pb_talon_abundance_no_sirv_no_ercc.tsv pb_talon_abundance.tsv

# create a whitelist for pb datasets
talon_filter_transcripts \
	--db $db \
	-a gencode_v29_SIRV \
	--maxFracA 0.5 \
	--o pb_whitelist.csv \
	--datasets pb_datasets

# create a filtered abundance file
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d pb_datasets \
	--whitelist pb_whitelist.csv \
	--o pb
cp pb_talon_abundance_filtered.tsv backup/
grep -v 'SIRV\|gSpikein_ERCC' pb_talon_abundance_filtered.tsv > pb_talon_abundance_filtered_no_sirv_no_ercc.tsv
head -1 pb_talon_abundance_filtered.tsv > pb_sirv_talon_abundance_filtered.tsv
grep SIRV pb_talon_abundance_filtered.tsv >> pb_sirv_talon_abundance_filtered.tsv
head -1 pb_talon_abundance_filtered.tsv > pb_ercc_talon_abundance_filtered.tsv
grep gSpikein_ERCC pb_talon_abundance_filtered.tsv >> pb_ercc_talon_abundance_filtered.tsv
mv pb_talon_abundance_filtered_no_sirv_no_ercc.tsv pb_talon_abundance_filtered.tsv

# create a filtered GTF
talon_create_GTF \
	--db $db \
	-b hg38_SIRV \
	--a gencode_v29_SIRV \
	--whitelist pb_whitelist.csv \
	-d pb_datasets \
	--o pb
cp pb_talon.gtf backup/
grep -v 'SIRV\|gSpikein_ERCC' pb_talon.gtf > pb_talon_no_sirv_no_ercc.gtf
grep SIRV pb_talon.gtf > pb_sirv_talon.gtf
grep gSpikein_ERCC pb_talon.gtf > pb_ercc_talon.gtf
mv pb_talon_no_sirv_no_ercc.gtf pb_talon.gtf


# use the talon_longest_ends module to assign transcript starts and ends
gtf=pb_talon.gtf
annot=GM12878_talon_read_annot.tsv
datasets=pb_datasets
talon_longest_end \
    -gtf $gtf \
    -read_annot $annot \
    --datasets $datasets \
    --mode both \
    --novelty all \
    -o pb_talon

## ONT ##

printf "ONT_GM12878_R1\nONT_GM12878_R2" > ont_datasets

# unfiltered PB abundance
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d ont_datasets \
	--o  ont
cp ont_talon_abundance.tsv backup/
grep -v 'SIRV\|gSpikein_ERCC' ont_talon_abundance.tsv > ont_talon_abundance_no_sirv_no_ercc.tsv
head -1 ont_talon_abundance.tsv > ont_sirv_talon_abundance.tsv
grep SIRV ont_talon_abundance.tsv >> ont_sirv_talon_abundance.tsv
head -1 ont_talon_abundance.tsv > ont_ercc_talon_abundance.tsv
grep gSpikein_ERCC ont_talon_abundance.tsv >> ont_ercc_talon_abundance.tsv
mv ont_talon_abundance_no_sirv_no_ercc.tsv ont_talon_abundance.tsv

# create a whitelist for ont datasets
talon_filter_transcripts \
	--db $db \
	-a gencode_v29_SIRV \
	--maxFracA 0.5 \
	--o ont_whitelist.csv \
	--datasets ont_datasets

# create a filtered abundance file
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d ont_datasets \
	--whitelist ont_whitelist.csv \
	--o ont
cp ont_talon_abundance_filtered.tsv backup/
grep -v 'SIRV\|gSpikein_ERCC' ont_talon_abundance_filtered.tsv > ont_talon_abundance_filtered_no_sirv_no_ercc.tsv
head -1 ont_talon_abundance_filtered.tsv > ont_sirv_talon_abundance_filtered.tsv
grep SIRV ont_talon_abundance_filtered.tsv >> ont_sirv_talon_abundance_filtered.tsv
head -1 ont_talon_abundance_filtered.tsv > ont_ercc_talon_abundance_filtered.tsv
grep gSpikein_ERCC ont_talon_abundance_filtered.tsv >> ont_ercc_talon_abundance_filtered.tsv
mv ont_talon_abundance_filtered_no_sirv_no_ercc.tsv ont_talon_abundance_filtered.tsv

# create a filtered GTF
talon_create_GTF \
	--db $db \
	-b hg38_SIRV \
	--a gencode_v29_SIRV \
	--whitelist ont_whitelist.csv \
	-d ont_datasets \
	--o ont
cp ont_talon.gtf backup/
grep -v 'SIRV\|gSpikein_ERCC' ont_talon.gtf > ont_talon_no_sirv_no_ercc.gtf
grep SIRV ont_talon.gtf > ont_sirv_talon.gtf
grep gSpikein_ERCC ont_talon.gtf > ont_ercc_talon.gtf
mv ont_talon_no_sirv_no_ercc.gtf ont_talon.gtf

# use the talon_longest_ends module to assign transcript starts and ends
gtf=ont_talon.gtf
annot=GM12878_talon_read_annot.tsv
datasets=pb_datasets
talon_longest_end \
    -gtf $gtf \
    -read_annot $annot \
    --datasets $datasets \
    --mode both \
    --novelty all \
    -o pb_talon

## PACBIO + ONT ##

printf "PB_GM12878_R1\nPB_GM12878_R2\nONT_GM12878_R1\nONT_GM12878_R2" > pb_ont_datasets

# unfiltered PB+ONT abundance
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d pb_ont_datasets \
	--o  pb_ont
cp pb_ont_talon_abundance.tsv backup/
grep -v 'SIRV\|gSpikein_ERCC' pb_ont_talon_abundance.tsv > pb_ont_talon_abundance_no_sirv_no_ercc.tsv
head -1 pb_ont_talon_abundance.tsv > pb_ont_sirv_talon_abundance.tsv
grep SIRV pb_ont_talon_abundance.tsv >> pb_ont_sirv_talon_abundance.tsv
head -1 pb_ont_talon_abundance.tsv > pb_ont_ercc_talon_abundance.tsv
grep gSpikein_ERCC pb_ont_talon_abundance.tsv >> pb_ont_ercc_talon_abundance.tsv
mv pb_ont_talon_abundance_no_sirv_no_ercc.tsv pb_ont_talon_abundance.tsv

# create a whitelist comprised of unique values in both the ONT and PB whitelists
cat pb_whitelist.csv > pb_ont_whitelist.csv 
cat ont_whitelist.csv >> pb_ont_whitelist.csv
cat pb_ont_whitelist.csv | sort | uniq > pb_ont_whitelist.csv

# create a filtered abundance file
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d pb_ont_datasets \
	--whitelist pb_ont_whitelist.csv \
	--o pb_ont
cp pb_ont_talon_abundance_filtered.tsv backup/
grep -v 'SIRV\|gSpikein_ERCC' pb_ont_talon_abundance_filtered.tsv > pb_ont_talon_abundance_filtered_no_sirv_no_ercc.tsv
head -1 pb_ont_talon_abundance_filtered.tsv > pb_ont_sirv_talon_abundance_filtered.tsv
grep SIRV pb_ont_talon_abundance_filtered.tsv >> pb_ont_sirv_talon_abundance_filtered.tsv
head -1 pb_ont_talon_abundance_filtered.tsv > pb_ont_ercc_talon_abundance_filtered.tsv
grep gSpikein_ERCC pb_ont_talon_abundance_filtered.tsv >> pb_ont_ercc_talon_abundance_filtered.tsv
mv pb_ont_talon_abundance_filtered_no_sirv_no_ercc.tsv pb_ont_talon_abundance_filtered.tsv

# create a filtered GTF
talon_create_GTF \
	--db $db \
	-b hg38_SIRV \
	--a gencode_v29_SIRV \
	--whitelist pb_ont_whitelist.csv \
	-d pb_ont_datasets \
	--o pb_ont
cp pb_ont_talon.gtf backup/
grep -v 'SIRV\|gSpikein_ERCC' pb_ont_talon.gtf > pb_ont_talon_no_sirv_no_ercc.gtf
grep SIRV pb_ont_talon.gtf > pb_ont_sirv_talon.gtf
grep gSpikein_ERCC pb_ont_talon.gtf > pb_ont_ercc_talon.gtf
mv pb_ont_talon_no_sirv_no_ercc.gtf pb_ont_talon.gtf



