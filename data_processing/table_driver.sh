
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

# create a filtered GTF
talon_create_GTF \
	--db $db \
	-b hg38_SIRV \
	--a gencode_v29_SIRV \
	--whitelist pb_whitelist.csv \
	--observed \
	-d pb_datasets \
	--o pb

## ONT ##
printf "ONT_GM12878_R1\nONT_GM12878_R2" > ont_datasets

# unfiltered PB abundance
talon_abundance \
	--db $db \
	--a gencode_v29_SIRV \
	-b hg38_SIRV \
	-d ont_datasets \
	--o  ont

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

# create a filtered GTF
talon_create_GTF \
	--db $db \
	-b hg38_SIRV \
	--a gencode_v29_SIRV \
	--whitelist ont_whitelist.csv \
	--observed \
	-d ont_datasets \
	--o ont

