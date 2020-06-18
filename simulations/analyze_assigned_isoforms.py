import pandas as pd
from collections import defaultdict

def get_fields(fields):
    attributes = {}
    description = fields.strip()
    description = [x.strip() for x in description.split(";")]
    for pair in description:
        if pair == "": continue
        pair = pair.replace('"', '')
        key, val = pair.split()
        attributes[key] = val
    # put in placeholders for important attributes (such as gene_id) if they
    # are absent
    if 'gene_id' not in attributes:
        attributes['gene_id'] = 'NULL'
    return attributes  

fname = 'talon/simulated_Gm12878_talon_read_annot.tsv'
ref_annot = '/Users/fairliereese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'
transcript_gene_map = {}
with open(ref_annot, 'r') as f:
	for line in f:
		if not line.startswith('#'):
			line = line.strip().split('\t')
			if line[2] == 'transcript':
				fields = line[-1]
				fields = get_fields(fields)
				gid = fields['gene_id']
				tid = fields['transcript_id'].split('.')[0]
				transcript_gene_map[tid] = gid

df = pd.read_csv(fname, sep='\t')

df = df[['read_name', 'dataset', 'transcript_ID', 'annot_transcript_id', 'annot_gene_id', 'transcript_novelty']]
df['read_gt_t'] = df.read_name.str.split(pat='_', n=1, expand=True)[0]
df['read_t'] = df.annot_transcript_id.str.split(pat='.', n=1, expand=True)[0]
df['read_gt_g'] = df.read_gt_t.map(transcript_gene_map, na_action='ignore')

df['same_gene'] = df.annot_gene_id == df.read_gt_g
df['same_transcript'] = df.read_t == df.read_gt_t

# print(df[['read_gt_t', 'read_gt_g', 'annot_gene_id', 'transcript_novelty', 'same_gene']].head(40))
# exit()

# which reads were simulated in the wrong orientation
df['read_sense'] = df.read_name.str.split(pat='_', expand=True)[4]

def display_results(n_df, p_df):
	total_reads = len(n_df.index)
	total_correct_t = len(n_df.loc[n_df.same_transcript == True].index)
	total_correct_g = len(n_df.loc[n_df.same_gene == True].index)
	print('Normal dataset has {} total reads'.format(total_reads))
	print('Normal dataset has {} correctly assigned transcripts'.format(total_correct_t))
	print('Normal dataset has {} correctly assigned genes'.format(total_correct_g))
	print('{}% of reads were assigned to the correct transcript'.format((total_correct_t/total_reads)*100))
	print('{}% of reads were assigned to the correct gene'.format((total_correct_g/total_reads)*100))
	total_reads = len(p_df.index)
	total_correct_t = len(p_df.loc[p_df.same_transcript == True].index)
	total_correct_g = len(p_df.loc[p_df.same_gene == True].index)
	print('Perfect dataset has {} total reads'.format(total_reads))
	print('Perfect dataset has {} correctly assigned transcripts'.format(total_correct_t))
	print('Perfect dataset has {} correctly assigned genes'.format(total_correct_g))
	print('{}% of reads were assigned to the correct transcript'.format((total_correct_t/total_reads)*100))
	print('{}% of reads were assigned to the correct gene'.format((total_correct_g/total_reads)*100))


def no_filter(df, perf, norm):
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	print('Unfiltered')
	display_results(n_df, p_df)


def ism_filter(df, perf, norm):
	df = df.loc[df.transcript_novelty != 'ISM']
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	print('Without TALON ISMs')
	display_results(n_df, p_df)


def antisense_filter(df, perf, norm):
	df = df.loc[df.read_sense != 'R']
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	print('Without NanoSim antisense reads')
	display_results(n_df, p_df)


def talon_filter(df, perf, norm):
	# use talon filter to remove reads assigned to a transcript that was filtered out
	wfile = 'talon/whitelist.csv'
	w_df = pd.read_csv(wfile, header=None, names=['gene_id', 'transcript_id'])
	whitelist = w_df.transcript_id.tolist()
	df = df.loc[df.transcript_ID.isin(whitelist)]
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	print('Without reads that did not pass talon filter')
	display_results(n_df, p_df)


def strict_talon_filter(df, perf, norm):
	# use talon filter to remove reads assigned to a transcript that was filtered out
	wfile = 'talon/whitelist.csv'
	w_df = pd.read_csv(wfile, header=None, names=['gene_id', 'transcript_id'])
	whitelist = w_df.transcript_id.tolist()
	df = df.loc[df.transcript_ID.isin(whitelist)]
	df = df.loc[df.transcript_novelty.isin(['Known', 'NIC', 'NNC'])]
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	print('Without reads that did not pass talon filter, only Known, NIC, NNC')
	display_results(n_df, p_df)


def t_dist_by_novelty(df, perf, norm):
	# only known transcripts can be correctly assigned
	df = df.loc[df.transcript_novelty == 'Known']
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	p_dist = p_df[['transcript_novelty', 'same_transcript',
			'transcript_ID']].groupby(['transcript_novelty',
				'same_transcript']).count()
	p_dist.rename({'transcript_ID':'counts', 
		'same_transcript': 'correct_transcript'}, axis=1, inplace=True)
	p_dist.reset_index(inplace=True)
	print('Perfect correct transcript assignment by novelty category')
	print(p_dist)
	p_dist['total'] = p_dist.apply(lambda x: \
		p_dist.loc[p_dist.transcript_novelty == x.transcript_novelty, 'counts'].sum(), axis=1)
	p_dist = p_dist.loc[p_dist.same_transcript == True]
	p_dist['percent_correct'] = p_dist.apply(lambda x: \
		(x.counts/x.total)*100, axis=1)
	print(p_dist)
	n_dist = n_df[['transcript_novelty', 'same_transcript',
			'transcript_ID']].groupby(['transcript_novelty',
				'same_transcript']).count()
	n_dist.rename({'transcript_ID': 'counts',
		'same_transcript': 'correct_transcript'}, axis=1, inplace=True)
	n_dist.reset_index(inplace=True)
	print('Normal correct transcript assignment by novelty category')
	print(n_dist)	
	n_dist['total'] = n_dist.apply(lambda x: \
		n_dist.loc[n_dist.transcript_novelty == x.transcript_novelty, 'counts'].sum(), axis=1)
	n_dist = n_dist.loc[n_dist.same_transcript == True]
	n_dist['percent_correct'] = n_dist.apply(lambda x: \
		(x.counts/x.total)*100, axis=1)
	print(n_dist)

def g_dist_by_novelty(df, perf, norm):
	p_df = df.loc[df.dataset.isin(perf_datasets)]
	n_df = df.loc[df.dataset.isin(norm_datasets)]
	print()
	p_dist = p_df[['transcript_novelty', 'same_gene',
			'transcript_ID']].groupby(['transcript_novelty',
				'same_gene']).count()
	p_dist.rename({'transcript_ID':'counts', 
		'same_gene': 'correct_transcript'}, axis=1, inplace=True)
	p_dist.reset_index(inplace=True)
	print('Perfect correct gene assignment by novelty category')
	print(p_dist)
	p_dist['total'] = p_dist.apply(lambda x: \
		p_dist.loc[p_dist.transcript_novelty == x.transcript_novelty, 'counts'].sum(), axis=1)
	p_dist = p_dist.loc[p_dist.same_gene == True]
	p_dist['percent_correct'] = p_dist.apply(lambda x: \
		(x.counts/x.total)*100, axis=1)
	print(p_dist)
	n_dist = n_df[['transcript_novelty', 'same_gene',
			'transcript_ID']].groupby(['transcript_novelty',
				'same_gene']).count()
	n_dist.rename({'transcript_ID': 'counts',
		'same_gene': 'correct_transcript'}, axis=1, inplace=True)
	n_dist.reset_index(inplace=True)
	print('Normal correct gene assignment by novelty category')
	print(n_dist)	
	n_dist['total'] = n_dist.apply(lambda x: \
		n_dist.loc[n_dist.transcript_novelty == x.transcript_novelty, 'counts'].sum(), axis=1)
	n_dist = n_dist.loc[n_dist.same_gene == True]
	n_dist['percent_correct'] = n_dist.apply(lambda x: \
		(x.counts/x.total)*100, axis=1)
	print(n_dist)

perf_datasets = ['rep1_perf', 'rep2_perf']
norm_datasets = ['rep1', 'rep2']

no_filter(df, perf_datasets, norm_datasets)
# how do these numbers change when accounting for ISM reads
ism_filter(df, perf_datasets, norm_datasets)
# and how do these numbers change when accounting for reverse strand reads
antisense_filter(df, perf_datasets, norm_datasets)
# what about when using reads from transcripts that pass the talon filter?
talon_filter(df, perf_datasets, norm_datasets)
# what about the talon filter + only known, nic, nnc?
strict_talon_filter(df, perf_datasets, norm_datasets)

# what novelty type categories are associated with the incorrectly versus 
# correctly assiged reads? Make a counts table.
t_dist_by_novelty(df, perf_datasets, norm_datasets)
g_dist_by_novelty(df, perf_datasets, norm_datasets)


