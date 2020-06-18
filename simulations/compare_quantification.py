import pandas as pd
import scipy.stats as stats

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

# get transcript id to gene id map
ref_annot = '/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'
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

ufilt_ab = 'talon/talon_abundance.tsv'
ufilt_df = pd.read_csv(ufilt_ab, sep='\t', usecols=[2,3,11,12])
filt_ab = 'talon/talon_abundance_filtered.tsv'
filt_df = pd.read_csv(filt_ab, sep='\t', usecols=[2,3,11,12])
perf_ufilt_ab = 'talon/perf_talon_abundance.tsv'
perf_ufilt_df = pd.read_csv(perf_ufilt_ab, sep='\t', usecols=[2,3,11,12])
perf_filt_ab = 'talon/perf_talon_abundance_filtered.tsv'
perf_filt_df = pd.read_csv(perf_filt_ab, sep='\t', usecols=[2,3,11,12])

def calc_tpm(df, r1, r2):
	# r1 tpm
	total_count = df[r1].sum()
	df['{}_tpm'.format(r1)] = (df[r1]/total_count)*1000000
	# r2 tpm
	total_count = df[r2].sum()
	df['{}_tpm'.format(r2)] = (df[r2]/total_count)*1000000
	return df 

def calc_transcript_tpm(df, r1, r2):
	df = df.copy(deep=True)
	# remove .# from each tid
	df['tid'] = df.annot_transcript_id.str.split(pat='.', n=1, expand=True)[0]
	df = calc_tpm(df, r1, r2)
	return df

def calc_gene_tpm(df, r1, r2):
	df = df.copy(deep=True)
	# groupby and sum to get gene expression
	df = df[['annot_gene_id', r1, r2]].groupby('annot_gene_id').sum()
	df.reset_index(inplace=True)
	df = calc_tpm(df, r1, r2)
	return df

# calculate gene tpm
talon_g_df = calc_gene_tpm(ufilt_df, 'rep1', 'rep2')
talon_perf_g_df = calc_gene_tpm(perf_ufilt_df, 'rep1_perf', 'rep2_perf')

# calculate transcript tpm
talon_t_df = calc_transcript_tpm(filt_df, 'rep1', 'rep2')
talon_perf_t_df = calc_transcript_tpm(perf_filt_df, 'rep1_perf', 'rep2_perf')

rep1 = 'rep1/rep1_headers.fasta'
rep2 = 'rep2/rep2_headers.fasta'
rep1_perf = 'rep1_perf/rep1_perf_headers.fasta'
rep2_perf = 'rep2_perf/rep2_perf_headers.fasta'


def compare_reps(control_headers, talon_t_df, talon_g_df, rep):
	print()
	print('Comparing {}'.format(rep))
	df = pd.read_csv(control_headers, header=None, names=['read_name'])
	df['read_t'] = df.read_name.str.split(pat='_', n=1, expand=True)[0]
	df['read_t'] = df.read_t.str.split(pat='>', n=1, expand=True)[1]
	df['read_g'] = df.read_t.map(transcript_gene_map, na_action='ignore')

	# calculate gene and transcript counts
	t_df = df[['read_t', 'read_g']].groupby(['read_t']).count()
	t_df.rename({'read_g': 'counts'}, axis=1, inplace=True)
	t_df.reset_index(inplace=True)
	g_df = df[['read_g', 'read_t']].groupby(['read_g']).count()
	g_df.rename({'read_t': 'counts'}, axis=1, inplace=True)
	g_df.reset_index(inplace=True)

	# calculate the tpm
	total_counts = t_df.counts.sum()
	t_df['tpm'] = (t_df.counts/total_counts)*1000000
	total_counts = g_df.counts.sum()
	g_df['tpm'] = (g_df.counts/total_counts)*1000000

	# merge the genes and correlate genes
	g_df = g_df.merge(talon_g_df, how='left', left_on='read_g', right_on='annot_gene_id').fillna(0)
	control = g_df['tpm'].tolist()
	talon = g_df['{}_tpm'.format(rep)].tolist()
	g_corr, g_pval = stats.pearsonr(control, talon) 

	# merge the genes and correlate transcripts
	t_df = t_df.merge(talon_t_df, how='left', left_on='read_t', right_on='tid').fillna(0)
	control = t_df['tpm'].tolist()
	talon = t_df['{}_tpm'.format(rep)].tolist()
	t_corr, t_pval = stats.pearsonr(control, talon) 
	
	# print stuff
	print('Gene correlation: {}, gene pval: {}'.format(g_corr, g_pval))
	print('Transcript correlation: {}, transcript pval: {}'.format(t_corr, t_pval))


compare_reps(rep1, talon_t_df, talon_g_df, 'rep1')
compare_reps(rep2, talon_t_df, talon_g_df, 'rep2')
compare_reps(rep1_perf, talon_perf_t_df, talon_perf_g_df, 'rep1_perf')
compare_reps(rep2_perf, talon_perf_t_df, talon_perf_g_df, 'rep2_perf')

