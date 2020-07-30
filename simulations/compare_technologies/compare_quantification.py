import pandas as pd
import argparse
import scipy.stats as stats


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-sim_files', dest='sim_files',
		help='Comma-separated list of simulated fasta header files')
	parser.add_argument('-sim_name', dest='sim_name',
		help='Name of simulated dataset, ie "Perfect" or "Normal"')
	parser.add_argument('-t_ab', dest='t_ab_file',
		help='Abundance file used for transcript quantification, formatted like TALON abundance files')
	parser.add_argument('-g_ab', dest='g_ab_file',
		help='Abundance file used for gene quantification, formatted like TALON abundance files')
	parser.add_argument('-tech_name', dest='tech_name', 
		help='Name of technology analyzed, ie TALON or SQANTI')
	parser.add_argument('-ref_gtf', dest='ref_gtf', 
		help='Reference GTF to grab gene id to transcript id mappings')
	parser.add_argument('-merge_reps', dest='merge_reps',
		help='Consider replicates individually or separately. Default=False',
		action='store_true')
	parser.add_argument('-datasets', dest='datasets',
		help='Comma-separated list of dataset names corresponding '
		'to the datsets in the same order as sim_files')
	args = parser.parse_args()
	return args

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

# get a dictionary that maps transcript id to gene id
# used to determine which gene each simulated read 
# comes from 
def get_transcript_gene_map(ref_annot):

	# get transcript id to gene id map
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
	return transcript_gene_map

# calculate gene expression in tpm from an input
# abundance file
# if merge_reps=True, returns dataframe with abundance t_df
# if merge_reps=False, returns a list of dfs, one for each 
# dataset ie [g_df dataset1, g_df dataset2...]
def calc_gene_exp(ab_file, merge_reps, datasets):
	df = pd.read_csv(ab_file, sep='\t')
	df = df[['annot_gene_id']+datasets]

	if merge_reps:
		df['total_counts'] = df[datasets].sum(axis=1)
		df = df[['annot_gene_id', 'total_counts']].groupby(['annot_gene_id']).sum()
		df.reset_index(inplace=True)
		total_counts = df.total_counts.sum()
		df['tpm'] = (df.total_counts/total_counts)*1000000
		g_dfs = df
	else: 
		g_dfs = []
		for d in datasets: 
			temp = df.copy(deep=True)
			temp = temp[['annot_gene_id', d]].groupby(['annot_gene_id']).sum()
			temp.reset_index(inplace=True)
			total_counts = temp[d].sum()
			temp['tpm'] = (temp[d]/total_counts)*1000000
			temp = temp[['annot_gene_id', 'tpm']]
			g_dfs.append(temp)

	return g_dfs

# calculate transcript expression in tpm from an input
# abundance file
# if merge_reps=True, returns dataframe with abundance t_df
# if merge_reps=False, returns a list of dfs, one for each 
# dataset ie [t_df dataset1, t_df dataset2...]
def calc_transcript_exp(ab_file, merge_reps, datasets):
	df = pd.read_csv(ab_file, sep='\t')
	df = df[['annot_transcript_id', 'transcript_novelty']+datasets]

	# first, reformat transcript_id because nanosim 
	# inexplicably removes the '.#' suffix from transcript ids
	df['mini_tid'] = df.annot_transcript_id.str.extract(r'([A-z]+[0-9]+).',
		expand=False)

	# also account for flair stupidity
	df['base_tid'] = df.annot_transcript_id.str.extract(r'([A-z]+[0-9.]+)',
		expand=False)

	# remove novel things, because simulated data shouldn't have 
	# any novel things
	df = df.loc[df.transcript_novelty == 'Known']
	df.drop('annot_transcript_id', axis=1, inplace=True)

	if merge_reps:
		df['total_counts'] = df[datasets].sum(axis=1)
		df = df[['mini_tid', 'total_counts']].groupby(['mini_tid']).sum()
		df.reset_index(inplace=True)
		total_counts = df.total_counts.sum()
		df['tpm'] = (df.total_counts/total_counts)*1000000
		t_dfs = df
	else: 
		t_dfs = []
		for d in datasets: 
			temp = df.copy(deep=True)
			temp = temp[['mini_tid', d]].groupby(['mini_tid']).sum()
			temp.reset_index(inplace=True)
			total_counts = temp[d].sum()
			temp['tpm'] = (temp[d]/total_counts)*1000000
			t_dfs.append(temp)
	return t_dfs

# calculates gene and transcript counts from a processed 
# simulated data headers file
def calc_ctrl_tpms(df):
	t_df = df[['read_t', 'read_g']].groupby(['read_t']).count()
	t_df.rename({'read_g': 'counts'}, axis=1, inplace=True)
	t_df.reset_index(inplace=True)

	g_df = df[['read_g', 'read_t']].groupby(['read_g']).count()
	g_df.rename({'read_t': 'counts'}, axis=1, inplace=True)
	g_df.reset_index(inplace=True)

	# convert counts to tpm
	total_counts = t_df.counts.sum()
	t_df['tpm'] = (t_df.counts/total_counts)*1000000

	total_counts = g_df.counts.sum()
	g_df['tpm'] = (g_df.counts/total_counts)*1000000

	return g_df, t_df

# calculates the gene and transcript-level expression in TPM
# from the simulation read header files
# if merge_reps is true, will return a tuple of dataframes (g_df, t_df)
# if merge_reps is false, will return a tuple of lists of dataframes,
# in the order of the simulated files ie
# ([g_df from sim1, g_df from sim2... ],[t_df from sim1, t_df from sim2... ])
def get_control_exp(sim_files, t_g_map, merge_reps):

	if merge_reps: 
		df = pd.DataFrame(columns=['read_t', 'read_g'])
	else:
		g_dfs = []
		t_dfs = []

	# get the unique isoform ids from each simulated dataset 
	for f in sim_files:
		temp = pd.read_csv(f, header=None, names=['read_name'])
		temp['read_t'] = temp.read_name.str.split(pat='_', n=1, expand=True)[0]
		temp['read_t'] = temp.read_t.str.split(pat='>', n=1, expand=True)[1]	
		temp = temp['read_t'].to_frame()
		temp.columns = ['read_t']
		temp['read_g'] = temp.read_t.map(t_g_map, na_action='ignore')

		# add to the parent df if merging
		if merge_reps:
			df = pd.concat([df,temp], axis=0)
		# if not go ahead and calculate the tpm already
		else: 
			g_df, t_df = calc_ctrl_tpms(temp)
			g_dfs.append(g_df)
			t_dfs.append(t_df)

	# calculate the tpm after merging
	if merge_reps:
		g_dfs, t_dfs = calc_ctrl_tpms(df)

	return g_dfs, t_dfs

def main():
	args = get_args()
	sim_files = args.sim_files
	t_ab_file = args.t_ab_file
	g_ab_file = args.g_ab_file
	tech_name = args.tech_name
	sim_name = args.sim_name
	ref_gtf = args.ref_gtf
	merge_reps = args.merge_reps
	datasets = args.datasets

	# get the simulated data filenames
	sim_files = sim_files.split(',')
	datasets = datasets.split(',')

	# get the reference gid to tid map from the annotation gtf
	t_g_map = get_transcript_gene_map(ref_gtf)

	# get ground truth gene and transcript expression from the 
	# simulated data
	ctrl_g_dfs, ctrl_t_dfs = get_control_exp(sim_files, t_g_map, merge_reps)

	# get technology's quantification of genes and transcripts
	t_dfs = calc_transcript_exp(t_ab_file, merge_reps, datasets)
	g_dfs = calc_gene_exp(g_ab_file, merge_reps, datasets)

	if merge_reps:
		# transcript correlation
		t_df = ctrl_t_dfs.merge(t_dfs, how='left',
			left_on='read_t', right_on='mini_tid',
			suffixes=['_ctrl', '_test']).fillna(0)
		control = t_df['tpm_ctrl'].tolist()
		test = t_df['tpm_test'].tolist()
		t_corr, t_pval = stats.pearsonr(control, test)

		# gene correlation
		g_df = ctrl_g_dfs.merge(g_dfs, how='left', 
			left_on='read_g', right_on='annot_gene_id',
			suffixes=['_ctrl', '_test']).fillna(0)
		control = g_df['tpm_ctrl'].tolist()
		test = g_df['tpm_test'].tolist()
		g_corr, t_pval = stats.pearsonr(control, test)

		print('Simulated dataset: {}, Technology: {}'.format(sim_name, tech_name))
		print('Gene correlation: {}, gene pval: {}'.format(g_corr, g_pval))
		print('Transcript correlation: {}, transcript pval: {}'.format(t_corr, t_pval))

	else: 
		for d, ctrl_t_df, ctrl_g_df, test_t_df, test_g_df in zip(datasets,ctrl_t_dfs,ctrl_g_dfs,t_dfs,g_dfs):
			# transcript correlation
			t_df = ctrl_t_df.merge(test_t_df, how='left',
				left_on='read_t', right_on='mini_tid',
				suffixes=['_ctrl', '_test']).fillna(0)
			control = t_df['tpm_ctrl'].tolist()
			test = t_df['tpm_test'].tolist()
			t_corr, t_pval = stats.pearsonr(control, test)

			# gene correlation
			g_df = ctrl_g_df.merge(test_g_df, how='left', 
				left_on='read_g', right_on='annot_gene_id',
				suffixes=['_ctrl', '_test']).fillna(0)
			control = g_df['tpm_ctrl'].tolist()
			test = g_df['tpm_test'].tolist()
			g_corr, g_pval = stats.pearsonr(control, test)

			print('Simulated dataset: {}, {}, Technology: {}'.format(sim_name, d, tech_name))
			print('Gene correlation: {}, gene pval: {}'.format(g_corr, g_pval))
			print('Transcript correlation: {}, transcript pval: {}'.format(t_corr, t_pval))

if __name__ == '__main__':
	main()