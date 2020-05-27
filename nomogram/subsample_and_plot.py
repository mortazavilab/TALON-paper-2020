import sys
import pandas as pd
import argparse 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--f', dest='infile', 
		help='gene and transcript id tsv from talon read_annot file')
	# parser.add_argument('--intervals', dest='intervals',
	# 	help='CSV integer read number intervals at which to sample')
	# parser.add_argument('--whitelist', dest='whitelist', 
	# 	help='CSV whitelist file of transcripts that have passed filtering')
	parser.add_argument('--prefix', dest='prefix',
		help='prefix for saved plots.')
	args = parser.parse_args()
	return args

# subsample and compute tpm
def subsample_transcripts(df1, df2, n):
	sub1 = df1.sample(n, random_state=1)
	sub2 = df2.sample(n, random_state=1)

	# compute counts for each rep separately
	sub1 = sub1.groupby(['gene_ID','transcript_ID','gene_novelty','transcript_novelty'])['transcript_ID'].count().to_frame()
	sub1.rename({'transcript_ID': 'rep1_counts'}, axis=1, inplace=True)
	sub1.reset_index(inplace=True)
	sub2 = sub2.groupby(['gene_ID','transcript_ID','gene_novelty','transcript_novelty'])['transcript_ID'].count().to_frame()
	sub2.rename({'transcript_ID': 'rep2_counts'}, axis=1, inplace=True)
	sub2.reset_index(inplace=True)

	# merge the replicates
	sub_df = sub1.merge(sub2, how='outer', on=['gene_ID','transcript_ID','gene_novelty','transcript_novelty']).fillna(0)
	
	# the famous talon filter!
	sub_df = sub_df.loc[((sub_df.rep1_counts>5)&(sub_df.rep2_counts>5))|(sub_df.transcript_novelty=='Known')]

	# only nncs/nics that pass the reproducibility filter
	sub_df = sub_df.loc[(sub_df.transcript_novelty=='Known')|(sub_df.transcript_novelty=='NIC')|(sub_df.transcript_novelty=='NNC')]

	# sum counts across reps and compute TPM
	sub_df['counts'] = sub_df.rep1_counts+sub_df.rep2_counts
	total_count = sub_df.counts.sum()
	sub_df['tpm'] = (sub_df.counts*1000000)/total_count
	sub_df.drop(['rep1_counts', 'rep2_counts'], axis=1, inplace=True)
	return sub_df

def compute_total_t_tpm(df1, df2):

	# compute counts for each rep separately
	df1 = df1.groupby(['gene_ID','transcript_ID','gene_novelty','transcript_novelty'])['transcript_ID'].count().to_frame()
	df1.rename({'transcript_ID': 'rep1_counts'}, axis=1, inplace=True)
	df1.reset_index(inplace=True)
	df2 = df2.groupby(['gene_ID','transcript_ID','gene_novelty','transcript_novelty'])['transcript_ID'].count().to_frame()
	df2.rename({'transcript_ID': 'rep2_counts'}, axis=1, inplace=True)
	df2.reset_index(inplace=True)

	# merge the replicates
	ab = df1.merge(df2, how='outer', on=['gene_ID','transcript_ID','gene_novelty','transcript_novelty']).fillna(0)
	
	# the famous talon filter!
	ab = ab.loc[((ab.rep1_counts>5)&(ab.rep2_counts>5))|(ab.transcript_novelty=='Known')]

	# only nncs/nics that pass the reproducibility filter
	ab = ab.loc[(ab.transcript_novelty=='Known')|(ab.transcript_novelty=='NIC')|(ab.transcript_novelty=='NNC')]

	# sum counts across reps and compute TPM
	ab['counts'] = ab.rep1_counts+ab.rep2_counts
	total_count = ab.counts.sum()
	ab['ab_tpm'] = (ab.counts*1000000)/total_count

	# get only relevant fields
	ab = ab[['transcript_ID', 'ab_tpm']]

	return ab

# subsample and compute tpm
def subsample_genes(df1, df2, n):
	sub1 = df1.sample(n, random_state=1)
	sub2 = df2.sample(n, random_state=1)

	# concatenate the different reps
	sub_df = pd.concat([sub1, sub2])

	# filter here
	# only known genes
	sub_df = sub_df.loc[sub_df.gene_novelty == 'Known']

	sub_df = sub_df.groupby(['gene_ID'])['gene_ID'].count().to_frame()
	sub_df.rename({'gene_ID': 'counts'}, axis=1, inplace=True)
	sub_df.reset_index(inplace=True)
	total_count = sub_df.counts.sum()
	sub_df['tpm'] = (sub_df.counts*1000000)/total_count
	return sub_df

def compute_total_g_tpm(df):
	ab = df.copy(deep=True)

	# only known genes
	ab = ab.loc[ab.gene_novelty == 'Known']

	# compute TPM
	ab = ab.groupby(['gene_ID'])['gene_ID'].count().to_frame()
	ab.rename({'gene_ID': 'counts'}, axis=1, inplace=True)
	ab.reset_index(inplace=True)
	total_count = ab.counts.sum()
	ab['ab_tpm'] = (ab.counts*1000000)/total_count

	# get only relevant fields
	ab = ab[['gene_ID', 'ab_tpm']]

	return ab

def gene(df, df1, df2, intervals, prefix):
	# aggregate abundance and compute tpm
	ab = compute_total_g_tpm(df)

	# bin transcripts by expression level. How many transcripts belong to each bin?
	# assign each transcript a bin based on its expression level as well
	# bins = [(0,1),(1,2),(2,5),(5,10),(10,50),(50,100),(100,500),(500,ab.ab_tpm.max())]
	bins = [0,5,10,50,100,500,ab.ab_tpm.max()+1]
	ab_tpm = ab.ab_tpm.values.tolist()
	ab_bins = np.digitize(ab_tpm, bins)
	ab['bin'] = ab_bins
	ab['bin_total'] = ab['bin'].map(ab['bin'].value_counts())

	# remove entries from bin 1
	ab = ab.loc[ab.bin != 1] 

	# create a bin df to keep track of how many transcripts belong to each bin
	bin_df = ab[['bin', 'bin_total']].groupby(['bin', 'bin_total']).count()
	bin_df.reset_index(inplace=True)

	# loop through each interval and subsample df. 
	plot_data = pd.DataFrame(columns=['reads','bin','perc_within_10'])
	for n in intervals:
		sub_df = subsample_genes(df1, df2, n)
		sub_df = sub_df.merge(ab, how='inner', on='gene_ID').fillna(0)
		sub_df['within_5'] = (sub_df.tpm >= sub_df.ab_tpm*.9)&(sub_df.tpm <= sub_df.ab_tpm*1.1)

		temp = sub_df[['bin', 'within_5', 'tpm']].groupby(['bin', 'within_5']).count()
		temp.rename({'tpm':'bin_count'},axis=1,inplace=True)
		temp.reset_index(inplace=True)
		temp = temp.loc[temp.within_5 == True]

		for b in ab.bin.unique().tolist():
			if b not in temp.bin.tolist():
				temp = temp.append({'bin':b,'within_5':True,'bin_count':0}, ignore_index=True)

		temp = temp.merge(bin_df, how='left', on='bin')
		temp['perc_within_10'] = (temp.bin_count/temp.bin_total)*100
		temp['reads'] = n
		plot_data = pd.concat([plot_data,temp[['reads', 'bin', 'perc_within_10']]])

	# convert to human-readable TPM values
	plot_data['TPM bin'] = plot_data.apply(lambda x: (bins[x.bin-1],bins[x.bin]), axis=1)
	ax = sns.lineplot(x='reads', y='perc_within_10', hue='TPM bin', marker='o', data=plot_data)
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
	plt.savefig('{}_gene_nomogram.png'.format(prefix))
	plt.clf()	

def transcript(df1, df2, intervals, prefix):
	# aggregate abundance over both replicates and compute tpm
	# from parent abundance count
	# ab = pd.read_csv(file, '\t')
	# ab_cols = ab.columns[11:]
	# ab['counts'] = ab[ab_cols].sum(axis=1)

	# total_count = ab.counts.sum()
	# ab['ab_tpm'] = (ab.counts*1000000)/total_count
	# ab = ab[['transcript_ID', 'ab_tpm']]
	ab = compute_total_t_tpm(df1, df2)

	# bin transcripts by expression level. How many transcripts belong to each bin?
	# assign each transcript a bin based on its expression level as well
	# bins = [(0,1),(1,2),(2,5),(5,10),(10,50),(50,100),(100,500),(500,ab.ab_tpm.max())]
	bins = [0,5,10,50,100,500,ab.ab_tpm.max()+1]
	ab_tpm = ab.ab_tpm.values.tolist()
	ab_bins = np.digitize(ab_tpm, bins)
	ab['bin'] = ab_bins
	ab['bin_total'] = ab['bin'].map(ab['bin'].value_counts())

	# remove entries from bin 1
	ab = ab.loc[ab.bin != 1]

	# create a bin df to keep track of how many transcripts belong to each bin
	bin_df = ab[['bin', 'bin_total']].groupby(['bin', 'bin_total']).count()
	bin_df.reset_index(inplace=True)

	# loop through each interval and subsample df. 
	plot_data = pd.DataFrame(columns=['reads','bin','perc_within_10'])
	for n in intervals:
		sub_df = subsample_transcripts(df1, df2, n)
		sub_df = sub_df.merge(ab, how='inner', on='transcript_ID').fillna(0)
		sub_df['within_5'] = (sub_df.tpm >= sub_df.ab_tpm*.9)&(sub_df.tpm <= sub_df.ab_tpm*1.1)

		temp = sub_df[['bin', 'within_5', 'tpm']].groupby(['bin', 'within_5']).count()
		temp.rename({'tpm':'bin_count'},axis=1,inplace=True)
		temp.reset_index(inplace=True)
		temp = temp.loc[temp.within_5 == True]

		for b in ab.bin.unique().tolist():
			if b not in temp.bin.tolist():
				temp = temp.append({'bin':b,'within_5':True,'bin_count':0}, ignore_index=True)

		temp = temp.merge(bin_df, how='left', on='bin')
		temp['perc_within_10'] = (temp.bin_count/temp.bin_total)*100
		temp['reads'] = n
		plot_data = pd.concat([plot_data,temp[['reads', 'bin', 'perc_within_10']]])

	# convert to human-readable TPM values
	plot_data['TPM bin'] = plot_data.apply(lambda x: (bins[x.bin-1],bins[x.bin]), axis=1)
	ax = sns.lineplot(x='reads', y='perc_within_10', hue='TPM bin', marker='o', data=plot_data)
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
	plt.savefig('{}_transcript_nomogram.png'.format(prefix))
	plt.clf()

def main():
	args = get_args()
	df = pd.read_csv(args.infile, sep='\t')
	df1 = df.loc[df.dataset.str.contains('R1')].copy(deep=True)
	df2 = df.loc[df.dataset.str.contains('R2')].copy(deep=True)

	# remove novel genes and transcripts
	# df = df.loc[(df.gene_novelty == 'Known')&(df.transcript_novelty == 'Known')]

	# whitelist = pd.read_csv(args.whitelist)
	# whitelist.columns = ['gene_ID', 'transcript_ID']

	min_reads = min(len(df1.index), len(df2.index))
	intervals = [i for i in range(250000, min_reads, 250000)]
	intervals.append(min_reads)

	# # parse the different intervals
	# intervals = args.intervals.split(',')
	# intervals = [int(i) for i in intervals]
	# intervals.append(len(df.index))

	transcript(df1, df2, intervals, args.prefix)
	gene(df, df1, df2, intervals, args.prefix)

if __name__ == '__main__':
	main()