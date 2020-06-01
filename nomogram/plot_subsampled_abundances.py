import os
import sys
import pandas as pd
import argparse 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import glob

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-dir', dest='indir', 
		help='directory to parse for subsampled shit')
	parser.add_argument('--prefix', dest='prefix',
		help='prefix for saved plots.')
	parser.add_argument('--filter', dest='filter_type',
		help="type of filter to use. can be 'known' or 'talon'")
	parser.add_argument('--max_reads', dest='maxfile',
		help='Max number of reads file')
	args = parser.parse_args()
	return args

def compute_t_tpm(df, full=False):
	counts_cols = [col for col in df.columns if 'sam' in col]
	df['counts'] = df[counts_cols].sum(axis=1)
	df.drop(counts_cols, axis=1, inplace=True)
	total_count = df.counts.sum()
	df['tpm'] = (df.counts*1000000)/total_count

	df = df[['transcript_ID', 'tpm']]

	if full:
		df.rename({'tpm': 'ab_tpm'}, axis=1, inplace=True)

	return df

def compute_g_tpm(df, full=False):
	counts_cols = [col for col in df.columns if 'sam' in col]
	df['counts'] = df[counts_cols].sum(axis=1)
	counts_cols = counts_cols+['transcript_ID', 'transcript_novelty']
	df.drop(counts_cols, axis=1, inplace=True)

	# gene groupby
	gb_cols = ['gene_ID', 'gene_novelty']
	df = df.groupby(gb_cols)['counts'].agg('sum').to_frame()
	df.reset_index(inplace=True)

	total_count = df.counts.sum()
	df['tpm'] = (df.counts*1000000)/total_count

	df = df[['gene_ID', 'tpm']]

	if full:
		df.rename({'tpm': 'ab_tpm'}, axis=1, inplace=True)

	return df

def filter(f_type, g_or_t, df):
	if g_or_t == 'gene':
		df = df.loc[df.gene_novelty == 'Known']
	elif g_or_t == 'transcript':
		if f_type == 'known':
			df = df.loc[df.transcript_novelty == 'Known']
		elif f_type == 'talon':
			df = df.loc[df.transcript_novelty.isin(['Known', 'NNC', 'NIC'])]
	return df

def get_read_num(elem):
	try:
		n = int(elem.split('_')[2])
	except:
		n = 8000000
	return n

def gene(ab, ufilt_files, read_nums, f_type, max_reads, prefix):

	# bin transcripts by expression level. How many transcripts belong to each bin?
	# assign each transcript a bin based on its expression level as well
	# bins = [(0,1),(1,2),(2,5),(5,10),(10,50),(50,100),(100,500),(500,ab.ab_tpm.max())]
	bins = [5,10,50,100,500,ab.ab_tpm.max()+1]
	ab_tpm = ab.ab_tpm.values.tolist()
	ab_bins = np.digitize(ab_tpm, bins)
	ab['bin'] = ab_bins
	ab['bin_total'] = ab['bin'].map(ab['bin'].value_counts())

	# remove entries from bin 0
	ab = ab.loc[ab.bin != 0]

	# create a bin df to keep track of how many transcripts belong to each bin
	bin_df = ab[['bin', 'bin_total']].groupby(['bin', 'bin_total']).count()
	bin_df.reset_index(inplace=True)

	# loop through each interval and analyze subsampled data
	plot_data = pd.DataFrame(columns=['reads','bin','perc_within_10'])
	for fname,n in zip(ufilt_files, read_nums):

		# grab abundances, filter transcripts, and compute TPM
		sub_df = pd.read_csv(fname,sep='\t',usecols=[0,1,8,9,11,12])
		sub_df = filter(f_type, 'gene', sub_df)
		sub_df = compute_g_tpm(sub_df)

		# merge with full and determine if each gene is within 10%
		# of tpm val calculated with all reads
		sub_df = sub_df.merge(ab, how='right', on='gene_ID').fillna(0)
		sub_df['within_10'] = (sub_df.tpm >= sub_df.ab_tpm*.9)&(sub_df.tpm <= sub_df.ab_tpm*1.1)
		sub_df = sub_df[['bin', 'within_10', 'tpm']].groupby(['bin','within_10']).count()
		sub_df.rename({'tpm':'bin_count'}, axis=1, inplace=True)
		sub_df.reset_index(inplace=True)
		sub_df = sub_df.loc[sub_df.within_10 == True]

		# if there are missing True entries, fill them in so we can still calculate
		for b in ab.bin.unique().tolist():
			if b not in sub_df.bin.tolist():
				sub_df = sub_df.append({'bin':b,'within_10':True,'bin_count':0},
					ignore_index=True)

		# merge with total bin numbers so we can calculate percentages
		sub_df = sub_df.merge(bin_df, how='right', on='bin')
		sub_df['perc_within_10'] = (sub_df.bin_count/sub_df.bin_total)*100

		# set up some more metadata for plotting purposes
		if n == 'full':
			sub_df['reads'] = max_reads
		else:
			sub_df['reads'] = n
		plot_data = pd.concat([plot_data,sub_df[['reads','bin','perc_within_10']]])

	# convert to human-readable TPM values
	plot_data['TPM bin'] = plot_data.apply(lambda x: (bins[x.bin-1],bins[x.bin]), axis=1)
	ax = sns.lineplot(x='reads', y='perc_within_10', hue='TPM bin', marker='o', data=plot_data)
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
	plt.savefig('{}_gene_nomogram.png'.format(prefix))
	plt.clf()	

def transcript(ab, filt_files, read_nums, f_type, max_reads, prefix):

	# bin transcripts by expression level. How many transcripts belong to each bin?
	# assign each transcript a bin based on its expression level as well
	# bins = [(0,1),(1,2),(2,5),(5,10),(10,50),(50,100),(100,500),(500,ab.ab_tpm.max())]
	bins = [5,10,50,100,500,ab.ab_tpm.max()+1]
	ab_tpm = ab.ab_tpm.values.tolist()
	ab_bins = np.digitize(ab_tpm, bins)
	ab['bin'] = ab_bins
	ab['bin_total'] = ab['bin'].map(ab['bin'].value_counts())

	# remove entries from bin 0
	ab = ab.loc[ab.bin != 0]

	# create a bin df to keep track of how many transcripts belong to each bin
	bin_df = ab[['bin', 'bin_total']].groupby(['bin', 'bin_total']).count()
	bin_df.reset_index(inplace=True)

	# loop through each interval and analyze subsampled data
	plot_data = pd.DataFrame(columns=['reads','bin','perc_within_10'])
	for fname,n in zip(filt_files, read_nums):

		# grab abundances, filter transcripts, and compute TPM
		sub_df = pd.read_csv(fname,sep='\t',usecols=[0,1,8,9,11,12])
		sub_df = filter(f_type, 'transcript', sub_df)
		sub_df = compute_t_tpm(sub_df)

		# merge with full and determine if each transcript is within 10%
		# of tpm val calculated with all reads
		sub_df = sub_df.merge(ab, how='right', on='transcript_ID').fillna(0)
		sub_df['within_10'] = (sub_df.tpm >= sub_df.ab_tpm*.9)&(sub_df.tpm <= sub_df.ab_tpm*1.1)
		sub_df = sub_df[['bin', 'within_10', 'tpm']].groupby(['bin','within_10']).count()
		sub_df.rename({'tpm':'bin_count'}, axis=1, inplace=True)
		sub_df.reset_index(inplace=True)
		sub_df = sub_df.loc[sub_df.within_10 == True]

		# if there are missing True entries, fill them in so we can still calculate
		for b in ab.bin.unique().tolist():
			if b not in sub_df.bin.tolist():
				sub_df = sub_df.append({'bin':b,'within_10':True,'bin_count':0},
					ignore_index=True)

		# merge with total bin numbers so we can calculate percentages
		sub_df = sub_df.merge(bin_df, how='right', on='bin')
		sub_df['perc_within_10'] = (sub_df.bin_count/sub_df.bin_total)*100

		# set up some more metadata for plotting purposes
		if n == 'full':
			sub_df['reads'] = max_reads
		else:
			sub_df['reads'] = n
		plot_data = pd.concat([plot_data,sub_df[['reads','bin','perc_within_10']]])

	# convert to human-readable TPM values
	plot_data['TPM bin'] = plot_data.apply(lambda x: (bins[x.bin-1],bins[x.bin]), axis=1)
	ax = sns.lineplot(x='reads', y='perc_within_10', hue='TPM bin', marker='o', data=plot_data)
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
	plt.savefig('{}_transcript_nomogram.png'.format(prefix))
	plt.clf()

def main():

	args = get_args()
	indir = args.indir
	prefix = args.prefix
	f_type = args.filter_type
	maxfile = args.maxfile

	# get max reads
	with open(maxfile, 'r') as beep:
		for line in beep:
			max_reads = int(line.strip())
			max_reads = int(max_reads/2)

	# get each sam file in the input dir
	filt_files = []
	for file in glob.glob(indir+'*filtered.tsv'):
		filt_files.append(file)
	filt_files.sort(key=get_read_num)

	# get each sam file in the input dir
	ufilt_files = []
	for file in glob.glob(indir+'*abundance.tsv'):
		ufilt_files.append(file)
	ufilt_files.sort(key=get_read_num)

	# parse out the read numbers from each thing
	read_nums = [int(i.split('_')[2]) for i in filt_files[:-1]]
	read_nums.append('full')

	# load the full abundance files to compute ground-truth TPM vals
	g_df = pd.read_csv(ufilt_files[-1], sep='\t', usecols=[0,1,8,9,11,12])
	t_df = pd.read_csv(filt_files[-1], sep='\t', usecols=[0,1,8,9,11,12])

	# filter things according to whichever filter we've decided on
	g_df = filter(f_type, 'gene', g_df)
	t_df = filter(f_type, 'transcript', t_df)

	# compute tpms
	g_df = compute_g_tpm(g_df, full=True)
	t_df = compute_t_tpm(t_df, full=True)

	# plot de plot
	gene(g_df, ufilt_files, read_nums, f_type, max_reads, prefix)
	transcript(t_df, filt_files, read_nums, f_type, max_reads, prefix)


if __name__ == '__main__':
	main()