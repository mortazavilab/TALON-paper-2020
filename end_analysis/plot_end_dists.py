import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import argparse 
import seaborn as sns

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', dest='infile', 
		help='Read annot file')
	parser.add_argument('-g', dest='gtf',
		help='Annotation GTF')
	parser.add_argument('--datasets', dest='datasets',
		help='CSV list of datasets to analyze', default=None)
	parser.add_argument('-o', dest='oprefix',
		help='Output plot prefix')
	args = parser.parse_args()
	return args

def main():
	args = get_args()
	infile = args.infile
	datasets = args.datasets
	oprefix = args.oprefix
	gtf = args.gtf

	df = pd.read_csv(infile, sep='\t')

	# only use requested datasets
	if datasets:
		beep = datasets.split(',')
		df = df.loc[df.dataset.isin(beep)]

	# only use Known transcripts
	df = df.loc[df.transcript_novelty == 'Known']

	# remove SIRVs and ERCCs
	df = df.loc[~df.chrom.str.contains('SIRV')]
	df = df.loc[~df.chrom.str.contains('ERCC')]


	# pair it down only to the columns we need
	df = df[['read_name', 'read_start', 'read_end', 'annot_transcript_id', 'strand']]
	df.rename({'annot_transcript_id':'transcript_id'}, axis=1, inplace=True)

	# read in gtf 
	gtf_df = pd.read_csv(gtf, sep='\t', comment='#',
		names=['feature_type', 'annot_start', 'annot_end', 'fields'],
		usecols=[2,3,4,8])
	gtf_df = gtf_df.loc[gtf_df.feature_type == 'transcript']
	gtf_df['transcript_id'] = gtf_df.fields.str.extract(r'transcript_id "([A-z]+[0-9.]+)', expand=False)
	# gtf_df[['beep', 'bop']] = gtf_df.fields.str.split(pat='transcript_id "', expand=True, n=1)
	# gtf_df[['transcript_id', 'plop']] = gtf_df.bop.str.split(pat='";', expand=True, n=1)
	gtf_df.drop(['fields', 'feature_type'], axis=1, inplace=True)

	# merge gtf with the read annot file 
	df = df.merge(gtf_df, how='left', on='transcript_id')

	# calculate distance between annotated and sequenced ends
	df['tss_diff'] = np.nan
	df['tes_diff'] = np.nan

	# forward strand
	df.loc[df.strand == '+', 'tss_diff'] = df['read_start'] - df['annot_start']
	df.loc[df.strand == '+', 'tes_diff'] = df['read_end'] - df['annot_end']

	# reverse strand 
	df.loc[df.strand == '-', 'tss_diff'] = df['annot_start'] - df['read_start']
	df.loc[df.strand == '-', 'tes_diff'] = df['annot_end'] - df['read_end']

	# plot de plot
	bins = [i for i in range(-800,801)]
	ax = sns.distplot(df.tss_diff, kde=False, bins=bins)
	ax.set(xlabel='Distance from annotated TSS (bp)')
	ax.set(ylabel='Number of reads')
	ax.set(xlim=(-800,800))
	ax.set(ylim=(0,75000))
	plt.savefig('{}_tss_dists.png'.format(oprefix))
	plt.clf()

	ax = sns.distplot(df.tes_diff, kde=False, bins=bins)
	ax.set(xlabel='Distance from annotated TES (bp)')
	ax.set(ylabel='Number of reads')
	ax.set(xlim=(-800,800))
	ax.set(ylim=(0,75000))
	plt.savefig('{}_tes_dists.png'.format(oprefix))
	plt.clf()

if __name__ == '__main__':
	main()