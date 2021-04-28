import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import argparse 
import seaborn as sns

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', dest='infile', 
		help='Read annot file')
	parser.add_argument('-a', dest='annot_ends',
		help='Ends from the reference annotation')
	parser.add_argument('--datasets', dest='datasets',
		help='CSV list of datasets to analyze', default=None)
	parser.add_argument('-xlim', dest='xlim',
		help='X limit to serve as boundary for all plot')
	parser.add_argument('-ylim', dest='ylim',
		help='Y limit for the all plot')
	parser.add_argument('-o', dest='oprefix',
		help='Output plot prefix')
	args = parser.parse_args()
	return args

def plot_hires_bins(df, end_type, oprefix):
	diff_col = '{}_diff'.format(end_type)
	hires_bin_col = '{}_hires_bin'.format(end_type)
	perc_col = 'perc_{}_hires'.format(end_type)

	hires_total = len(df.loc[(df[diff_col] > -500)&(df[diff_col] <= 500)].index)
	hires_bins = [i for i in range(-500,0,50)]
	hires_bins += [-1, 0, 1]
	hires_bins += [i for i in range(50,550,50)]

	df[hires_bin_col] = pd.cut(df[diff_col], bins=hires_bins)
	hires_df = df[[hires_bin_col, diff_col]].groupby(hires_bin_col).count()
	hires_df.rename({diff_col: 'counts'}, axis=1, inplace=True)
	hires_df.reset_index(inplace=True)
	hires_df[perc_col] = (hires_df.counts/hires_total)*100
	hires_df.counts.fillna(0, inplace=True)
	hires_df[perc_col].fillna(0, inplace=True)

	ax = sns.barplot(x=hires_bin_col, y=perc_col,
			data=hires_df, color='#009E73', saturation=1)
	ax.set(xlabel='Distance from annotated {} (bp)'.format(end_type.upper()))
	ax.set(ylabel='Percentage of Known Reads within 500 bp')
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
	ax.set(ylim=(0,100))
	fname = '{}_{}_dists_hires.png'.format(oprefix, end_type)
	print('Saving plot {}'.format(fname))
	plt.savefig(fname,
		bbox_inches='tight')
	plt.clf()

# def plot_all_bins(df, end_type, oprefix):
# 	diff_col = '{}_diff'.format(end_type)
# 	bin_col = '{}_bin'.format(end_type)
# 	perc_col = 'perc_{}'.format(end_type)

# 	total = len(df.index)
# 	bins = [i for i in range(-500,0,50)]
# 	bins += [-1, 0, 1]
# 	bins += [i for i in range(50,550,50)]
# 	bins += [df[diff_col].max()]
# 	bins = [df[diff_col].min()] + bins

# 	df[bin_col] = pd.cut(df[diff_col], bins=bins)
# 	bin_df = df[[bin_col, diff_col]].groupby(bin_col).count()
# 	bin_df.rename({diff_col: 'counts'}, axis=1, inplace=True)
# 	bin_df.reset_index(inplace=True)
# 	bin_df[perc_col] = (bin_df.counts/total)*100
# 	bin_df.counts.fillna(0, inplace=True)
# 	bin_df[perc_col].fillna(0, inplace=True)

# 	ax = sns.barplot(x=bin_col, y=perc_col,
# 			data=bin_df, color='#009E73', saturation=1)
# 	ax.set(xlabel='Distance from annotated {} (bp)'.format(end_type.upper()))
# 	ax.set(ylabel='Percentage of Known Reads')
# 	ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
# 	ax.set(ylim=(0,100))
# 	fname = '{}_{}_dists.png'.format(oprefix, end_type)
# 	print('Saving plot {}'.format(fname))
# 	plt.savefig(fname,
# 		bbox_inches='tight')
# 	plt.clf()

def plot_all_bins(df, end_type, xlim, ylim, oprefix):
	diff_col = '{}_diff'.format(end_type)

	bins = [i for i in range(-1*xlim,xlim+1)]
	ax = sns.distplot(df[diff_col], kde=False, bins=bins,
		color='#009E73')
	ax.set(xlabel='Distance from annotated {} (bp)'.format(end_type.upper()))
	ax.set(ylabel='Number of reads')
	ax.set(xlim=(-1*xlim,xlim))
	# ax.set(ylim=(0,ylim))
	fname = '{}_{}_dists.png'.format(oprefix, end_type)
	print('Saving plot {}'.format(fname))
	plt.savefig(fname,
		bbox_inches='tight')
	plt.clf()

def main():
	args = get_args()
	infile = args.infile
	datasets = args.datasets
	oprefix = args.oprefix
	annot_ends = args.annot_ends
	xlim = int(args.xlim)
	ylim = int(args.ylim)

	df = pd.read_csv(infile, sep='\t')

	# only use requested datasets
	if datasets:
		beep = datasets.split(',')
		df = df.loc[df.dataset.isin(beep)]

	# only use Known transcripts
	df = df.loc[df.transcript_novelty == 'Known']

	# remove SIRVs and ERCCs and chrM
	df = df.loc[~df.chrom.str.contains('SIRV')]
	df = df.loc[~df.chrom.str.contains('ERCC')]
	df = df.loc[df.chrom != 'chrM']

	# only look at multiexonic things, which are transcripts
	# we are more confident about
	df = df.loc[df.n_exons > 1]

	# pair it down only to the columns we need
	df = df[['read_name', 'read_start', 'read_end', 'annot_transcript_id', 'strand']]

	# read in annotation ends 
	annot_df = pd.read_csv(annot_ends, sep='\t')

	# Merge together reads with GENCODE start data
	df = pd.merge(df, annot_df, on='annot_transcript_id', 
					how='left')

	# Compute TSS/TES distances
	# Negative distance is upstream, positive is downstream
	df['tss_diff'] = df.read_start - df.TSS_pos
	df['tes_diff'] = df.read_end - df.TES_pos
	df.loc[df['strand']=='-', ['tss_diff', 'tes_diff']] *= -1

	# plot de plot
	plot_hires_bins(df, 'tss', oprefix)
	plot_hires_bins(df, 'tes', oprefix)
	plot_all_bins(df, 'tss', xlim, ylim, oprefix)
	plot_all_bins(df, 'tes', xlim, ylim, oprefix)

if __name__ == '__main__':
	main()

