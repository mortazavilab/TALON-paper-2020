import pandas as pd
import argparse 
import math

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-gtf', dest='gtf',
		help='GTF from stringtie with transcript isoforms')
	parser.add_argument('-ref_tid_field', dest='ref_tid_field', 
		help='field name where transcript id from the ref annot')
	parser.add_argument('-ref_gid_field', dest='ref_gid_field',
		help='field name where gene id from the ref annot')
	args = parser.parse_args()
	return args

def process_st_gtf(f, ref_id_field, no_cov_filt):
	df = pd.read_csv(f,
		sep='\t', comment='#', usecols=[0,2,8],
		names=['chrom', 'entry_type', 'fields'])

	# remove sirvs/erccs
	df = df[~df.chrom.str.contains('SIRV')]
	df = df[~df.chrom.str.contains('ERCC')]

	# we only care about transcripts
	df = df.loc[df.entry_type == 'transcript']
	df.drop('entry_type', axis=1, inplace=True)

	# get the transcript id
	df['ref_id'] = df.fields.str.extract(r'{} "([A-z0-9.-]+)"'.format(ref_id_field),
		expand=False)
	df.ref_id.fillna('Novel', inplace=True)
	df.drop('fields', axis=1, inplace=True)

	# get the gene id

	return df

def main():
	args = get_args()
	gtf = args.gtf
	ref_id_field = args.ref_id_field
	no_cov_filt = args.no_cov_filt

	df = process_st_gtf(gtf, ref_id_field, no_cov_filt)

	known = len(df.loc[df.ref_id.str.contains('SIRV')].index)
	novel = len(df.loc[~df.ref_id.str.contains('SIRV')].index)

	print('StringTie2 found {} known SIRV isoforms'.format(known))
	print('StringTie2 found {} novel SIRV isoforms'.format(novel))

if __name__ == '__main__':
	main()