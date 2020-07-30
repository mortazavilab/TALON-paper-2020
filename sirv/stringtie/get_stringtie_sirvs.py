import pandas as pd
import argparse 

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-gtfs', dest='gtfs',
		help='CSV list of GTFs from stringtie with transcript isoforms')
	parser.add_argument('-ref_id_field', dest='ref_id_field', 
		help='field name where transcript id from the ref annot')
	args = parser.parse_args()
	return args

def process_st_gtf(f, ref_id_field):
	df = pd.read_csv(f,
		sep='\t', comment='#', usecols=[0,2,8],
		names=['chrom', 'entry_type', 'fields'])

	# we only care about transcripts
	df = df.loc[df.entry_type == 'transcript']
	df.drop('entry_type', axis=1, inplace=True)

	# we only care about SIRVs
	df = df.loc[df.chrom.str.contains('SIRV')]
	df.drop('chrom', axis=1, inplace=True)

	# # we only care about things that are expressed
	df['coverage'] = df.fields.str.extract(r'cov "([0-9.]+)', expand=False)
	df.coverage = df.coverage.astype('float64')
	len1 = len(df.index)
	df = df.loc[df.coverage != 0]
	len2 = len(df.index)
	df.drop('coverage', axis=1, inplace=True)

	if len1 - len2 != 0:
		print('Removed {} 0-coverage transcripts'.format(abs(len1-len2)))
	else:
		print('All transcripts present have above 0 coverage')

	# get the transcript id
	df['ref_id'] = df.fields.str.extract(r'{} "([A-z0-9.-]+)"'.format(ref_id_field),
		expand=False)
	df.ref_id.fillna('Novel', inplace=True)
	df.drop('fields', axis=1, inplace=True)

	return df

def main():
	args = get_args()
	gtfs = args.gtfs
	ref_id_field = args.ref_id_field

	gtfs = gtfs.split(',')
	
	first_df = True
	for g in gtfs:
		temp = process_st_gtf(g, ref_id_field)
		if first_df:
			df = temp
			first_df = False
		else:
			df = df.merge(temp, how='outer', on='ref_id')

	known = len(df.loc[df.ref_id.str.contains('SIRV')].index)
	novel = len(df.loc[~df.ref_id.str.contains('SIRV')].index)

	print('StringTie2 found {} known SIRV isoforms'.format(known))
	print('StringTie2 found {} novel SIRV isoforms'.format(novel))

if __name__ == '__main__':
	main()