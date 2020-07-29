import pandas as pd

def main():
	df = pd.read_csv('PB125_PB126_FLNC_stringtieMerge.out.gtf',
		sep='\t', comment='#', usecols=[0,2,8],
		names=['chrom', 'entry_type', 'fields'])

	# we only care about transcripts
	df = df.loc[df.entry_type == 'transcript']
	df.drop('entry_type', axis=1, inplace=True)

	# we only care about SIRVs
	df = df.loc[df.chrom.str.contains('SIRV')]
	df.drop('chrom', axis=1, inplace=True)

	# we only care about things that are expressed
	df['coverage'] = df.fields.str.extract(r'cov "([0-9.]+)', expand=False)
	df.coverage = df.coverage.astype('float64')
	df = df.loc[df.coverage != 0]
	df.drop('coverage', axis=1, inplace=True)

	# get the transcript id
	df['transcript_id'] = df.fields.str.extract(r'transcript_id "([A-z0-9.-]+)"',
		expand=False)

	# how many known transcripts do we have
	known = len(df.loc[~df.transcript_id.str.contains('MSTRG')].index)

	# how many novel transcripts do we have
	novel = len(df.loc[df.transcript_id.str.contains('MSTRG')].index)

	print('StringTie2 found {} known SIRV isoforms'.format(known))
	print('StringTie2 found {} novel SIRV isoforms'.format(novel))

if __name__ == '__main__':
	main()