import pandas as pd
import argparse 

def process_st_gtf(f):
df = pd.read_csv(f,
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

df = df['transcript_id'].to_frame()

	return df

def main():
	fs = ['PB125_abundanceMerged_stringtie.out.gtf', 'PB126_abundanceMerged_stringtie.out.gtf']
	
	df1 = process_st_gtf(fs[0])
	df2 = process_st_gtf(fs[1])

	df = df1.merge(df2, how='outer', on='transcript_id')
	# df.to_csv('all_expressed_stringtie_sirvs.tsv', sep='\t')

	known = len(df.loc[df.transcript_id.str.contains('SIRV')].index)
	novel = len(df.loc[~df.transcript_id.str.contains('SIRV')].index)

	print('StringTie2 found {} known SIRV isoforms'.format(known))
	print('StringTie2 found {} novel SIRV isoforms'.format(novel))

	known_bois = df.loc[~df.transcript_id.str.contains('MSTRG'), 'transcript_id']
	known_bois.to_csv('known_transcripts_for_gaby.tsv', sep='\t')

	novel_bois = df.loc[df.transcript_id.str.contains('MSTRG'), 'transcript_id']
	novel_bois.to_csv('novel_transcripts_for_gaby.tsv', sep='\t')	

	# # why does sqanti find fewer fsms than stringtie claims to report?
	# df = pd.read_csv('PB125_PB126_FLNC_stringtieMerge.out.gtf',
	# 	sep='\t', comment='#',
	# 	names=['chrom', 'dumb1','entry_type', 'start', 'stop',
	# 		   'dumb2', 'strand', 'dumb3', 'fields'])
	# df.drop(['dumb1', 'dumb2', 'dumb3'], axis=1, inplace=True)

	# # only look at SIRVs
	# df = df.loc[df.chrom.str.contains('SIRV')]

	# # we only care about things that are expressed
	# df['coverage'] = df.fields.str.extract(r'cov "([0-9.]+)', expand=False)
	# df.coverage = df.coverage.astype('float64')
	# df = df.loc[df.coverage != 0]
	# df.drop('coverage', axis=1, inplace=True)

	# # get the transcript id
	# df['transcript_id'] = df.fields.str.extract(r'transcript_id "([A-z0-9.-]+)"',
	# 	expand=False)

	# # remove novel transcripts
	# df = df.loc[~df.transcript_id.str.contains('MSTRG')]

	# # only keep transcripts I guess
	# df = df.loc[df.entry_type == 'transcript']

	# # do the same to the annotation GTF
	# annot_file = '/Users/fairliereese/mortazavi_lab/ref/gencode.v29/gencode.v29.SIRV.ERCC.annotation.gtf'
	# annot = pd.read_csv(annot_file,
	# 	sep='\t', comment='#', 
	# 	names=['chrom', 'dumb1','entry_type', 'start', 'stop',
	# 		   'dumb2', 'strand', 'dumb3', 'fields'])
	# annot.drop(['dumb1', 'dumb2', 'dumb3'], axis=1, inplace=True)

	# # only look at SIRVs
	# annot = annot.loc[annot.chrom.str.contains('SIRV')]

	# # only keep transcripts I guess
	# annot = annot.loc[annot.entry_type == 'transcript']

	# # get transcript ids
	# annot['transcript_id'] = annot.fields.str.extract(r'transcript_id "([A-z0-9.-]+)"',
	# 	expand=False)



if __name__ == '__main__':
	main()