import argparse
import pandas as pd 

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r1', dest='rep1', 
		help='Illumina Rep1 Kallisto expression')
	parser.add_argument('-r2', dest='rep2', 
		help='Illumina Rep1 Kallisto expression')
	args = parser.parse_args()
	return args

def read_illumina_ab(fname):
	df1 = pd.read_csv(fname, sep='|', usecols=[0])
	df1.columns = ['target_id']
	df2 = pd.read_csv(fname, sep='\t', usecols=[3,4])
	df = pd.concat([df1, df2], axis=1)
	return df

def combine_abs(rep1, rep2):
	df = rep1.merge(rep2, on='target_id')
	df['tpm'] = (df.tpm_x+df.tpm_y)/2
	df['est_counts'] = df.est_counts_x+df.est_counts_y
	df = df[['target_id', 'est_counts', 'tpm']]
	df = df.loc[df.est_counts > 0]
	return df

def main():
	args = get_args()
	df1 = read_illumina_ab(args.rep1)
	df2 = read_illumina_ab(args.rep2)
	df = combine_abs(df1, df2)

	df.to_csv('expression_abundance.tsv', sep='\t', index=False)

if __name__ == '__main__':
	main()