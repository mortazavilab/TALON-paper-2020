import pandas as pd
import scipy.stats as stats
import argparse

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-sim_files', dest='sim_files',
		help='Comma-separated list of simulated fasta header files')
	parser.add_argument('-sim_name', dest='sim_name',
		help='Name of simulated dataset, ie "Perfect" or "Normal"')
	parser.add_argument('-ab', dest='ab_file',
		help='Abundance file, formatted like TALON abundance files')
	parser.add_argument('-tech_name', dest='tech_name', 
		help='Name of technology analyzed, ie TALON or SQANTI')
	parser.add_argument('-merge_reps', dest='merge_reps',
		help='Consider replicates individually or separately. Default=False',
		action='store_true')
	parser.add_argument('-datasets', dest='datasets',
		help='Comma-separated list of dataset names corresponding '
		'to the datsets in the same order as sim_files')
	args = parser.parse_args()
	return args

# if we're merging returns the number of unique isoforms 
# present across the simulated files. 
# if not, returns a list of the number of isoforms present 
# in each simulated file ie [num isoforms rep 1, num isoforms rep 2]
def count_total_isoforms(sim_files, merge_reps):

	if merge_reps:
		df = pd.DataFrame(columns=['read_t'])
	else:
		num_isoforms = []

	# get the unique isoform ids from each simulated dataset 
	for f in sim_files:
		temp = pd.read_csv(f, header=None, names=['read_name'])
		temp['read_t'] = temp.read_name.str.split(pat='_', n=1, expand=True)[0]
		temp['read_t'] = temp.read_t.str.split(pat='>', n=1, expand=True)[1]	
		temp = temp['read_t'].to_frame()
		temp.drop_duplicates(inplace=True)
		temp.columns = ['read_t']

		if merge_reps:
			# add to the parent df
			df = pd.concat([df,temp], axis=0)
		else:
			num_isoforms.append(len(temp.read_t.unique()))

	# get the number of isoforms from the merged table 
	if merge_reps:
		# make sure to remove duplicates that may have arisen 
		# across replicates
		df.drop_duplicates(inplace=True)
		df.sort_values(by='read_t', inplace=True)

		# count and return the number of unique isoforms that should
		# be present in the dataset
		num_isoforms = len(df.read_t.unique())

	return num_isoforms

# if we're merging replicates, returns a tuple of integers of the
# (num known, num novel) isoforms detected. if not merging, 
# returns a tuple of lists corresponding to each dataset of the
# ([num known dataset 1, num known dataset 2...],
# [num novel dataset 1, num novel dataset 2 ...])
def count_detected_isoforms(ab_file, datasets, merge_reps):

	# if merging
	if merge_reps:
		# known detections
		df = pd.read_csv(ab_file, sep='\t')

		# account for weird FLAIR tids
		df['base_tid'] = df.annot_transcript_id.str.extract(r'([A-z]+[0-9.]+)', expand=False)

		df = df.loc[(df.transcript_novelty == 'Known')&(df[datasets].any(axis=1))]
		num_known = len(df.base_tid.unique())

		# novel detections
		df = pd.read_csv(ab_file, sep='\t')
		df = df.loc[(df.transcript_novelty != 'Known')&(df[datasets].any(axis=1))]
		num_novel = len(df.index)
	
	# consider each dataset separately
	else: 
		num_known = []
		num_novel = []
		for d in datasets:

			# known detections
			df = pd.read_csv(ab_file, sep='\t')
			df['base_tid'] = df.annot_transcript_id.str.extract(r'([A-z]+[0-9.]+)', expand=False)
			df = df.loc[(df.transcript_novelty == 'Known')&(df[d] != 0)]
			num_known.append(len(df.base_tid.unique()))

			# novel detections
			df = pd.read_csv(ab_file, sep='\t')
			df = df.loc[(df.transcript_novelty == 'Novel')&(df[d] != 0)]
			num_novel.append(len(df.index))

	return num_known, num_novel

def main():
	args = get_args()
	sim_files = args.sim_files
	ab_file = args.ab_file
	tech_name = args.tech_name
	sim_name = args.sim_name
	merge_reps = args.merge_reps
	datasets = args.datasets

	# get the simulated data filenames and corresponding dataset names
	sim_files = sim_files.split(',')
	datasets = datasets.split(',')

	# count the total number of isoforms from the simulated data
	num_isoforms = count_total_isoforms(sim_files, merge_reps)

	# count the total number of isoforms detected in the 
	# analyzed data
	num_known, num_novel = count_detected_isoforms(ab_file, datasets, merge_reps)

	if not merge_reps: 
		for d, n, nk, nn in zip(datasets, num_isoforms, num_known, num_novel): 
			print()
			print('Simulated dataset: {}, {}, Technology: {}'.format(sim_name, d, tech_name))
			print('Detected {:,} out of {:,} known isoforms'.format(nk, n))
			print('Detected {:,} novel isoforms'.format(nn))
	else: 
		print('Simulated dataset: {}, Technology: {}'.format(sim_name, tech_name))
		print('Detected {:,} out of {:,} known isoforms'.format(num_known, num_isoforms))
		print('Detected {:,} novel isoforms'.format(num_novel))

if __name__ == '__main__': main()