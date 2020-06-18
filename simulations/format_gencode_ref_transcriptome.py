import argparse

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-rt', dest='rt')
	parser.add_argument('-ofile', dest='ofile')
	args = parser.parse_args()
	return args

def main():

	args = get_args()
	with open(args.rt, 'r') as infile:
		with open(args.ofile, 'w') as outfile:
			for line in infile:
				if line.startswith('>'):
					header = line.split('|')[0]
					header = header.split('_')[0]+'\n'
					outfile.write(header)
				else:
					outfile.write(line)

if __name__ == '__main__': main()