import argparse

def parse_args():
	# argument parse config stuff
	parser = argparse.ArgumentParser(description=\
	    'Separates a gtf file into novelty classifications and'+\
	    ' generates a text file with track information to paste into genome browser')
	parser.add_argument('-gtf', help='gtf file')