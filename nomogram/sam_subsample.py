import random
import sys
from optparse import OptionParser
def getOptions():
    parser = OptionParser()
    parser.add_option("--sam", dest = "sam",
                      help = "Input SAM file", metavar = "FILE", type = str)
    parser.add_option("--n", dest = "n",
                      help = "Number of reads to sample per cell",type = int)
    parser.add_option("--o", dest = "outfile",
                      help = "Output filename", metavar = "FILE", type = str)
    (options, args) = parser.parse_args()
    return options
def subsample_reads(sam, N, outfile, seed = 30):
    random.seed(seed)
    header = []
    lines = []
    n_lines = 0
    with open(sam, 'r') as f:
        for l in f:
            l = l.strip()
            if l.startswith("@"):
                header.append(l)
            else:
                lines.append(l)
                n_lines += 1
    if N > n_lines:
        print("File does not have enough lines.")
        sys.exit(1)
    with open(outfile, 'w') as o:
        for line in header:
            o.write(line + "\n")
        for entry in random.sample(lines, N):
            o.write(entry + "\n")
def main():
    options = getOptions()
    subsample_reads(options.sam, options.n, options.outfile)
if __name__ == '__main__':
    main()