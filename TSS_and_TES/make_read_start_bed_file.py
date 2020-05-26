import pandas as pd
from optparse import OptionParser
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import subprocess

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "TALON read annot file", metavar = "FILE", type = str)
    parser.add_option("--datasets", dest = "datasets", default = None,
        help = "Optional: comma-delim list of datasets to include", type = str)
    parser.add_option("--includeSpikes", dest = "spikes", action = "store_true",
                      default = False,
                      help = "Whether to include SIRV and ERCC reads. Default = False")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for outfile",
        metavar = "FILE", type = str)

    (options, args) = parser.parse_args()
    return options


def starts2bed(data, outprefix):
    """ Converts start positions from SAM read annot file to BED format.
        Returns name of outfile """

    bed_file = outprefix + "_known_read_starts.bed"

    bed_data = data[["chrom", "read_start"]].copy()
    bed_data["end"] = bed_data["read_start"]
    bed_data["read_start"] -= 1
    bed_data["name"] = data["read_name"]
    bed_data["score"] = "."
    bed_data["strand"] = data["strand"]

    bed_data.to_csv(bed_file, sep = '\t', header = False, index = False)
    return bed_file

def main():
    options = getOptions()

    data = pd.read_csv(options.infile, sep = '\t', header = 0)

    # Limit to known transcripts
    data = data.loc[data.transcript_novelty == "Known"]

    # Remove chrM reads
    data = data.loc[data.chrom != "chrM"]

    # Filter datasets (optional)
    if options.datasets != None:
        datasets = options.datasets.split(",")
        data = data[data['dataset'].isin(datasets)]

    # Remove spikes unless requested to keep
    if options.spikes == False:
         data = data[~data.chrom.str.contains("SIRV")]
         data = data[~data.chrom.str.contains("ERCC")]
 
    # Create BED file of read starts
    start_bed = starts2bed(data, options.outprefix) 
   
if __name__ == '__main__':
    main() 
