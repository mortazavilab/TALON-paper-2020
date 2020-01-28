import pandas as pd
import numpy as np
from optparse import OptionParser

def get_options():
    parser = OptionParser(description = ("Plot"))

    parser.add_option("--f", dest = "read_annot",
                      help = ("TALON read_annot file"))
    parser.add_option("--a", dest = "frac_a",
                      help = ("File of read IDs and fraction As in sequence "
                              "following the transcript end. To provide more "
                              "than one file, comma-delimit them."))
    parser.add_option("--outprefix", dest = "outprefix",
                      help = ("Prefix for outfile"))

    (options, args) = parser.parse_args()
    return options


def main():
    options = get_options()
    frac_a_files = options.frac_a.split(",")

    # Concatenate fraction A dfs
    frames = [ pd.read_csv(f, sep='\t', header = 0) for f in frac_a_files ]
    frac_a = pd.concat(frames)

    read_annot = pd.read_csv(options.read_annot, sep='\t', header = 0)

    # Merge the data.
    data = pd.merge(read_annot, frac_a, how = 'left', 
                    on = "read_name").dropna()[["gene_ID", "transcript_ID", 
                                                "transcript_novelty", "fraction_As"]]
    data = data.loc[data.transcript_novelty != "Genomic"]

    # Keep rows where transcript is known OR where fraction As is <= 0.5
    # Use the max fracA value for transcripts seen more than once
    data = data.groupby(["gene_ID", "transcript_ID","transcript_novelty"]).max()
    data = data.reset_index()
    data = data.loc[(data.transcript_novelty == 'Known') | (data.fraction_As <= 0.5)]
    whitelist = data[["gene_ID", "transcript_ID", "transcript_novelty"]].drop_duplicates()
    print(whitelist.groupby(['transcript_novelty']).count())

    # Make a TALON filter file
    fname = options.outprefix + "_frac_A_whitelist.csv"
    whitelist.to_csv(fname, header = False, sep=",", index = False)

if __name__ == '__main__':
    main()
