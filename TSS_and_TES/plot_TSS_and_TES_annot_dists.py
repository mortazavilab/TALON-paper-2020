import pandas as pd
from optparse import OptionParser
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "TALON read annot file", metavar = "FILE", type = str)
    parser.add_option("--ref", dest = "ref_sites",
        help = "TSV file of GENCODE starts and ends for each transcript", 
        metavar = "FILE", type = str)
    parser.add_option("--datasets", dest = "datasets", default = None,
        help = "Optional: comma-delim list of datasets to include", type = str)
    parser.add_option("--includeSpikes", dest = "spikes", action = "store_true",
                      default = False,
                      help = "Whether to include SIRV and ERCC reads. Default = False")
    parser.add_option("--xmax", dest = "xmax", type = int,
                      help = "Max x value for plots.", default = 5000)
    parser.add_option("--ymax", dest = "ymax", type = int,
                      help = "Max y value for plots.", default = 1000000)
    parser.add_option("--o", dest = "outprefix", help = "Prefix for outfile",
        metavar = "FILE", type = str)

    (options, args) = parser.parse_args()
    return options

def plot_histogram(data, xvar, label, xmax, ymax, fname):
    x = pd.Series(data[xvar], name=label)
    plt.xlim(-1*xmax,xmax)
    plt.ylim(0,ymax)
    ax = sns.distplot(x, kde = False, color = "dodgerblue",
                  bins = np.arange(min(x), max(x), 5))

    med = round(np.median(x), 1)

    style = dict(size=12, color='black')
    plt.axvline(med, linestyle = '--', color = 'lightgrey')

    ax.text(med + 50, ymax*7/8, "Median: " + str(med) + " bp", **style)

    plt.tight_layout()
    plt.savefig(fname, dpi = 600, bbox_inches='tight')
    plt.close()

def main():
    options = getOptions()

    data = pd.read_csv(options.infile, sep = '\t', header = 0)
    ref_sites = pd.read_csv(options.ref_sites, sep = '\t', header = 0)

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

    # Merge together reads with GENCODE start data
    data = data[["annot_transcript_id", "strand", "read_start", "read_end"]]
    data = pd.merge(data, ref_sites, on = "annot_transcript_id", 
                    how = "left")

    # Compute TSS/TES distances
    # Negative distance is upstream, positive is downstream
    data['TSS_dist'] = data.read_start - data.TSS_pos
    data['TES_dist'] = data.read_end - data.TES_pos
    data.loc[data['strand']=='-', ['TSS_dist', 'TES_dist']] *= -1

    # Plot
    plot_histogram(data, "TSS_dist", "Distance from TSS annotated to intron chain (bp)", 
                   options.xmax, options.ymax, 
                   options.outprefix + "_TSS_dist_known.png")
    plot_histogram(data, "TES_dist", "Distance from TES annotated to intron chain (bp)",
                   options.xmax, options.ymax,
                   options.outprefix + "_TES_dist_known.png")

if __name__ == '__main__':
    main()
