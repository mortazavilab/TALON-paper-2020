import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
import numpy as np

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
                      help = "EdgeR results table")
    parser.add_option("--platform", dest = "platform",
                      help = "Name of long-read sequencing platform",
                      default = "PacBio")
    parser.add_option("--col", dest = "colname",
                      help = "Name of column in file to use for length",
                      default = "median_length")
    parser.add_option("--ymax", dest = "ymax",
                      help = "Max value for y axis.",
                      default = 10)
    parser.add_option("--o", dest = "outprefix",
                      help = "Prefix for plot filename", metavar = "FILE",
                      type = "string")

    (options, args) = parser.parse_args()
    return options

def violin_plot(data, colname, ymax, fname):
    """ Plot a violin plot with the length of each read by novelty category"""

    sns.set_context("paper", font_scale=1.3)
    #ax = sns.stripplot(x='transcript_novelty', y='read_length', data=data, color="grey", jitter = True)

    #data.transcript_novelty = data.transcript_novelty.astype("category")
    #data.transcript_novelty.cat.set_categories(cat_order, inplace=True)
    #data = data.sort_values(["transcript_novelty"])
    ax = sns.boxplot(x='DE_type', y=colname, data=data)

    #ax = sns.violinplot(x='DE_type', y=colname, legend = False,
    #                    data=data,
    #                    #order=cat_order,
    #                    linewidth = 1,
    #                    inner = 'box', cut = 0)

    # Calculate number of obs per group & position labels
    nobs = list(data.groupby("DE_type").size())
    nobs = [str(x) for x in nobs]
    nobs = ["n=" + i for i in nobs]

    # Add it to the plot
    ypos = data.groupby(['DE_type'])[colname].max().dropna().values
    pos = range(len(nobs))
    for tick,label in zip(pos,ax.get_xticklabels()):
        ax.text(pos[tick], ypos[tick] + 500, nobs[tick],
        horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.legend().set_visible(False)
    ax.set_yscale("log", basey=10)
    plt.xlabel("")
    plt.ylabel("Gene length (%s of transcripts)" % colname)
    #ymin = min(data.groupby(['transcript_novelty'])['read_length'].min().values)
    #plt.ylim(0, max(ypos)*1.1)
    plt.tight_layout()
    plt.savefig(fname, dpi = 600, bbox_inches='tight')
    plt.close()

def label_DE(row):
    if row["status"] == "not_sig":
        return "Not DE"
    elif row["status"] == "significant":
        if row["logFC"] > 0:
            return "higher"
        elif row["logFC"] < 0:
            return "lower"
        else:
            raise ValueError("Expected logFC to be nonzero")

    else:
        raise ValueError("Unexpected value in 'status' column")

def main():
    options = getOptions()

    data = pd.read_csv(options.infile, sep = "\t", header = 0)

    # Categorize data as higher in long-read platform, higher in Illumina, or
    # neither 
    higher = "Higher in %s" % (options.platform)
    lower = "Higher in Illumina"

    data['DE_type'] = data.apply(label_DE, axis=1)
    data = data.replace("higher", higher)
    data = data.replace("lower", lower)
    print(data.groupby("DE_type").size())

    fname = options.outprefix + "DE_gene_length_by_" + options.colname + ".png"
    violin_plot(data, options.colname, options.ymax, fname)

    print(data.loc[data[options.colname].idxmax()])

if __name__ == '__main__':
    main()
