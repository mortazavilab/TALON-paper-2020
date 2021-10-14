import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
import numpy as np

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
                      help = "TALON read annot file")
    parser.add_option("--datasets", dest = "datasets", default = None,
                      help = ("Comma-delimited list of datasets to use. "
                               "Default behavior is to use all datasets."))
    parser.add_option("--whitelist", dest = "whitelist", default = None,
                      help = "CSV TALON whitelist to filter transcripts")
    parser.add_option("--platform", dest = "platform",
                      help = "Sequencing platform",
                      default = "PacBio")
    parser.add_option("--ymax", dest = "ymax",
                      help = "Max value for y axis.",
                      default = 10)
    parser.add_option("--o", dest = "outprefix",
                      help = "Prefix for plot filename", metavar = "FILE",
                      type = "string")

    (options, args) = parser.parse_args()
    return options


def violin_plot(data, platform, fname):
    """ Plot a violin plot with the length of each read by novelty category"""

    sns.set_context("paper", font_scale=1.3)
    #ax = sns.stripplot(x='transcript_novelty', y='read_length', data=data, color="grey", jitter = True)

    color_dict = {'Known': '#009E73', 
                      'ISM': '#0072B2',
                      'NIC': '#D55E00', 
                      'NNC': '#E69F00',
                      'Genomic': '#F0E442',
                      'Antisense': 'gray',
                      'Intergenic': '#CC79A7'}
    cat_order = [ "Known", "ISM", "NIC", "NNC", "Antisense", "Intergenic"]#, "Genomic"]

    data.transcript_novelty = data.transcript_novelty.astype("category")
    data.transcript_novelty.cat.set_categories(cat_order, inplace=True)
    data = data.sort_values(["transcript_novelty"])

    ax = sns.violinplot(x='transcript_novelty', y="read_length", legend = False,
                        data=data, palette=color_dict,
                        order=cat_order,
                        linewidth = 1, saturation = 0.8, color = 'black',
                        inner = 'box', cut = 0)

    # Calculate number of obs per group & position labels
    nobs = data['transcript_novelty'].value_counts().values
    nobs = [str(x) for x in nobs.tolist() if x != 0]
    nobs = ["n=" + i for i in nobs]

    print(data.groupby(['transcript_novelty'])['read_length'].min()) 
    # Add it to the plot
    ypos = data.groupby(['transcript_novelty'])['read_length'].max().dropna().values
    print(data.groupby(['transcript_novelty'])['read_length'].max().dropna())
    pos = range(len(nobs))
    for tick,label in zip(pos,ax.get_xticklabels()):
        ax.text(pos[tick], ypos[tick] + 200, nobs[tick],
        horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.legend().set_visible(False)
    plt.xlabel("")
    plt.ylabel(platform + " read length")
    #ymin = min(data.groupby(['transcript_novelty'])['read_length'].min().values)
    plt.ylim(0, max(ypos)*1.1)
    plt.tight_layout()
    plt.savefig(fname, dpi = 600, bbox_inches='tight')
    plt.close()

def filter_by_datasets(data, datasets):
    """ Given a list of datasets, limit the data frame to reads coming from
        only thoes datasets """

    # Validate dataset names first
    permitted_datasets = set(list(data.dataset))
    for d in datasets:
        if d not in permitted_datasets:
            raise ValueError("Invalid dataset name: %s" % (d))

    data = data[data['dataset'].isin(datasets)]
    return data

def filter_by_whitelist(data, whitelist_file):
    """ Limit reads to those that are assigned to a whitelist transcript"""

    whitelist = pd.read_csv(whitelist_file, sep=",", header = None)
    whitelist.columns = ["gene_ID", "transcript_ID"]

    data = pd.merge(data, whitelist, on = "transcript_ID", how = "right").dropna()

    return data

def main():
    options = getOptions()
   
    data = pd.read_csv(options.infile, sep = "\t", header = 0)

    # Limit to datasets of interest
    if options.datasets != None:
        datasets = options.datasets.split(",")
        data = filter_by_datasets(data, datasets)
        

    # Remove SIRV/ERCCs
    filter = data['chrom'].str.contains("ERCC")
    data = data[~filter]
    filter = data['chrom'].str.contains("SIRV")
    data = data[~filter]

    # Now apply the whitelist
    if options.whitelist != None:
        data = filter_by_whitelist(data, options.whitelist)

    # Sort the reads by novelty category?

    # Now plot the read length by novelty category
    if datasets:
        fname = options.outprefix + "_".join(datasets) + options.platform + \
                "_read_len_by_novelty.png"
    else:
        fname = options.outprefix + options.platform + \
                "_read_len_by_novelty.png"
    
    violin_plot(data, options.platform, fname)

if __name__ == '__main__':
    main()
