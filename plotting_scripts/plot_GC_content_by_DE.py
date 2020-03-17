import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
import gzip
from statannot import add_stat_annotation

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
                      help = "EdgeR results table")
    parser.add_option("--platform", dest = "platform",
                      help = "Name of long-read sequencing platform",
                      default = "PacBio")
    parser.add_option("--s", dest = "fasta",
                      help = "FASTA file of sequences fro known transcripts")
    parser.add_option("--ymax", dest = "ymax", type = float,
                      help = "Max value for y axis.")
    parser.add_option("--o", dest = "outprefix",
                      help = "Prefix for plot filename", metavar = "FILE",
                      type = "string")

    (options, args) = parser.parse_args()
    return options

def violin_plot(data, colname, ymax, fname):
    """ Plot a violin plot with the length of each read by novelty category"""

    sns.set_context("paper", font_scale=1.3)
    ax = sns.stripplot(x='DE_type', y=colname, data=data, color="black", 
                       alpha = 0.5, size = 1.5, jitter = True)

    ax = sns.boxplot(x='DE_type', y=colname, data=data, palette = "Blues")

    add_stat_annotation(ax, data=data, x='DE_type', y=colname,
                    box_pairs=[("Higher in Illumina", "Higher in PacBio"),
                               ("Higher in Illumina", "Not DE"),
                               ("Higher in PacBio", "Not DE")],
                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

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
        ax.text(pos[tick], ypos[tick] + ypos[tick]*0.1, nobs[tick],
        horizontalalignment='center', size='x-small', color='black', weight='semibold')

    ax.legend().set_visible(False)
    plt.xlabel("")
    plt.ylabel("GC percentage of gene")
    #ymin = min(data.groupby(['transcript_novelty'])['read_length'].min().values)
    plt.ylim(0, 100)
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

def compute_all_GCs(fasta):
    """ For each fasta transcript:
          1) Extract gene name
          2) Compute GC content of sequence
          3) Record gene name, transcript ID, and GC content in pandas df.
    """
    gene_names = []
    transcript_IDs = []
    GC_content = []

    try:
        with gzip.open(fasta, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                split_ID = (record.id).split("|")
                gene_name = split_ID[5]
                transcript_ID = split_ID[0]

                gene_names.append(gene_name)
                transcript_IDs.append(transcript_ID)
                GC_content.append(GC(record.seq))
               
    except Exception as e:
        print(e)
        sys.exit("Problem reading fasta sequence file. Expecting gzipped file")

    # Convert lists into a pandas data frame
    df = pd.DataFrame({"gene": gene_names, 
                       "transcript_ID": transcript_IDs, 
                       "GC": GC_content})
  
    return df

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

    # Read fasta file and compute GC content of entries 
    GC_table = compute_all_GCs(options.fasta)
    median_GC = GC_table.groupby("gene").median()

    # Merge data
    data = pd.merge(data, median_GC, on = "gene", how = "left")

    fname = options.outprefix + "DE_median_gene_GC-content.png"
    violin_plot(data, "GC", options.ymax, fname)


if __name__ == '__main__':
    main()
