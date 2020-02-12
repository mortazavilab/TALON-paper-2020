import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from optparse import OptionParser

def get_options():
    parser = OptionParser(description = ("Plot"))

    parser.add_option("--f", dest = "read_annot",
                      help = ("TALON read_annot file"))
    parser.add_option("--omitGenomic", dest ="omit_genomic", action='store_true',
                  help = ("If this option is set, reads from the Genomic "
                          "novelty category will be omitted in the output."),
                   default = False)
    parser.add_option("--outprefix", dest = "outprefix",
                      help = ("Prefix for outfile"))

    (options, args) = parser.parse_args()
    return options


def main():
    options = get_options()

    data = pd.read_csv(options.read_annot, sep='\t', header = 0)
    data['percent_As'] = data.fraction_As*100

    if options.omit_genomic == True:
        data = data.loc[data.transcript_novelty != 'Genomic']
        cat_order = [ 'Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic']
    else:
        cat_order = [ 'Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic', 
                      'Genomic' ]

    # Plot fraction As by novelty category
    fname = options.outprefix + "_percentA_distribution_by_novelty_type.png"

    novelty_colors = {'Known': '#009E73', 
                      'ISM': '#0072B2',
                      'NIC': '#D55E00', 
                      'NNC': '#E69F00',
                      'Genomic': '#F0E442',
                      'Antisense': '#000000',
                      'Intergenic': '#CC79A7'}

    g = sns.FacetGrid(data, col="transcript_novelty", hue = "transcript_novelty",
                      palette = novelty_colors, hue_order = cat_order)
    g = g.map(sns.distplot, "percent_As")
    g.set_ylabels('Density')
    g.set_xlabels('Percent As after read end (10 bp)')

    # Try to put an n label on the plots
    #for ax, title in zip(g.axes.flat, col_order):

    plt.savefig(fname, dpi = 600, bbox_inches='tight')


if __name__ == '__main__':
    main()
