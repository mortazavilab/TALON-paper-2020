import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from optparse import OptionParser

def get_options():
    parser = OptionParser(description = ("Plot"))

    parser.add_option("--f", dest = "read_annot",
                      help = ("TALON read_annot file"))
    parser.add_option("--a", dest = "frac_a",
                      help = ("File of read IDs and fraction As in sequence "
                              "following the transcript end."))
    parser.add_option("--outprefix", dest = "outprefix",
                      help = ("Prefix for outfile"))

    (options, args) = parser.parse_args()
    return options

def main():
    options = get_options()

    read_annot = pd.read_csv(options.read_annot, sep='\t', header = 0)
    frac_a = pd.read_csv(options.frac_a, sep='\t', header = 0)

    # Merge the data
    data = pd.merge(read_annot, frac_a, how = 'right', on = "read_name")
    data = data.loc[data.transcript_novelty != "Genomic"]

    # Plot fraction As by novelty category
    fname = options.outprefix + "_fracA_distribution_by_novelty_type.png"

    novelty_colors = {'Known': '#009E73', 
                      'ISM': '#0072B2',
                      'NIC': '#D55E00', 
                      'NNC': '#E69F00',
                      'Genomic': '#F0E442',
                      'Antisense': '#000000',
                      'Intergenic': '#CC79A7'}
    cat_order = [ 'Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic' ]

    g = sns.FacetGrid(data, col="transcript_novelty", hue = "transcript_novelty",
                      palette = novelty_colors, hue_order = cat_order)
    g = g.map(sns.distplot, "fraction_As").set_ylabels('Density')

    plt.savefig(fname, dpi = 600, bbox_inches='tight')

if __name__ == '__main__':
    main()
