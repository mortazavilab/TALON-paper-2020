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

def plot_fraction_internal_primed_per_category(data, outprefix):
    """ Input: TALON read_annot file """

    cat_order = [ 'Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic',
                      'Genomic' ]
    novelty_colors = {'Known': '#009E73',
                      'ISM': '#0072B2',
                      'NIC': '#D55E00',
                      'NNC': '#E69F00',
                      'Genomic': '#F0E442',
                      'Antisense': '#000000',
                      'Intergenic': '#CC79A7'}

    # Manipulate data to get counts per category and percentages within novelty groups
    data['internal_primed'] = [ x  > 0.5 for x in list(data.fraction_As) ]
    primed_per_cat = data.groupby(["transcript_novelty", "internal_primed"]).size()
    primed_per_cat = primed_per_cat.reset_index()
    primed_per_cat.columns = ["transcript_novelty", "internal_primed", "count"]
    totals = primed_per_cat[["transcript_novelty", "count"]].groupby(["transcript_novelty"]).sum().reset_index()
    totals.columns = ["transcript_novelty", "total"]
    primed_per_cat = primed_per_cat.merge(totals, how = 'left', on = "transcript_novelty")
    primed_per_cat.columns = ["transcript_novelty", "internal_primed", "count", "total"]
    primed_per_cat['percent'] = primed_per_cat['count']*100.0/primed_per_cat['total']
 
    # Reorder
    primed_per_cat['order'] = primed_per_cat['transcript_novelty'].apply(lambda x: \
                                             cat_order.index(x))
    primed_per_cat = primed_per_cat.sort_values(by='order', ascending = True) 
    print(primed_per_cat)

    fname = outprefix + "_internalPrimed_reads_by_novelty_type.png"
    
    # Plot stacked bar chart
    n_cats = len(set(list(primed_per_cat.transcript_novelty)))
    x_range = np.arange(1, 10*n_cats, 10) # x locations
    width = 8       # the width of the bars
    plt.rcParams.update({'font.size': 12, 'font.family': 'sans-serif',
                         'hatch.linewidth': 5})
    fig, ax = plt.subplots()
 
    # Plot the bars one category at a time
    range_index = 0
    cats = []
    offset = max(primed_per_cat['total'])*0.01
    for novelty in cat_order:
        subdata = primed_per_cat.loc[primed_per_cat['transcript_novelty'] == novelty]
        if len(subdata) == 0:
            continue

        curr_range = x_range[range_index:range_index+1]
        range_index += 1

        ip = plt.bar(curr_range, 
                      subdata.loc[subdata.internal_primed == True]['count'], 
                      width = width, hatch='//', 
                      color = 'white', edgecolor=[novelty_colors[novelty]])


        other = plt.bar(curr_range, 
                        subdata.loc[subdata.internal_primed == False]['count'],
                        width = width, 
                        color = 'white', edgecolor=[novelty_colors[novelty]], 
                        bottom = subdata.loc[subdata.internal_primed == True]['count'])
    
        percent_ip = round(list(subdata.loc[subdata.internal_primed == True]['percent'])[0],1)
        plt.text(curr_range, 
                 subdata.loc[subdata.internal_primed == False]['total'] + offset, 
                 str(percent_ip) + "%", 
                 color = novelty_colors[novelty], ha="center")
        cats.append(ip)
      
    #plt.title('Number of reads with internal priming by novelty category')
    plt.xlabel('Category')
    plt.ylabel('Number of reads')
    plt.xticks(x_range, [])
    #plt.ylim([0,ymax])
    plt.tight_layout() 
    plt.savefig(fname, dpi = 600, bbox_inches='tight')

def main():
    options = get_options()

    data = pd.read_csv(options.read_annot, sep='\t', header = 0)

    if options.omit_genomic == True:
        data = data.loc[data.transcript_novelty != 'Genomic']

    # Plot fraction As by novelty category
    plot_fraction_internal_primed_per_category(data, options.outprefix)

if __name__ == '__main__':
    main()
