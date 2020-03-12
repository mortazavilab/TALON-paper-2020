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
    parser.add_option('--maxFracA', dest='a_frac', default=0.5,
                  help='The maximum fraction of As to be considered not internally primed. Default=0.5')
    parser.add_option("--outprefix", dest = "outprefix",
                      help = ("Prefix for outfile"))

    (options, args) = parser.parse_args()
    return options

def plot_fraction_internal_primed_per_category(data, outprefix, a_frac):
    """ Input: TALON read_annot file """

    cat_order = [ 'Known', 'ISM', 'Prefix ISM', 'Suffix ISM', 
                  'NIC', 'NNC', 'Antisense', 'Intergenic',
                      'Genomic' ]
    novelty_colors = {'Known': '#009E73',
                      'ISM': '#0072B2',
                      'Prefix ISM': '#56B4E9',
                      'Suffix ISM': '#698BAC',
                      'NIC': '#D55E00',
                      'NNC': '#E69F00',
                      'Genomic': '#F0E442',
                      'Antisense': '#000000',
                      'Intergenic': '#CC79A7'}

    # Manipulate data to get counts per category and percentages within novelty groups
    data['internal_primed'] = [ float(x)  > float(a_frac) for x in list(data.fraction_As) ]
    primed_per_cat = data.groupby(["transcript_novelty", "ISM_subtype", "internal_primed"]).size()
    primed_per_cat = primed_per_cat.reset_index()
    primed_per_cat.columns = ["transcript_novelty", "ISM_subtype", "internal_primed", "count"]
    totals = primed_per_cat[["transcript_novelty", "ISM_subtype", "count"]].groupby(["transcript_novelty", "ISM_subtype"]).sum().reset_index()


    # add "both" ISM counts to both the prefix and suffix categories for totals
    both_count = totals.loc[(totals.transcript_novelty == 'ISM') & (totals.ISM_subtype == 'Both'), 'count'].tolist()[0]
    pre_count = totals.loc[(totals.transcript_novelty == 'ISM') & (totals.ISM_subtype == 'Prefix'), 'count'].tolist()[0]
    suf_count = totals.loc[(totals.transcript_novelty == 'ISM') & (totals.ISM_subtype == 'Suffix'), 'count'].tolist()[0]
    totals.loc[(totals.transcript_novelty == 'ISM') & (totals.ISM_subtype == 'Prefix'), 'count'] = pre_count+both_count
    totals.loc[(totals.transcript_novelty == 'ISM') & (totals.ISM_subtype == 'Suffix'), 'count'] = suf_count+both_count
    totals = totals[totals.ISM_subtype != 'Both']

    # add "both" ISM counts to both the prefix and suffix categories for primed_per_cat
    primed_per_cat.set_index(['transcript_novelty', 'ISM_subtype', 'internal_primed'], inplace=True)

    def get_count(df, ind1, ind2, ind3):
      try: 
        count = df.loc[ind1, ind2, ind3]['count']
      except:
        count = 0
      return count

    both_true_count = get_count(primed_per_cat, 'ISM', 'Both', True)
    both_false_count = get_count(primed_per_cat, 'ISM', 'Both', False)
    pre_true_count = get_count(primed_per_cat, 'ISM', 'Prefix', True)
    pre_false_count = get_count(primed_per_cat, 'ISM', 'Prefix', False)
    suf_true_count = get_count(primed_per_cat, 'ISM', 'Suffix', True)
    suf_false_count = get_count(primed_per_cat, 'ISM', 'Suffix', False)


    primed_per_cat.loc['ISM', 'Prefix', True]['count'] = pre_true_count+both_true_count
    primed_per_cat.loc['ISM', 'Prefix', False]['count'] = pre_false_count+both_false_count
    primed_per_cat.loc['ISM', 'Suffix', True]['count'] = suf_true_count+both_true_count
    primed_per_cat.loc['ISM', 'Suffix', False]['count'] = suf_false_count+both_false_count

    primed_per_cat.reset_index(inplace=True)
    primed_per_cat = primed_per_cat[primed_per_cat.ISM_subtype != 'Both']


    # rename ISM columns to make more sense
    prefix_name = 'Prefix ISM'
    suffix_name = 'Suffix ISM'

    primed_per_cat.loc[(primed_per_cat.transcript_novelty=='ISM')&(primed_per_cat.ISM_subtype=='Prefix'),
      'transcript_novelty']=prefix_name
    primed_per_cat.loc[(primed_per_cat.transcript_novelty=='ISM')&(primed_per_cat.ISM_subtype=='Suffix'),
      'transcript_novelty']=suffix_name
    primed_per_cat.drop('ISM_subtype', axis=1, inplace=True)

    totals.loc[(totals.transcript_novelty=='ISM')&(totals.ISM_subtype=='Prefix'), 'transcript_novelty']=prefix_name
    totals.loc[(totals.transcript_novelty=='ISM')&(totals.ISM_subtype=='Suffix'), 'transcript_novelty']=suffix_name
    totals.drop('ISM_subtype', axis=1, inplace=True)

    print(primed_per_cat)

    totals.columns = ["transcript_novelty", "total"]
    primed_per_cat = primed_per_cat.merge(totals, how = 'left', on = "transcript_novelty")
    primed_per_cat.columns = ["transcript_novelty", "internal_primed", "count", "total"]
    primed_per_cat['percent'] = primed_per_cat['count']*100.0/primed_per_cat['total']
 
    # Reorder
    primed_per_cat['order'] = primed_per_cat['transcript_novelty'].apply(lambda x: \
                                             cat_order.index(x))
    primed_per_cat = primed_per_cat.sort_values(by='order', ascending = True) 

    fname = outprefix + "_internalPrimed_reads_by_novelty_type.png"
    
    # Plot stacked bar chart
    n_cats = len(set(list(primed_per_cat.transcript_novelty)))
    x_range = np.arange(1, 10*n_cats, 10) # x locations
    width = 6.2       # the width of the bars
    plt.rcParams.update({'font.size': 15, 'font.family': 'sans-serif',
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
                      subdata.loc[subdata.internal_primed == True, 'count'], 
                      width = width, hatch='//', 
                      color = 'white', edgecolor=[novelty_colors[novelty]])


        other = plt.bar(curr_range, 
                        subdata.loc[subdata.internal_primed == False, 'count'],
                        width = width, 
                        color = 'white', edgecolor=[novelty_colors[novelty]], 
                        bottom = subdata.loc[subdata.internal_primed == True, 'count'])
         
        try:
          percent_ip = round(list(subdata.loc[subdata.internal_primed == True, 'percent'])[0],1)
        except: percent_ip = 0.0
        plt.text(curr_range, 
                 subdata.loc[subdata.internal_primed == False, 'total'] + offset, 
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
    plot_fraction_internal_primed_per_category(data, options.outprefix, options.a_frac)

if __name__ == '__main__':
    main()
