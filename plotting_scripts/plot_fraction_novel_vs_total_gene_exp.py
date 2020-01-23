import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from optparse import OptionParser

def get_options():
    parser = OptionParser(description = ("Plot the fraction novel reads against "
                                         "the total read count in a dataset"))

    parser.add_option("--f", dest = "infile",
                      help = ("TALON unfiltered abundance file"))
    parser.add_option("--dataset", dest = "dataset",
                      help = ("Name of dataset to use"))
    parser.add_option("--outprefix", dest = "outprefix",
                      help = ("Prefix for outfile. Dataset name will be added "
                              "automatically."))

    (options, args) = parser.parse_args()
    return options

def main():
    options = get_options()
    dataset = options.dataset
    data = pd.read_csv(options.infile, sep="\t", header = 0)

    # For each gene, compute the total reads annotated to it as well as the
    # fraction of reads mapped to novel models
    gene_info = data[["annot_gene_id", dataset]].groupby("annot_gene_id").sum()
    gene_info.rename(columns={dataset: "total_reads"}, inplace=True)

    novel_per_gene = data.loc[data.transcript_novelty != "Known"].groupby("annot_gene_id").sum()
    novel_per_gene.rename(columns={dataset: "novel_reads"}, inplace=True)
    gene_info = gene_info.merge(novel_per_gene['novel_reads'], 
                                on = "annot_gene_id", how = 'left').fillna(0)
    gene_info['fraction_novel'] = 100.0*gene_info["novel_reads"]/gene_info["total_reads"]
    print(gene_info) 

    # Plot
    fname = options.outprefix + "_" + dataset + "_gene-expr_vs_frac-novel-reads.png"
    sns.scatterplot(gene_info['total_reads'], gene_info['fraction_novel'], markers = {'s':6})
    plt.xlabel("Total read count for gene")
    plt.ylabel("Fraction of reads annotated as novel")


    plt.savefig(fname, dpi = 600, bbox_inches='tight')

if __name__ == '__main__':
    main()
