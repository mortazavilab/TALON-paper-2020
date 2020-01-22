# Takes an exon-only GTF file and identifies cases where a gene has transcripts on
# more than one strand. Split those transcripts into different genes

from optparse import OptionParser
import pandas as pd
import csv

def get_options():
    parser = OptionParser(description=("Takes an exon-only GTF file and "
                                       "identifies cases where a gene has "
                                       "transcripts on more than one strand. "
                                       "Split those transcripts into different "
                                       "genes."))
    parser.add_option("--f", dest = "gtf_file",
                      help = "GTF file containing exon entries only")
    parser.add_option("--outfile", dest = "outfile",
                      help = "Outfile name")

    (options, args) = parser.parse_args()
    return options

def main():

    options = get_options()

    # Read GTF
    gtf = pd.read_csv(options.gtf_file, sep = "\t", header=None)
    gtf.columns = ["chrom", "source", "feat", "start", "end", "misc1",
                   "strand", "misc2", "description"]

    # Create gene_id, transcript_id columns
    gene_ids = [] 
    transcript_ids = []   
    for values in gtf['description']:
        gene_ids.append(values.split('gene_id "')[-1].split('"')[0]) 
        transcript_ids.append(values.split('transcript_id "')[-1].split('"')[0])

    gtf['gene_id'] = gene_ids
    gtf['transcript_id'] = transcript_ids

    # Sub-number genes by gene ID 
    gene_tab = gtf[['gene_id', 'strand']].drop_duplicates()
    gene_tab["subgene"] = gene_tab.groupby(['gene_id']).cumcount()+1

    # Now merge this back into the GTF data frame
    gtf = pd.merge(gtf, gene_tab, how = 'left', on = ['gene_id', 'strand'])
    gtf['gene_id'] = gtf.gene_id + ['.']*len(gtf) + gtf["subgene"].map(str)
    gtf['transcript_id'] = gtf.gene_id + ['-']*len(gtf) + gtf["transcript_id"]

    # Write to file
    description = []
    for index, row in gtf.iterrows():
        description.append('gene_id "%s"; transcript_id "%s";' % \
                           (row['gene_id'], row['transcript_id']))
    gtf.description = description
    gtf = gtf[ ["chrom", "source", "feat", "start", "end", "misc1",
                "strand", "misc2", "description"] ]
    
    gtf.to_csv(options.outfile, sep = '\t', header = False, index = False,
               quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
