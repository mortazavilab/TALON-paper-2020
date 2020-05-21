from optparse import OptionParser

def getOptions():
    parser = OptionParser()

    parser.add_option("--gtf", dest = "gtf",
        help = "GENCODE GTF annotation", metavar = "FILE", type = str)
    parser.add_option("--o", dest = "outfile", help = "Name for outfile",
        metavar = "FILE", type = str)

    (options, args) = parser.parse_args()
    return options

def parse_out_transcript_ID(gtf_metadata):
    """ Parse GTF metadata in order to extract the transcript ID """
    transcript_ID = (gtf_metadata.split('transcript_id "')[1]).split('";')[0]
    return(transcript_ID)

def get_annot_starts_and_ends(gtf, outfile):
    """ Given a GTF file, creates a tab-delimited outfile with the name, TSS,
        and TES of each transcript. Skips chrM."""

    outfile = open(outfile, 'w')
    outfile.write("\t".join(["annot_transcript_id", "TSS_pos", "TES_pos"]) + "\n")

    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            line = line.split("\t")

            # Skip things that are not transcripts
            if line[2] != "transcript":
                continue

            # Skip circular chromosome
            if line[0] == "chrM":
                continue

            gtf_start = int(line[3])
            gtf_end = int(line[4])
            strand = line[6]
            meta = line[-1]

            transcript_ID = parse_out_transcript_ID(meta)

            if strand == "+":
                TSS = gtf_start
                TES = gtf_end
            elif strand == "-":
                TSS = gtf_end
                TES = gtf_start
            else:
                raise ValueError("Invalid strand")

            entry = [transcript_ID, TSS, TES]
            outfile.write("\t".join([str(x) for x in entry]) + "\n")

    outfile.close()

def main():
    options = getOptions()

    get_annot_starts_and_ends(options.gtf, options.outfile)

if __name__ == '__main__':
    main()
