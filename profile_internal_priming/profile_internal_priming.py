import pyfaidx

def fetch_seq(chrom=str, start=int, stop=int, strand=str, genome=pyfaidx.Fasta, 
              indexing = 0):
    """ Given a genomic interval, return the sequence with respect to the 
        strand supplied. 
        If 1-based indexing is specified, then 1 will be subtracted from the
        position to convert to the Python indexing. """

    if start > stop:
        raise ValueError("Start must be less than or equal to stop")

    if indexing != 0:
        if indexing == 1:
            start -= 1
        else:
            raise ValueError("Valid indexing modes include: 1 or 0")

    seq = genome[chrom][start:stop]

    if strand == "-":
        seq = seq.reverse.complement

    return str(seq)







#    n = str(seq).count('A')
#    return n
