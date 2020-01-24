import sys
sys.path.append("..")
import profile_internal_priming as prof_ip

def test_fetch_range_plus():
    """
       Example:
                 1 2 3 4 5 6
                 | | | | | |
                 T A C G T C 
       desired length = 2
       If the transcript ends at pos 4 on the + strand, the correct range is 5,6
    """
    assert prof_ip.fetch_range_after_transcript(4, '+', 2) == (5,6)

def test_fetch_range_minus():
    """
       Example:
                 1 2 3 4 5 6
                 | | | | | |
                 T A C G T C
       desired length = 2
       If the transcript ends at pos 4 on the - strand, the correct range is 2,3
    """
    assert prof_ip.fetch_range_after_transcript(4, '-', 2) == (2,3)
