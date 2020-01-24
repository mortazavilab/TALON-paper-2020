import sys
sys.path.append("..")
import profile_internal_priming as prof_ip

def test_frac_as():
    """ Compute the fraction of As in the sequence, making sure we don't have
        int rounding """

    assert prof_ip.compute_frac_As("AAAAAA") == 1
    assert prof_ip.compute_frac_As("AATG") == 0.5
    assert prof_ip.compute_frac_As("ACTGACTGG") == 2.0/9.0
