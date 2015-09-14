#
# Copyright John Reid 2007, 2008
#

"""
Code to model PSSMs (gapped and ungapped) using HMMs.
"""



# from background import * # only need this stuff for old style background models
from find_pssms import *
from gapped_pssm_traits import *
from model_builder import *
# from pssm import * # only need this stuff for old style background models
from pssm_traits import *
from sequence import *


def base_information_content(base_dist, log_bg = None, M = 4):
    '''
    The information content of a distribution
    '''
    from math import log
    import numpy

    if None == log_bg:
        log_bg = numpy.log(numpy.ones(M) / M)

    return sum(
            (p > 0 and p*(log(p) - log_b or 0.0))
            for p, log_b in zip(base_dist, log_bg)
    )



def information_content(pssm, bg=None, M=4):
    '''
    The information content of a pssm

    See U{http://rsat.ulb.ac.be/rsat/tutorials/tut_PSSM.html}
    '''
    from math import log
    import numpy

    if None == bg:
        # assume uniform background
        bg = numpy.ones(M) / M
    log_bg = numpy.log(bg)

    return sum(
            base_information_content(base, log_bg, M)
            for base in pssm
    )



def _safe_x_log_x(x):
    """
    @return: x * log(x), with zeros where x is zero
    """
    result = numpy.log(x) * x
    result[0==x] = 0.
    return result



def entropy(emissions, probs):
    "@return: The entropy of the emissions at all positions weighted by the probs."
    return -(probs * (_safe_x_log_x(emissions).sum(axis=1))).sum()
