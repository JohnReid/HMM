#
# Copyright John Reid 2007, 2008
#

"""
Code for reverse complements, conversion from strings to numpy arrays and vice versa.
"""


import numpy.random

def base_to_index(base):
    if 'a' == base or 'A' == base: return 0
    if 'c' == base or 'C' == base: return 1
    if 'g' == base or 'G' == base: return 2
    if 't' == base or 'T' == base: return 3
    return 4

def index_to_base(base):
    if 0 == base: return 'a'
    if 1 == base: return 'c'
    if 2 == base: return 'g'
    if 3 == base: return 't'
    if 4 == base: return 'n'
    raise RuntimeError('Unknown base: %d' % base)

def seq_to_numpy(seq):
    result = numpy.empty(len(seq), dtype=int)
    for i, s in enumerate( seq ):
        result[ i ] = base_to_index( s )
    return result

def seqs_to_numpy(seqs):
    return [ seq_to_numpy(seq) for seq in seqs ]

def numpy_to_seq(seq):
    return ''.join(index_to_base(b) for b in seq)

def complement(b):
    return (4 == b) and 4 or (3 - b)

def rev_comp(seq):
    result = numpy.empty_like(seq)
    for i, b in enumerate(seq):
        result[-i-1] = complement(b)
    return result

def num_bases(seqs):
    return sum(
            [
                    len( seq )
                    for seq in seqs
            ]
    )

def num_known_bases(seqs):
    return sum(
            [
                    sum( seq != 4 )
                    for seq in seqs
            ]
    )

def random_sequence(length):
    return numpy.random.randint( 4, size = length )

def random_sequences(num_seqs, average_length):
    return [
            random_sequence( length = numpy.random.poisson( lam = average_length ) )
            for i in xrange( num_seqs )
    ]
