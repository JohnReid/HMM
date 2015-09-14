#
# Copyright John Reid 2008
#

"""
Code to test count k-mer functionality.
"""

if '__main__' == __name__:
    import hmm
    from gapped_pssms.sequence import seq_to_numpy


    sequences = [
      'aacc',
      'ggtt',
      'aaaaaaaaaaaa',
    ]
    seqs = map(seq_to_numpy, sequences)
    print 'Sequences:'
    print '\n'.join(str(s) for s in seqs)
    print

    top_k_mers = hmm.top_mers_by_sequence_membership(
      seqs,
      k=3,
      n=10
    )
    print '\n'.join(str(o) for o in top_k_mers)
