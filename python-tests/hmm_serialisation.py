#
# Copyright John Reid 2008
#

import _hmm, numpy, os

def assert_matrices_close(a1, a2):
    assert numpy.abs((a1-a2)).sum() < 1e-5

def assert_same(m1, m2):
    assert m1.M == m2.M
    assert m1.P == m2.P
    assert m1.N == m2.N
    assert_matrices_close(m1.A, m2.A)
    assert_matrices_close(m1.B, m2.B)

for markov_order in xrange(3):
    for binary in (True, False):
        filename = 'model.hmm'
        m = _hmm.build_fully_connected_hmm(3, 4 ** (markov_order+1), markov_order)
        m.serialise(filename, binary)
        m_copy = _hmm.Model.deserialise(filename, binary)
        os.remove(filename)
        assert_same(m, m_copy)
