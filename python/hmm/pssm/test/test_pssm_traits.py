#
# Copyright John Reid 2007, 2008
#

"""
Code to test PSSM traits.
"""


import hmm, hmm.pssm, unittest

class TestPssmTraits(unittest.TestCase):
    def test_traits(self):
        from hmm.pssm import create_background_model, PssmTraits, seq_to_numpy
        from infpy import check_is_close_2

        p_binding_site = .01
        num_background_states = 2
        emission_dists = [
          [ 1., 0., 0., 0. ],
          [ 0., 1., 0., 0. ],
          [ 0., 1., 0., 0. ],
          [ 1., 0., 0., 0. ],
          [ 0., 0., 1., 0. ],
          [ 0., 0., 0., 1. ],
          [ 0., 0., 0., 1. ],
          [ 0., 0., 0., 1. ],
          [ 0., 0., 1., 0. ],
          [ 0., 1., 0., 0. ],
          [ 1., 0., 0., 0. ],
          [ 0., 1., 0., 0. ],
          [ 0., 0., 0., 1. ],
        ]
        K = len(emission_dists)
        test_seq = 'accagtttgcact' # matches dist above
        test_seq_order_0 = seq_to_numpy(test_seq)

        # for various different orders
        for order in [1, 2]:

            # build a model of distribution above
            traits = PssmTraits(K, p_binding_site, order, num_background_states, create_background_model, emission_dists=emission_dists)
            model = traits.new_model()
            converted = hmm.model_states_2_model(model)
            B = converted.B

            # check the reverse complement states are correct
            for n in xrange(model.N):
                for o in xrange(model.M):
                    rev_comp_state, rev_comp_obs = traits.get_non_reverse_complement(n,o)
                    assert check_is_close_2(B[rev_comp_state,rev_comp_obs], B[n,o]), ('%d,%d %d,%d: %f %f' % (rev_comp_state,rev_comp_obs,n,o,B[rev_comp_state,rev_comp_obs],B[n,o]))

            # check viterbi gives correct result
            test_seq_order_n = converted.converter.to_order_n(test_seq_order_0)
            LL, states = converted.viterbi(test_seq_order_n)
            for i, state in enumerate(states):
                assert state == num_background_states+i


if __name__ == '__main__':
    unittest.main()
