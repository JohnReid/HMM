#
# Copyright John Reid 2007, 2008
#

"""
Code to test gapped PSSM traits.
"""



import hmm, hmm.pssm, unittest

class TestGappedPssmTraits(unittest.TestCase):
    def test_traits(self):
        from hmm.pssm import create_background_model, seq_to_numpy
        from infpy import check_is_close_2

        for order in [0, 1, 2]:
            num_background_mosaics = 4
            traits = hmm.pssm.GappedPssmTraits(
                    K=7,
                    p_binding_site=.1,
                    background_order=order,
                    num_background_mosaics=num_background_mosaics,
                    background_model_creator=create_background_model
            )
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
            m=hmm.as_model(traits.new_model(emission_dists))

            # check we have parameterised the reverse complement gaps correctly
            assert m.get_transition_parameterisation(29,28).idx == m.get_transition_parameterisation(14,15).idx
            assert m.get_transition_parameterisation(29,28).idx != m.get_transition_parameterisation(4,5).idx

            # check the reverse complement states are correct
            B = m.B
            for n in xrange(m.N):
                for o in xrange(m.M):
                    rev_comp_state, rev_comp_obs = traits.get_non_reverse_complement(n,o)
                    assert check_is_close_2(
                            B[rev_comp_state,rev_comp_obs],
                            B[n,o]
                    ), (
                            '%d,%d %d,%d: %f %f' % (
                                    rev_comp_state,rev_comp_obs,n,o,B[rev_comp_state,rev_comp_obs],B[n,o]
                            )
                    )

            # check viterbi gives correct result
            test_seq = 'acgtgat' # matches dist above
            test_seq_order_0 = hmm.pssm.seq_to_numpy(test_seq)
            test_seq_order_n = m.converter.to_order_n(test_seq_order_0)
            LL, states = m.viterbi(test_seq_order_n)
            for i, state in enumerate(states):
                assert (state-num_background_mosaics)/2 == i


if __name__ == '__main__':
    unittest.main()
