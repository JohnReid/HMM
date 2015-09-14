#
# Copyright John Reid 2007, 2008
#

"""
Traits of PSSMs modelled by HMMs.
"""

import hmm, numpy, pickle
from traits import BaseTraits

class PssmTraits(BaseTraits):
    '''
    Builds HMM models to represent PSSMs
    '''

    def name(self):
        return 'pssm'

    def __init__(
            self,
            K,
            p_binding_site,
            background_order,
            num_background_mosaics,
            background_model_creator,
            emission_dists = None,
            prior_strength = 0.0,
            interesting_ic_per_base = .3
    ):
        BaseTraits.__init__(
                self,
                K,
                p_binding_site,
                background_order,
                num_background_mosaics,
                background_model_creator,
                prior_strength = prior_strength
        )
        self.model_builder = hmm.pssm.ModelBuilder(self.order)
        self.background_states = set(xrange(num_background_mosaics))
        self.reverse_complements = set()
        self.reverse_complements = dict(
          (state, state - self.K)
          for state in
                xrange(num_background_mosaics + K, num_background_mosaics + 2*K)
        )
        self.emission_dists = emission_dists
        self.prior_strength = prior_strength
        self.interesting_ic_per_base = interesting_ic_per_base
        if None != emission_dists and self.K != len(emission_dists):
            raise RuntimeError('K and len(emission_dists) mismatch')

    def new_model(self, emission_dists = None, dirichlet_strength = 10.0):
        model = self.background_model_creator(self.order, self.num_background_mosaics)

        # have we been given any emission dist?
        if None == emission_dists:
            # no - were we initialised with some distribution?
            if None == self.emission_dists:
                # no - create a random one
                emission_dists = [ hmm.dirichlet_draw(numpy.ones(4)*dirichlet_strength) for k in xrange(self.K) ]
            else:
                # yes - use it
                emission_dists = self.emission_dists

        # add the states for the positive strand orientation
        positive_states = [
          self.model_builder.add_order_0_parameterised_state(
                        model,
                        pi=self.p_binding_site,
                        emission_dist=dist)
                for dist in emission_dists
        ]

        # add the states for the negative strand orientation
        negative_states = [
                self.model_builder.add_order_0_rev_comp_state(
                        model,
                        positive_state,
                        pi=self.p_binding_site)
                for positive_state in positive_states
        ]
        negative_states.reverse()

        # connect the background to the both first states with equal prob
        param = model.add_parameter(self.p_binding_site/2)
        for bg in xrange(self.num_background_mosaics):
            model.states[bg].add_successor(positive_states[0], param)
            model.states[bg].add_successor(negative_states[0], param)

        # connect the states in the pssm together
        one_param = model.add_parameter(1.0)
        for k in xrange(self.K-1):
            positive_states[k].add_successor(positive_states[k+1], one_param)
            negative_states[k].add_successor(negative_states[k+1], one_param)

        # connect the last states back to the background
        one_param = model.add_parameter(1.0/self.num_background_mosaics)
        for bg in xrange(self.num_background_mosaics):
            positive_states[-1].add_successor(model.states[bg], one_param)
            negative_states[-1].add_successor(model.states[bg], one_param)

        return model

    def N(self):
        return self.num_background_mosaics + 2 * self.K

    def M(self):
        return self.model_builder.converter.order_n_size

    def pssm_dist(self, model):
        return hmm.as_model(model).B[self.num_background_mosaics:self.num_background_mosaics+self.K,:4]

    def shift_left(self, model):
        # print 'Shifting left'
        alt_model = hmm.Model(model) # generate a copy
        B = alt_model.B
        # shift the dists
        for k in xrange(self.num_background_mosaics, self.num_background_mosaics + self.K-1):
            alt_model.set_emission_distribution(k, B[k+1])
        # set the last to be uniform
        alt_model.set_emission_distribution(self.num_background_mosaics + self.K-1, numpy.ones(alt_model.M) / 4)
        alt_model.normalise()
        return alt_model

    def shift_right(self, model):
        # print 'Shifting right'
        alt_model = hmm.Model(model) # generate a copy
        B = alt_model.B
        # shift the dists
        for k in xrange(self.num_background_mosaics+1, self.num_background_mosaics + self.K):
            alt_model.set_emission_distribution(k, B[k-1])
        # set the last to be uniform
        alt_model.set_emission_distribution(self.num_background_mosaics, numpy.ones(alt_model.M) / 4)
        alt_model.normalise()
        return alt_model

    def alternative_generators(self, model):
        # look at information content of first and last base in pssm distribution
        # to choose whether to return shift_left or shift_right
        pssm_dist = self.pssm_dist(model)
        if hmm.pssm.information_content(pssm_dist[0:]) > hmm.pssm.information_content(pssm_dist[-1:]):
            return [ self.shift_right ]
        else:
            return [ self.shift_left ]

    def is_model_interesting(self, model):
        '''
        Takes a learnt model and determines if it is interesting enough.
        '''

        # do we have more than .5 information content per base?
        return hmm.pssm.information_content(self.pssm_dist(model)) > self.interesting_ic_per_base * self.K

    def create_test_data_generating_model(self, dirichlet_strength=.05):
        data_generating_model = hmm.as_model(self.new_model(dirichlet_strength = .05))
        for i in xrange(data_generating_model.N):
            data_generating_model.set_initial(i, 0.0)
        for i in xrange(self.num_background_mosaics):
            data_generating_model.set_initial(i, 1.0 / self.num_background_mosaics)
        data_generating_model.normalise()
        return data_generating_model

if '__main__' == __name__:
    from hmm.pssm import create_background_model
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
    for order in [0, 1, 2]:

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
