#
# Copyright John Reid 2007, 2008
#

"""
Traits of Gapped PSSMs modelled by HMMs. I{Is this superceded now?}
"""

import hmm, numpy, pickle
from traits import BaseTraits

class GappedPssmTraits(BaseTraits):
    '''
    Builds HMM models to represent gapped PSSMs
    '''

    def name(self):
        return 'gapped_pssm'

    def __init__(
            self,
            K,
            p_binding_site,
            background_order,
            num_background_mosaics,
            background_model_creator,
            emission_dists = None,
            p_gap = 0.35,
            gap_threshold = 0.1,
            ic_threshold = 0.5,
            prior_strength = 0.0
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
        self.K_prime = 2 * self.K - 1
        self.model_builder = hmm.pssm.ModelBuilder(self.order)
        self.background_states = range(num_background_mosaics)
        self.reverse_complements = dict(
          (state, state - self.K_prime)
          for state in
                xrange(num_background_mosaics + self.K_prime, num_background_mosaics + 2*self.K_prime)
        )
        self.emission_dists = emission_dists
        self.p_gap = p_gap
        self.gap_threshold = gap_threshold
        self.ic_threshold = ic_threshold

    def kth(self, k, rev_comp=False):
        "If the bases of the pssm are 0..K-1, return the idx of the k'th state"
        if k > self.K or k < 0:
            raise RuntimeError('k must be between 0 and K-1')
        result = self.num_background_mosaics + 2*k
        if rev_comp:
            result += 2*self.K-1
        return result

    def kth_gap(self, k, rev_comp=False):
        "If the bases of the pssm are 0..K-1, return the idx of the k'th gap state"
        if k > self.K-1 or k < 0:
            raise RuntimeError('k must be between 0 and K-2')
        result = self.num_background_mosaics + 2*k + 1
        if rev_comp:
            result += 2*self.K-1
        return result

    def gap_transition_prior(self, k, prior_strength, p_gap=.5):
        'Creates a matrix for use as a prior in various training algorithms'
        A = numpy.zeros((self.N(), self.N()))
        gap_prior = p_gap * prior_strength
        non_gap_prior = (1.0-p_gap) * prior_strength
        kth = self.kth(k)
        kth_gap = self.kth_gap(k)
        kth_plus_1 = self.kth(k+1)
        A[kth, kth_gap] = gap_prior
        A[kth, kth_plus_1] = non_gap_prior
        A[kth+self.K_prime, self.kth_gap+self.K_prime] = gap_prior
        A[kth+self.K_prime, kth_plus_1+self.K_prime] = non_gap_prior
        return A

    def emission_dist_including_n_mer(self, n_mer, strength=9.0, offset=0, dirichlet_strength = 10.0):
        '''
        Get an emission dist that incorporates the n-mer

        Draws a distribution from a dirichlet of given strength (default 10) then adds the n-mer in the middle
        (default # copies of n-mer = 9)
        '''
        emission_dists = [ hmm.dirichlet_draw(numpy.ones(4)*dirichlet_strength) for k in xrange(self.K_prime) ]
        n = len(n_mer)
        offset = (self.K - n)/2
        for i, base in enumerate(n_mer):
            idx = 2*(i+offset)
            emission_dists[idx].setflags(write=True)
            emission_dists[idx][base] += strength
            emission_dists[idx] /= emission_dists[idx].sum()
        return emission_dists

    def new_model(self, emission_dists = None, dirichlet_strength = 10.0, generate_uniform=False):
        model = self.background_model_creator(self.order, self.num_background_mosaics)

        # are we generating a uniform emission model?
        if generate_uniform:
            emission_dists = [ .25 * numpy.ones(4) for k in xrange(self.K_prime) ]
        else:
            # have we been given any emission dist?
            if None == emission_dists:
                # no - were we initialised with some distribution?
                if None == self.emission_dists:
                    # no - create a random one
                    emission_dists = [ hmm.dirichlet_draw(numpy.ones(4)*dirichlet_strength) for k in xrange(self.K_prime) ]
                else:
                    # yes - use it
                    emission_dists = self.emission_dists
        if len(emission_dists) != self.K_prime:
            raise RuntimeError('Wrong number of emissions')

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
        gap_params = [ model.add_parameter(self.p_gap) for i in xrange(self.K-1) ]
        non_gap_params = [ model.add_parameter(1.0-self.p_gap) for i in xrange(self.K-1) ]
        for k in xrange(self.K-1):
            positive_states[2*k+1].add_successor(positive_states[2*k+2], one_param)
            negative_states[2*k+1].add_successor(negative_states[2*k+2], one_param)
            positive_states[2*k].add_successor(positive_states[2*k+1], gap_params[k])
            negative_states[2*k].add_successor(negative_states[2*k+1], gap_params[self.K-2-k])
            positive_states[2*k].add_successor(positive_states[2*k+2], non_gap_params[k])
            negative_states[2*k].add_successor(negative_states[2*k+2], non_gap_params[self.K-2-k])

        # connect the last states back to the background
        one_param = model.add_parameter(1.0/self.num_background_mosaics)
        for bg in xrange(self.num_background_mosaics):
            positive_states[-1].add_successor(model.states[bg], one_param)
            negative_states[-1].add_successor(model.states[bg], one_param)

        return model

    def N(self):
        return self.num_background_mosaics + 2 * self.K_prime

    def M(self):
        return self.model_builder.converter.order_n_size

    def prior(self, model):
        tmp_model = hmm.as_model(self.new_model( generate_uniform=True ))
        prior = hmm.ModelPrior(self.N(), self.M())
        prior.A = tmp_model.A * self.prior_strength
        prior.B = tmp_model.B * self.prior_strength
        prior.pi = tmp_model.pi * self.prior_strength
        return prior

    def pssm_dist(self, model):
        return hmm.as_model(model).B[self.num_background_mosaics:self.num_background_mosaics+self.K_prime,:4]

    def information_content_per_base(self, model):
        model = hmm.as_model(model)
        emissions = self.pssm_dist(model)
        result = 0.0
        expected_num_bases = 0.0
        for k in xrange(self.K):
            result += hmm.pssm.base_information_content(emissions[2*k])
            expected_num_bases += 1.0
            if k < self.K - 1:
                p_gap = self.p_gap_for_model(model, k)
                expected_num_bases += p_gap
                result += p_gap * hmm.pssm.base_information_content(emissions[2*k+1])
        return result / expected_num_bases

    def write_logo(self, model, f, rev_comp=False):
        import hmm.pssm.logo as logo
        model = hmm.as_model(model)
        emissions = self.pssm_dist(model)
        transparencies = []
        pssm_dist = []
        for k in xrange(self.K):
            pssm_dist.append(emissions[2*k])
            transparencies.append(1.0)
            if k < self.K - 1:
                p_gap = self.p_gap_for_model(model, k)
                if p_gap > self.gap_threshold:
                    pssm_dist.append(emissions[2*k+1])
                    transparencies.append(p_gap)
        if rev_comp:
            pssm_dist.reverse()
            for i, emission in enumerate(pssm_dist):
                pssm_dist[i] = emission[::-1]
            transparencies.reverse()
        image = logo.pssm_as_image(pssm_dist, transparencies=transparencies)
        image.save(f, "PNG")
        return image

    def p_gap_for_model(self, model, k):
        return hmm.as_model(model).A[self.kth(k), self.kth_gap(k)]

    def log_info(self, model, log):
        model = hmm.as_model(model)
        BaseTraits.log_info(self, model, log)
        log.info('Adjusted IC/base = %.4f' % self.information_content_per_base(model))
        A = model.A
        for k in xrange(self.K):
            p_gap = self.p_gap_for_model(model, k)
            if p_gap > self.gap_threshold:
                log.info('p(gap at %d) = %.4f' % (k, p_gap))

    def alternative_generators(self, model):
        return [ ]

    def is_model_interesting(self, model):
        '''
        Takes a learnt model and determines if it is interesting enough.
        '''
        return self.information_content_per_base(model) > self.ic_threshold

    def create_test_data_generating_model(self, dirichlet_strength=.05):
        data_generating_model = hmm.as_model(self.new_model(dirichlet_strength = .05))
        for k in xrange(self.K-1):
            data_generating_model.set_transition(self.kth(k),self.kth_gap(k),0.0) # set all gaps to 0.0
        data_generating_model.set_transition(self.kth(self.K/2),self.kth_gap(self.K/2),0.5) # except for gap at midpoint
        for i in xrange(data_generating_model.N):
            data_generating_model.set_initial(i, 0.0)
        for i in xrange(self.num_background_mosaics):
            data_generating_model.set_initial(i, 1.0 / self.num_background_mosaics)
        data_generating_model.normalise()
        return data_generating_model

    def background_to_binding_transition_prior(self, p_binding_site, p_alt_background):
        '''
        Return a transition matrix that can be used as a prior for the transition
        between background states and binding site states
        '''
        result = numpy.zeros((self.N(), self.N()))
        if self.num_background_mosaics > 1:
            p_remain_in_bg = 1.0 - p_binding_site/2 - p_alt_background/(self.num_background_mosaics-1)
        else:
            p_remain_in_bg = 1.0 - p_binding_site
        for bg1 in self.background_states:
            for bg2 in self.background_states:
                if bg1 == bg2:
                    result[bg1,bg2] = p_remain_in_bg
                else:
                    result[bg1,bg2] = p_alt_background/(self.num_background_mosaics-1)
            result[bg1, self.kth(0,rev_comp=False)] = p_binding_site/2
            result[bg1, self.kth(0,rev_comp=True)] = p_binding_site/2
        return result


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
    K = 7
    test_seq = 'acgtgat' # matches dist above
    test_seq_order_0 = hmm.pssm.seq_to_numpy(test_seq)

    # for various different orders
    for order in [0, 1, 2]:

        # build a model of distribution above
        traits = GappedPssmTraits(K, p_binding_site, order, num_background_states, create_background_model, emission_dists=emission_dists)
        model = traits.new_model()
        converted = hmm.model_states_2_model(model)
        B = converted.B
        hmm.graph_as_svg(converted, 'gapped_pssm')


        # check the reverse complement states are correct
        for n in xrange(model.N):
            for o in xrange(model.M):
                rev_comp_state, rev_comp_obs = traits.get_non_reverse_complement(n,o)
                assert check_is_close_2(B[rev_comp_state,rev_comp_obs], B[n,o]), ('%d,%d %d,%d: %f %f' % (rev_comp_state,rev_comp_obs,n,o,B[rev_comp_state,rev_comp_obs],B[n,o]))

        # check viterbi gives correct result
        test_seq_order_n = converted.converter.to_order_n(test_seq_order_0)
        LL, states = converted.viterbi(test_seq_order_n)
        for i, state in enumerate(states):
            assert (state-num_background_states)/2 == i
