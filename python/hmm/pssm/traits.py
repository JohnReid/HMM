#
# Copyright John Reid 2007, 2008
#

"""
Base traits for gapped and ungapped PSSM HMM models.
"""

import hmm, numpy, pickle

class BaseTraits(object):
    '''
    Base for things that build HMM models to represent PSSMs
    '''

    def __init__(
            self,
            K,
            p_binding_site,
            background_order,
            num_background_mosaics,
            background_model_creator,
            prior_strength
    ):
        self.K = K
        self.p_binding_site = p_binding_site
        self.order = background_order
        self.background_model_creator = background_model_creator
        self.num_background_mosaics = num_background_mosaics
        self.prior_strength = prior_strength

    def write_logo(self, model, f):
        import hmm.pssm.logo as logo
        dist = self.pssm_dist(model)
        image = logo.pssm_as_image(dist)
        image.save(f, "PNG")
        return image

    def log_info(self, model, log):
        model = hmm.as_model(model)
        log.info('IC: %.3f' % hmm.pssm.information_content(self.pssm_dist(model)))
        log.info('p(binding site): %.6f' % self.p_binding_site_for_model(model))

    def prior(self, model):
        prior = hmm.ModelPrior(self.N(), self.M())
        prior.A = model.A * self.prior_strength
        prior.B = model.B * self.prior_strength
        prior.pi = model.pi * self.prior_strength
        return prior

    def binding_sites_from_states(self, states):
        return numpy.array([(state not in self.background_states and 1 or 0) for state in states])

    def get_p_binding(self, model):
        return 2.0 * hmm.as_model(model).A[0,self.num_background_mosaics]

    def set_p_binding(self, model, p):
        model.set_transition(0, self.num_background_mosaics, p/2)
        model.set_transition(0, 0, 1.0 - p)
        model.normalise()
