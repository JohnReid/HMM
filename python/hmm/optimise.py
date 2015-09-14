#
# Copyright John Reid 2007
#

'''
Code that helps to interface HMMs to various optimisation techniques
'''


import numpy

def parameter_transform(model):
    '''
    The logarithm of the model's parameters.
    '''
    return numpy.log(model.theta)

def parameterise_model(model, params):
    '''
    Takes the logarithm of the parameters and reparameterises the model with them
    '''
    for i, p in enumerate(params):
        model.theta[i] = numpy.exp(p)
    model.normalise()


class LL_calculator(object):
    '''
    Callable object that calculates -LL of the sequences it is initialised with
    given model parameters
    '''

    def __init__(model, seqs):
        self.model = hmm.Model(model)
        self.seqs = seqs
        self.known_bases = sum(self.model.converter.num_known_bases_order_n(seq) for seq in seqs)

    def __call__(params):
        parameterise_model(self.model, params)
        LL = sum(self.model.forward(seq)[0] for seq in self.seqs)
        return -LL / self.known_bases
