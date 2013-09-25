#
# Copyright John Reid 2007
#

"""
Support for numerical optimisation on HMMs
"""

import math, scipy.optimize, numpy


_largest = 2.**1023
_log_largest = math.log(_largest)
_smallest = -_largest

def sigmoid( x ):
    """
    y = 1 / ( 1 + exp( -x ) )
    """
    if _largest <= x:
        return 1.0
    if -_log_largest >= x:
        return 0.0
    return 1.0 / ( 1.0 + math.exp( -x ) )

def inverse_sigmoid( y ):
    """
    x = ln( y / ( 1 - y ) )
    """
    if 1.0 == y:
        return _largest
    if 0.0 == y:
        return -_log_largest
    return math.log( y ) - math.log( 1 - y )


class OptimisableModel(object):
    def __init__(self, model, param_indices, seqs, callback=None):
        self.model = model
        self.param_indices = param_indices
        self.seqs = seqs
        self.callback = callback

    def __call__(self, params):
        'The optimiser calls this'
        try:
            self.set_model_params(params)
        except:
            return 1e300 #if we cannot set the parameters return a large value
        LL = sum(self.model.LL(seq) for seq in self.seqs)
        #print 'LL=%f' % LL
        if None != self.callback:
            self.callback(LL)
        return -LL

    def get_model_params(self):
        return numpy.array([inverse_sigmoid(self.model.theta[idx]) for idx in self.param_indices])

    def set_model_params(self, params):
        assert len(self.param_indices) == len(params)
        for p, idx in zip(params, self.param_indices):
            self.model.theta[idx] = sigmoid(p)
        self.model.normalise()

if '__main__' == __name__:
    import hmm.pssm, numpy
    import scipy.optimize as O
    model = hmm.as_model(hmm.pssm.create_background_model(0, 1))
    param_idx_map = [ model.get_emission_parameterisation(0, k).idx for k in xrange(4) ]
    seqs = [
      numpy.array([ 0,1,2,3,0,1,2]),
      numpy.array([ 0,1,2,3,0,1,2,3]),
    ]
    optimisable = OptimisableModel(model, param_idx_map, seqs)
    xopt, fopt, iter, funcalls, warnflag = O.fmin(
            optimisable,
            optimisable.get_model_params(),
            full_output=True
    )
    optimisable.set_model_params(xopt)
