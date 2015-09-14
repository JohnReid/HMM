#
# Copyright John Reid 2007
#

"""
Code to generate synthetic data
"""

import random, numpy

class TestData(object):
    def __init__(self, seqs, binding_sites):
        self.seqs = seqs
        self.binding_sites = binding_sites

    def num_bases(self):
        return self.binding_sites.size

    def num_binding_bases(self):
        return self.binding_sites.sum()




def test_data_from_model(num_seqs, length, model, model_traits, p_unknown=0.02):
    "Generate so many sequences of given length from the model"
    seqs = []
    true_states = []
    for i in xrange(num_seqs):
        states, seq = model.sample(length)
        seq.dtype = numpy.int32
        seqs.append(seq)
        true_states.append(states)
        seq.setflags(write = True)
        for t in xrange(len(seq)):
            if random.random() < p_unknown: seq[t] = model.converter.order_0_size # set as unknown
    true_binding = numpy.array([[s not in model_traits.background_states for s in true_state] for true_state in true_states])
    return TestData(seqs, true_binding)
