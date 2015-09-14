#
# Copyright John Reid 2007
#


"""
To test problem with viterbi training
"""


import hmm, pickle

model, seqs = pickle.load(open('viterbi_problem.pickle'))
def callback(LL):
    print LL
LL, iterations = model.viterbi_training(
        seqs,
        #prior = self.traits.prior(model),
        tolerance = 1e-99,
        callback = callback,
        max_iterations = 0
)
