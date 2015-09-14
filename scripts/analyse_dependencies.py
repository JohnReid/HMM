#
# Copyright John Reid 2007
#

"""
Code to analyse dependencies between positions in sites
"""


import hmm, hmm.pssm, pickle

def load_sites(f):
    if isinstance(f, str):
        f = open(f)
    for l in f:
        if l.startswith('#'):
            continue
        site, sequence_idx, start, end, rev_comp = l.strip().split(';')
        yield site, sequence_idx, start, end, rev_comp

fragment = 'T00140'
pssm = 2

model_filename = 'models/%s/gapped_pssm/pssm_%d_model.pickle' % (fragment, pssm)
model = pickle.load(open(model_filename))

from sequences import *
seqs=seqs_for_fragment('T00140')
seq_2_order_0 = seqs[2]
seq_2_order_n = model.converter.to_order_n(seq_2_order_0)
for site, start, i in hmm.pssm.find_sites(seq_2_order_n, seq_2_order_0, model, xrange(4)):
    print site
    print start
    print i
viterbi_seq_2 = model.viterbi(seq_2)[1]

traits = hmm.pssm.GappedPssmTraits(12, .01, 2, 4, None)
learner = hmm.pssm.PssmLearner([seq_2_order_0], traits)
total_sites, sites, counts, model = learner.analyse_sites(model)

raise ''

for i in xrange(model.N):
    print i
    model.set_initial(i, 0.001)
model.set_initial(4,1.0)
#model.set_transition(3,4,.9)
model.normalise()

sites_filename = 'models/%s/gapped_pssm/pssm_%d_sites.txt' % (fragment, pssm)
sites = [ 'aaaa' + site[0] for site in load_sites(sites_filename) ]

numpy_sites = [ hmm.pssm.seq_to_numpy(seq) for seq in sites ]
numpy_order_n_sites = [ model.converter.to_order_n(s)[4:] for s in numpy_sites ]
viterbi_states = [ model.viterbi(s)[1] for s in numpy_order_n_sites ]
