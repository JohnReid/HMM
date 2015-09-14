#
# Copyright John Reid 2007
#

"""
Evaluate different sizes and orders of background models
"""

import hmm, hmm.pssm, numpy, os, time, random

from fragments import *
from sequences import *
from train import *

def train_and_validate(
        model,
        training_sequences,
        validation_sequences
):
    # convert to order n sequences
    training_sequences = [ model.converter.to_order_n(s) for s in training_sequences ]
    validation_sequences = [ model.converter.to_order_n(s) for s in validation_sequences ]

    # train on test
    def bw_callback(LL): pass; #print 'LL: %.3f' % LL
    tolerance = 1e-4 * num_known_bases(training_sequences)
    model.baum_welch(
            training_sequences,
            tolerance = tolerance,
            callback = bw_callback
    )

    # test on validation
    LL = sum( [ model.forward(s)[0] for s in validation_sequences ] )
    known_bases = sum( [ model.converter.num_known_bases_order_n(s) for s in validation_sequences ] )
    return LL / known_bases


def index_filter(predicate, iterable):
    for i, x in enumerate(iterable):
        if predicate(i, x):
            yield x

def k_fold_cross_validation(X, K, randomise = False):
    """
    Yields (training, validation) which are themselves iterables over a partition
    of X into a (k-1)/K*len(X) training set and a len(X)/K validation set

    E.g.
    X = [i for i in xrange(97)]
    for training, validation in k_fold_cross_validation(X, K=7):
            for x in X: assert (x in training) ^ (x in validation), x
    """
    if randomise: from random import shuffle; X=list(X); shuffle(X)
    for k in xrange(K):
        training = list( index_filter( lambda i,x: i % K != k, X ) )
        validation = list( index_filter( lambda i,x: i % K == k, X ) )
        yield training, validation

try:
    import cookbook
    cookbook.make_current_process_nice()
except:
    print 'Could not set process priority'


max_order = 3
max_num_mosaics = 6
num_test_sets = 5
LLs = numpy.zeros( (max_order+1, max_num_mosaics) )
K = 3
max_num_seqs = 600

if True: # pre-calculated
    LLs = numpy.array([[-1.38182726, -1.36265837, -1.36025689, -1.35821224, -1.35728875,
    -1.3567072 ],
   [-1.36292533, -1.33791934, -1.33418864, -1.33341424, -1.33235487,
    -1.33175728],
   [-1.35585918, -1.3342389 , -1.33090382, -1.3302001 , -1.32931403,
    -1.32898274],
   [-1.35164935, -1.33259588, -1.32956143, -1.32924632, -1.32886239,
    -1.3284187 ]])
else:
    try: seqs
    except:
        seqs = get_all_sequences()
        from random import shuffle
        shuffle(seqs)
        seqs = seqs[:max_num_seqs]
    for order in xrange(max_order+1):
        model_builder = hmm.pssm.ModelBuilder(order)
        for num_mosaics in xrange(1, max_num_mosaics+1):
            LL = 0.0
            for training, validation in k_fold_cross_validation(seqs, K):
                model_by_states = model_builder.create_background_mosaic_model(num_mosaics, 0.01, 100.0)
                model = hmm.model_states_2_model(model_by_states)
                LL_per_base = train_and_validate(model, training, validation)
                print 'Order: %d     # mosaics: %d      LL/known base: %.3f' % (order, num_mosaics, LL_per_base)
                LL += LL_per_base
            LLs[order, num_mosaics-1] = LL

import pylab
fig = pylab.figure()
pylab.xlabel('# mosaics')
pylab.ylabel('LL/base')
for order, l in enumerate(LLs):
    pylab.plot(xrange(1,len(l)+1), l, label='order: %d'%order)
pylab.legend(loc=4)
pylab.savefig('LLs.png')
