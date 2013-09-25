#
# Copyright John Reid 2007
#


"""
Implementation of PSSM finding algorithm - supercedes that in test_synthetic.py
"""


import hmm, hmm.pssm, numpy, os, os.path, logging, random, sys, math, traceback

from sequences import *
from find_motif import *
from background_models import *
from process_util import *


root_dir = 'fragments'
logger = setup_process(root_dir)

K = 10
background_order = 0
num_background_mosaics = 1
exp_sites_per_seq = 1.0

for fragment in all_fragments:

    # create directory for results
    directory = os.path.join(root_dir, '%s' % fragment)
    if not os.access(directory, os.R_OK):
        os.makedirs(directory)

    # set up logging
    logger = logging.getLogger('find_in_fragments.%s' % fragment)
    logger.addHandler(logging.FileHandler(os.path.join(directory,'log.txt'), 'w'))
    logger.setLevel(logging.INFO)
    logger.info('**************** %s *****************' % fragment)

    hmm.seed_rng(1)
    random.seed(1)

    # get the sequences
    seqs = seqs_for_fragment(fragment)
    num_bases = hmm.pssm.num_bases(seqs)
    num_known_bases = hmm.pssm.num_known_bases(seqs)
    logger.info('%d/%d (%d%%) known bases in %d sequences' % (num_known_bases, num_bases, (100*num_known_bases/num_bases), len(seqs)))

    p_binding_site = exp_sites_per_seq*len(seqs)/float(num_bases) # one binding site per sequence on average
    logger.info('%.2f expected sites per sequence gives p(binding site)=%f' % (exp_sites_per_seq, p_binding_site))
    def per_fragment_background_model(order, N):
        return hmm.as_state_model(global_background_model_cache().load((order, N, fragment)))
    traits = hmm.pssm.GappedPssmTraits(
            K=K,
            p_binding_site=100*p_binding_site,
            background_order=background_order,
            num_background_mosaics=num_background_mosaics,
            prior_strength = 1.0,
            p_gap = 0.16,
            #background_model_creator=background_model
            background_model_creator=per_fragment_background_model
    )

    # run the algorithm
    test_seqs = TestSequences(seqs)
    prior = hmm.ModelPrior(traits.N(), traits.M())
    prior.A = prior.A + 100000 * traits.background_to_binding_transition_prior(10*p_binding_site, 0.01)
    #import pdb; pdb.set_trace()
    algorithm = PssmFindingHmmTrainingAlgorithm(
            traits,
            use_viterbi=True,
            directory=directory,
            logger=logger,
            tolerance=.001,
            prior=prior
    )
    motif = algorithm(test_seqs)

    prediction = algorithm.make_prediction(motif, test_seqs)
    num_binding_bases = sum(sum(bs) for bs in prediction)
    logger.info('# binding bases: %d; about %d sites' % (num_binding_bases, num_binding_bases/K))
    #import pdb; pdb.set_trace()
