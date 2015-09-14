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


root_dir = 'synthetic'
logger = setup_process(root_dir)

generator = SyntheticDataGenerator(
  K = 10,
  background_order = 0,
  num_background_mosaics = 1,
  background_model_creator = background_model
)
N = 200
lengths = [
        #20,
        #50,
        100,
        #500,
        #1000
]
dirichlet_strengths = [
        0.01,
        0.05,
        0.10,
        0.20,
        0.50,
]
num_seeds = 10

rocs = []

def roc_statistics(rocs):
    for statistic in [
            infpy.RocCalculator.positive_predictive_value,
            infpy.RocCalculator.performance_coefficient
    ]:
        values = [ statistic(roc) for roc in rocs ]
        yield statistic, (numpy.mean(values), st)
        mean =
        stddev =

for T in lengths:
    for dirichlet_strength in dirichlet_strengths:
        rocs.append([])
        for seed in xrange(1,num_seeds+1):

            # create directory for results
            directory = os.path.join(root_dir, '%d_%d' % (T, seed))
            if not os.access(directory, os.R_OK):
                os.makedirs(directory)

            # set up logging
            logger = logging.getLogger('find_in_synthetic.%d.%d' % (T, seed))
            logger.addHandler(logging.FileHandler(os.path.join(directory,'log.txt'), 'w'))
            logger.setLevel(logging.INFO)
            logger.info('**************** length: %d; seed: %d *****************' % (T, seed))

            # get the sequences
            traits, seqs, true_binding_sites = generator.generate(N, T, seed, directory, logger, dirichlet_strength=dirichlet_strength)
            num_bases = hmm.pssm.num_bases(seqs)
            num_known_bases = hmm.pssm.num_known_bases(seqs)
            logger.info('%d/%d (%d%%) known bases in %d sequences' % (num_known_bases, num_bases, (100*num_known_bases/num_bases), len(seqs)))

            # run the algorithm
            test_seqs = TestSequences(seqs)
            #import pdb; pdb.set_trace()
            algorithm = PssmFindingHmmTrainingAlgorithm(
                    traits,
                    use_viterbi=False,
                    directory=directory,
                    logger=logger,
                    tolerance=1.,
            )
            motif = algorithm(test_seqs)

            prediction = algorithm.make_prediction(motif, test_seqs)
            num_binding_bases = sum(sum(bs) for bs in prediction)
            logger.info('# binding bases: %d; about %d sites' % (num_binding_bases, num_binding_bases/generator.K))
            #import pdb; pdb.set_trace()

            # evaluate results
            predicted_binding_sites = algorithm.make_prediction(motif, test_seqs)
            nucleotide_roc = calculate_roc(true_binding_sites, predicted_binding_sites)
            logger.info('Nucleotide level roc:')
            logger.info(str(nucleotide_roc))
            site_roc = binding_site_level_roc(true_binding_sites, predicted_binding_sites)
            logger.info('Binding site level roc:')
            logger.info(str(site_roc))

            rocs[-1].append((nucleotide_roc, site_roc))
