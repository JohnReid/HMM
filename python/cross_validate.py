#
# Copyright John Reid 2007
#

"""
Code to cross-validate MEME and our algorithm on test sequences
"""

from run_meme import *
from process_util import *
from fragments import *
from background_models import *
from find_motif import *


def training_set_filename(fragment='T00594', idx=0):
    return 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-trainingset-%d.fasta' % (fragment, fragment, idx)

def test_set_filename(fragment='T00594', idx=0):
    return 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-testset-%d.fasta' % (fragment, fragment, idx)

def meme_filename(fragment='T00594', idx=0):
    return 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-meme-output-%d.txt' % (fragment, fragment, idx)

class CrossValidationTest(object):
    def __init__(self, fragment='T00594', idx=0):
        self.fragment=fragment
        self.idx=idx
        self.training_set_filename = 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-trainingset-%d.fasta' % (fragment, fragment, idx)
        self.test_set_filename = 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-testset-%d.fasta' % (fragment, fragment, idx)
        self.meme_filename = 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-meme-output-%d.txt' % (fragment, fragment, idx)
        self.meme_logo_filename = 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-meme-logo-%d.png' % (fragment, fragment, idx)
        self.hmm_logo_filename = 'C:/Analysis/GappedPssms/fragments/test-sets/%s/%scorrectedRM-hmm-logo-%d.png' % (fragment, fragment, idx)

    def training_seqs(self):
        if not hasattr(self, '_training_seqs'):
            self._training_seqs = [ seq_to_numpy(seq) for seq in convert_seqs(test.training_set_filename) if len(seq) ]
        return self._training_seqs

    def test_seqs(self):
        if not hasattr(self, '_test_seqs'):
            self._test_seqs = [ seq_to_numpy(seq) for seq in convert_seqs(test.test_set_filename) if len(seq) ]
        return self._test_seqs


def model_LL(model, test):
    seqs = test.test_seqs()
    LL = sum([model.forward(model.converter.to_order_n(seq))[0] for seq in seqs])
    logging.getLogger().info('LL=%f' % LL)
    return LL

def meme_cross_validate(test):
    "Returns LL of test sequences"
    from TAMO.MD.MEME import Meme
    import hmm.pssm, logging
    if not os.access(test.meme_filename, os.R_OK):
        logging.getLogger().info('Running MEME')
        run_meme(test.training_set_filename, test.meme_filename, maxiter=150, minw=None, maxw=None, extra_args = ['-text', '-w 12'])
    meme = Meme(test.meme_filename)
    motif = meme.motifs[0]
    dist = meme_motif_dist(motif)
    K = len(dist)
    logging.getLogger().info('MEME found PSSM of size: K=%d' % K)
    p_binding_site_estimate = len(motif.seqs)/float(sum(len(seq) for seq in test.training_seqs()))
    logging.getLogger().info('Estimating p(binding)=%f' % p_binding_site_estimate)
    traits = hmm.pssm.PssmTraits(
      K=K,
      p_binding_site=p_binding_site_estimate,
      background_order=background_order,
      num_background_mosaics=num_background_mosaics,
      background_model_creator=background_model,
      emission_dists=dist,
    )
    model=hmm.as_model(traits.new_model())
    traits.write_logo(model, test.meme_logo_filename)
    return model_LL(model, test)

def hmm_cross_validate(
        test,
        algorithm,
):
    "Returns LL ratio between MEME and HMM over test sequences"
    import hmm.pssm, logging
    motif = algorithm(TestSequences(test.training_seqs()))
    model = algorithm.model_from_motif(motif)
    gapped_pssm_traits.write_logo(model, test.hmm_logo_filename)
    return model_LL(model, test)

if '__main__' == __name__:
    directory = 'cross-validation'
    logger = setup_process(directory)

    #test = CrossValidationTest()
    #print meme_cross_validate(test)

    K=12
    p_binding_site=.01
    background_order=2
    num_background_mosaics=4
    prior_strength = 1.0
    background_model_creator=background_model
    gapped_pssm_traits = hmm.pssm.GappedPssmTraits(
            K=K,
            p_binding_site=p_binding_site,
            background_order=background_order,
            p_gap = 0.1,
            num_background_mosaics=num_background_mosaics,
            prior_strength = prior_strength,
            background_model_creator=background_model
    )

    for fragment in all_fragments:
        for idx in xrange(5):
            logging.getLogger().info('*************** %s: %d ****************' % (fragment, idx))
            test = CrossValidationTest(fragment, idx)
            #meme_LL = meme_cross_validate(test)
            #logging.getLogger().info('MEME LL: %f' % (meme_LL))
            algorithm = PssmFindingAlgorithm(gapped_pssm_traits, directory='%s/%s_%d' % (directory, fragment, idx), logger=logger)
            hmm_LL = hmm_cross_validate(test, algorithm)
            logging.getLogger().info('HMM LL: %f' % (hmm_LL))
            #logging.getLogger().info('LL ratio=%f' % (hmm_LL - meme_LL))
