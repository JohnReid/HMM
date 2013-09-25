#
# Copyright John Reid 2007
#


"""
Generates synthetic data for tests
"""


import numpy, hmm, random, os, cPickle
from sequences import *

def test_seq_from_model(length, model, traits):
    "Returns seq, binding_sites"
    states, seq = model.sample(length)
    seq.dtype = numpy.int32
    return seq, traits.binding_sites_from_states(states)


def test_data_from_model(num_seqs, length, model, traits):
    "Returns two lists seqs, binding_sites"
    seqs = []
    binding_sites = []
    while len(seqs) < num_seqs:
        seq, sites = test_seq_from_model(length, model, traits)
        seqs.append(seq)
        binding_sites.append(sites)
    return seqs, binding_sites

def place_random_unknowns(seqs, unknown_base, p_unknown=.02):
    import random
    for seq in seqs:
        seq.setflags(write=True)
        for t in xrange(len(seq)):
            if random.random() < p_unknown:
                seq[t] = unknown_base # set as unknown

def fasta_filename(directory, T, seed):
    return os.path.join(directory, 'T-%d_seed-%d.fa' % (T, seed))

def binding_sites_filename(directory, T, seed):
    return os.path.join(directory, 'T-%d_seed-%d-sites.pickle' % (T, seed))

def save_synthetic_data(directory, T, seed, seqs, true_binding_sites):
    # save the data
    write_seqs_as_fasta(seqs, fasta_filename(directory, T, seed))
    cPickle.dump(true_binding_sites, open(binding_sites_filename(directory, T, seed), 'w'))

def load_synthetic_data(directory, T, seed):
    return (
seqs_to_numpy(convert_seqs(fasta_filename(directory, T, seed))),
cPickle.load(open(binding_sites_filename(directory, T, seed)))
    )

class SyntheticDataGenerator(object):
    def __init__(
            self,
            K,
            background_order,
            num_background_mosaics,
            background_model_creator
    ):
        self.K = K
        self.background_order = background_order
        self.num_background_mosaics = num_background_mosaics
        self.background_model_creator = background_model_creator

    def generate(
            self,
            N,
            T,
            seed,
            directory,
            logger,
            p_binding_site = None,
            p_unknown = 0.02,
            dirichlet_strength = 0.05
    ):
        hmm.seed_rng( T+seed )
        random.seed( T+seed )

        if None == p_binding_site:
            p_binding_site = 1.0/T # one binding site per sequence on average by default

        # traits that determine how our model maps to a PSSM
        traits = hmm.pssm.GappedPssmTraits(
                K=self.K,
                p_binding_site = p_binding_site,
                background_order = self.background_order,
                num_background_mosaics = self.num_background_mosaics,
                background_model_creator = self.background_model_creator
        )

        # generate some data
        data_generating_model = traits.create_test_data_generating_model(dirichlet_strength=dirichlet_strength)
        traits.write_logo(data_generating_model, os.path.join(directory,'data_generating.png'))
        traits.write_logo(data_generating_model, os.path.join(directory,'data_generating_rev_comp.png'), rev_comp=True)
        logger.info('Using background mosaic model of order %d and %d states' % (self.background_order, self.num_background_mosaics))
        logger.info('Implanting gapped pssm of length %d (without gap)' % self.K)
        logger.info('Generating %d sequences of length %d' %(N, T))
        seqs, true_binding_sites = test_data_from_model(num_seqs=N, length=T, model=data_generating_model, traits=traits)
        logger.info('%d bases are in binding sites' % (sum(tbs.sum() for tbs in true_binding_sites)))
        logger.info('Will replace %d%% of the sequences with \'N\'s' % int(100*p_unknown))
        place_random_unknowns(seqs, data_generating_model.converter.order_0_size, p_unknown=p_unknown)
        logger.info('Left with %d known bases' % hmm.pssm.num_known_bases(seqs))

        save_synthetic_data(directory, T, seed, seqs, true_binding_sites)

        return traits, seqs, true_binding_sites

if '__main__' == __name__:
    import hmm.pssm, numpy, os, os.path, logging, random, sys, math, traceback, cPickle
    from background_models import global_background_model_cache
    from sequences import *

    root_dir = 'synthetic'
    if not os.access(root_dir, os.R_OK):
        os.makedirs(root_dir)
    logging.shutdown()
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().addHandler(logging.FileHandler(os.path.join(root_dir, 'log.txt'), 'w'))

    try:
        import cookbook
        cookbook.make_current_process_nice()
    except:
        print 'Could not set process priority'

    def background_model(order, N):
        return hmm.as_state_model(global_background_model_cache().load((order, N, 'T00594')))

    generator = SyntheticDataGenerator(
      K = 10,
      background_order = 0,
      num_background_mosaics = 1,
      background_model_creator = background_model
    )
    N = 200
    lengths = [
            20,
            #50,
            #100,
            #500,
            #1000
    ]
    num_seeds = 20

    for T in lengths:
        for seed in xrange(1,num_seeds+1):
            directory = os.path.join(root_dir, '%d_%d' % (T, seed))
            if not os.access(directory, os.R_OK):
                os.makedirs(directory)

            logger = logging.getLogger('find_motif.%d.%d' % (T, seed))
            logger.addHandler(logging.FileHandler(os.path.join(directory,'log.txt'), 'w'))
            logger.setLevel(logging.INFO)

            logger.info('**************** T=%d, seed=%d *****************' % (T,seed))
            seqs, true_binding_sites = generator.generate(N, T, seed, directory, logger)
