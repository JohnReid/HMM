#
# Copyright John Reid 2007
#


"""
Implementation of PSSM finding algorithm - supercedes that in test_synthetic.py
"""


import hmm.pssm, numpy, os, os.path, logging, random, sys, math, traceback

from count_mers import *
from background_models import background_model
from generate_synthetic import *
from process_util import *
from optimisation import *


class PssmFindingAlgorithm(object):
    def __init__(
            self,
            traits,
            directory='pssms',
            logger=None,
            n_mer_strength=0.0,
    ):
        self.traits = traits
        self.directory = directory
        self.logger = logger and logger or logging.getLogger()
        self.n_mer_strength = n_mer_strength

        #make the directory for the output
        if not os.access(directory, os.R_OK):
            os.makedirs(directory)
        if not os.access(self.file('detail'), os.R_OK):
            os.makedirs(self.file('detail'))

    def count_n_mers(self, seqs, n=5):
        "Find the most common n-mer in the sequences and make the model traits use it to initialise models"
        n_mer, n_mer_count, LL = most_significan_n_mer(seqs.as_order[0], n)
        #for mer, count in top_k_n_mers(seqs.as_order[0], 5, 6):
        #       self.logger.info('n-mer count: %d: %s' % (count, str(mer)))
        self.logger.info('Using n-mer %s (%d) to initialise motif' % (str(n_mer), n_mer_count))
        self.traits.emission_dists = self.traits.emission_dist_including_n_mer(n_mer, strength=self.n_mer_strength)

    def new_model(self, seqs):
        self.count_n_mers(seqs)
        return hmm.model_states_2_model(self.traits.new_model())

    def LL_callback(self, model):
        LLs = []
        def callback(LL):
            LLs.append(LL)
            iteration = len(LLs)
            self.traits.write_logo(model, self.file('detail/logo_%d.png' % iteration))
            if len(LLs) > 1:
                pass ; self.logger.info('LL diff: %d, %f' % (iteration, LLs[-1]-LLs[-2]))
            stop = self.stopping_condition(LLs, model)
            #print 'stop=', stop
            return not stop
        return LLs, callback

    def stopping_condition(self, LLs, model):
        return False

    def file(self, name):
        return os.path.join(self.directory, name)

    def model_from_motif(self, motif):
        return motif[0]

    def make_prediction(self, motif, seqs):
        model = self.model_from_motif(motif)
        return [
                self.traits.binding_sites_from_states(model.viterbi(seq)[1])
                for seq
                in seqs.as_order[model.converter.n]
        ]

    def __call__(self, seqs):
        # details about the sequences
        num_known_bases = seqs.num_known_bases()
        self.logger.info('Finding pssm in %d sequences (totalling %d known bases)' % (len(seqs.seqs), num_known_bases))

        model, LLs = self.run(seqs)

        self.graph_LLs(LLs)
        self.traits.write_logo(model, self.file('logo.png'))
        return (model, LLs)



class PssmFindingHmmTrainingAlgorithm(PssmFindingAlgorithm):

    def __init__(
            self,
            traits,
            directory='pssms',
            max_iterations=0,
            tolerance=1.0,
            logger=None,
            use_viterbi=False,
            n_mer_strength=0.0,
            prior=None
    ):
        PssmFindingAlgorithm.__init__(self, traits, directory, logger, n_mer_strength=n_mer_strength)
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.prior = prior
        if use_viterbi:
            self.training_algorithm = hmm.Model.viterbi_training
        else:
            self.training_algorithm = hmm.Model.baum_welch

    def stopping_condition(self, LLs, model):
        "Returns true if we should stop training"
        #we need at least 6 LLs
        if 6 > len(LLs):
            #self.logger.info('Not enough LLs')
            return False

        #have we reached our tolerance?
        diffs = self.LL_diffs(LLs)
        if diffs[-1] > self.tolerance:
            #self.logger.info('Not reached tolerance')
            return False

        #make sure last 5 diffs are all decreasing
        for i in xrange(1,5):
            if diffs[-i] > diffs[-i-1]:
                #self.logger.info('diffs not decreasing - %f > %f' % (diffs[-i], diffs[-i-1]))
                return False

        return True


    def run(self, seqs):
        model = self.new_model(seqs)

        LLs, callback = self.LL_callback(model)
        self.logger.info('Running %s with tolerance of %f' % (str(self.training_algorithm), self.tolerance))
        if self.max_iterations > 0:
            logger.info('Max iterations: %d' % self.max_iterations)
        if None != self.prior:
            prior = self.prior
        else:
            prior = self.traits.prior(model)
        #import pdb; pdb.set_trace()
        LL, iterations = self.training_algorithm(
                model,
                seqs.as_order[model.converter.n],
                prior = prior,
                tolerance = -1e-99,
                callback = callback,
                max_iterations = self.max_iterations
        )
        self.logger.info('LL: %f' % LL)
        self.logger.info('# iterations: %d' % len(LLs))
        return model, LLs

    @staticmethod
    def LL_diffs(LLs):
        return [ (LLs[i+1]-LLs[i]) for i in xrange(len(LLs)-1) ]

    def graph_LLs(self, LLs):
        import pylab
        pylab.clf()
        pylab.plot(LLs)
        pylab.savefig(self.file('LLs.png'))
        pylab.clf()
        pylab.plot(self.LL_diffs(LLs))
        pylab.savefig(self.file('LL_diffs.png'))


class OptimisingPssmFinder(PssmFindingAlgorithm):
    def __init__(
            self,
            traits,
            directory='optimising',
            logger=None,
            #optimizer=scipy.optimize.fmin,
            #optimizer=scipy.optimize.fmin_powell,
            optimizer=scipy.optimize.anneal,
            n_mer_strength=0.0
    ):
        PssmFindingAlgorithm.__init__(self, traits, directory, logger, n_mer_strength=n_mer_strength)
        self.optimizer = optimizer

    def build_param_indices(self, model):
        'Return a list of those parameter indices that we want to optimise'
        param_idx_map = set()
        def add_param(idx):
            if(idx):
                param_idx_map.add(idx.idx)
        for i in xrange(model.N): #for each state, i
            # add the interesting emission parameters
            if i not in self.traits.background_states:
                for k in xrange(model.M): # for each emission, k
                    add_param(model.get_emission_parameterisation(i, k))
            # add the interesting transition parameters
            to_states = []
            for j in xrange(model.N): # for each state, j
                if i not in self.traits.background_states or j not in self.traits.background_states:
                    idx = model.get_transition_parameterisation(i, j)
                    if idx:
                        to_states.append(idx)
            if len(to_states) > 1:
                for idx in to_states:
                    add_param(idx)
        return [idx for idx in param_idx_map]

    def new_model(self, seqs):
        return hmm.model_states_2_model(self.traits.new_model())

    def __call__(
      self,
      seqs
    ):
        model = self.new_model(seqs)
        LLs, callback = self.LL_callback(model)
        optimisable = OptimisableModel(model, self.build_param_indices(model), seqs.as_order[model.converter.n], callback)
        outputs = self.optimizer(
                optimisable,
                optimisable.get_model_params(),
                full_output=True
        )
        best_params = outputs[0]
        negative_LL = outputs[1]
        optimisable.set_model_params(best_params)
        return model, -negative_LL

def find_motifs(seqs, algorithm, is_interesting, num_to_find, max_tries):
    """
    Find several motifs in the sequences using the given algorithm
    """
    motifs = []
    while num_to_find > len(motifs) and max_tries > 0:
        motif = algorithm.find_motif(seqs)
        if is_interesting(motif, seqs):
            motifs.append(motif)
            blank_sites_for_motif(seqs, motif)
        max_tries -= 1
    return motifs


class GetTestDataSites(object):
    "Gets the sites from the sequences and given positions"
    def __init__(self):
        self.site = []
        self.sites = []

    def add_sites_from(self, seq, binding_sites):
        for base, indicator in zip(seq, binding_sites):
            if not indicator:
                self._site_completed()
            else:
                self.site.append(base)
        self._site_completed()

    def _site_completed(self):
        if len(self.site):
            self.sites.append(numpy.array(self.site))
            self.site = []

def calculate_sites(seqs, binding_sites):
    sites = GetTestDataSites()
    for seq, binding_sites in zip(seqs, binding_sites):
        sites.add_sites_from(seq, binding_sites)
    return sites.sites

def calculate_roc(true_binding_sites, predicted_binding_sites):
    from infpy.roc import RocCalculator, update_roc
    roc = RocCalculator()
    for predicted, truth in zip(predicted_binding_sites, true_binding_sites):
        update_roc( roc, zip(truth,predicted) )
    return roc

class MarkovOrderDict(dict):
    def __init__(self, seqs, alphabet_size=4):
        self.seqs = seqs
        self.alphabet_size = alphabet_size
    def __missing__(self, n):
        v = [
                hmm.MarkovOrderConverter(self.alphabet_size, n).to_order_n(s)
                for s in self.seqs
        ]
        self[n] = v
        return v

class TestSequences(object):
    def __init__(self, seqs):
        self.as_order = MarkovOrderDict(seqs)
        self.seqs = seqs
    def num_known_bases(self):
        return hmm.pssm.num_known_bases(self.as_order[0])

def format_string(i):
    return _format_strings[i%len(_format_strings)]

_colors = [
        'b',
        'g',
        'r',
        'c',
        'm',
        'y',
        'k',
]
_markers = [
        's',
        '<',
        'o',
        '>',
        'v',
        'd',
        '^',
        'p',
        'h',
        '8',
]
def plot_rocs(all_rocs):
    import pylab
    scatter = pylab.figure()
    for i, (T, rocs) in enumerate(all_rocs.iteritems()):
        sensitivities = numpy.array([roc.sensitivity() for roc in rocs])
        specificities = numpy.array([roc.specificity() for roc in rocs])
        pylab.scatter(
                sensitivities,
                specificities,
                marker=_markers[i],
                color=_colors[i],
                label='%d'%T
        )
    pylab.legend()
    pylab.axis(xmax=1.0, ymax=1.0, figure=scatter)

if '__main__' == __name__:

    root_dir = 'pssms'
    logger = setup_process(root_dir)

    N = 200
    lengths = [
            #20,
            #50,
            100,
            #500,
            #1000
    ]
    num_seeds = 20
    K = 10
    background_order = 2
    num_background_mosaics = 4

    class DictOfLists(dict):
        def __missing__(self, k):
            self[k] = []
            return self[k]
    all_rocs = DictOfLists()

    for T in lengths:
        for seed in xrange(1,num_seeds+1):
            directory = os.path.join(root_dir, '%d_%d' % (T, seed))
            if not os.access(directory, os.R_OK):
                os.makedirs(directory)

            logger = logging.getLogger('find_motif.%d.%d' % (T, seed))
            logger.addHandler(logging.FileHandler(os.path.join(directory,'log.txt'), 'w'))
            logger.setLevel(logging.INFO)
            logger.info('**************** T=%d, seed=%d *****************' % (T,seed))
            logger.info('background model order=%d, # mosaics=%d' % (background_order,num_background_mosaics))

            try:
                hmm.seed_rng( T+seed )
                random.seed( T+seed )

                # traits that determine how our model maps to a PSSM
                p_binding_site = 1.0/T # one binding site per sequence on average
                traits = hmm.pssm.GappedPssmTraits(
                        K=K,
                        p_binding_site=p_binding_site,
                        background_order=background_order,
                        num_background_mosaics=num_background_mosaics,
                        prior_strength = 1.0,
                        background_model_creator=background_model
                )

                # generate some data
                data_generating_model = traits.create_test_data_generating_model()
                traits.write_logo(data_generating_model, os.path.join(directory,'data_generating.png'))
                traits.write_logo(data_generating_model, os.path.join(directory,'data_generating_rev_comp.png'), rev_comp=True)
                seqs, true_binding_sites = test_data_from_model(num_seqs=N, length=T, model=data_generating_model, traits=traits)
                place_random_unknowns(seqs, data_generating_model.converter.order_0_size)
                test_seqs = TestSequences(seqs)

                # find pssm
                #algorithm = PssmFindingHmmTrainingAlgorithm(traits, directory=directory, logger=logger)
                algorithm = OptimisingPssmFinder(traits, directory=directory, logger=logger)
                motif = algorithm(test_seqs)

                # evaluate results
                predicted_binding_sites = algorithm.make_prediction(motif, test_seqs)
                roc = calculate_roc(true_binding_sites, predicted_binding_sites)
                logger.info(str(roc))
                all_rocs[T].append(roc)

            except KeyboardInterrupt:
                logger.error(traceback.format_exc())
                raise

            #except:
            #       logger.error(traceback.format_exc())
            #       sys.exc_clear()

    plot_rocs(all_rocs)
