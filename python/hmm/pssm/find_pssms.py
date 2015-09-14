#
# Copyright John Reid 2007,2008
#

"""
Older code to learn HMM PSSMs from sequences.
"""


import hmm, numpy, os, pickle, logging, time

def find_sites(order_n_seq, model, background_states):
    """
    Breaks a sequence up into binding sites. Each site is a maximal subsequence
    that is not background

    Yields (site, states, start) where the site is the subsequence, states are the
    viterbi states and start is the index in the original sequence the site is at
    """
    states  = model.viterbi(order_n_seq)[1]
    site = []
    site_states = []
    for i, (state, x) in enumerate(zip(states, order_n_seq)):
        if state in background_states:
            if len(site):
                yield site, site_states, i
                site = []
                site_states = []
        else:
            site.append(model.converter.make_order_0(x))
            site_states.append(state)
    if len(site):
        yield site, site_states, i

def find_non_rev_comp_sites(order_n_seq, model, background_states, rev_comp_states):
    """
    Same as find_sites but flips reverse_complement sites

    Yields (site, states, start, is_rev_comp)
    """
    for site, states, start in find_sites(order_n_seq, model, background_states):
        if states[0] not in rev_comp_states:
            yield site, states, start, False
        else:
            site.reverse()
            states.reverse()
            site = [ hmm.pssm.complement(b) for b in site ]
            states = [ rev_comp_states[state] for state in states ]
            yield site, states, start, True



class PssmLearner(object):
    '''
    Learns pssms from given sequences
    '''

    def __init__(
            self,
            seqs,
            pssm_traits,
            num_models=1,
            logger=None
    ):
        '''
        Initialise a PssmLearner object

        seqs: the sequences that we try to find the binding sites in.
        pssm_traits: describes how to build a HMM representation of a binding site
        num_models: # of models to learn concurrently
        '''
        self.order_0_seqs = seqs
        self.pssm_traits = pssm_traits
        self.converter = hmm.MarkovOrderConverter(4, pssm_traits.order)
        self.num_models = num_models
        self.LL_history = []

        #convert the sequences to the correct order for the model
        self.order_n_seqs = [ self.converter.to_order_n(seq) for seq in self.order_0_seqs ]
        for seq in self.order_n_seqs:
            seq.setflags(write=True)

        if None == logger:
            self.logger = logging.getLogger()
        else:
            self.logger = logger

    def bw_callback(self, LL):
        self.LL_history.append(LL)
        #self.logger.info('LL: %.3f' % LL)

    def train_model(self, model, use_baum_welch=True):
        'Returns LL, num_iterations'
        # train them
        start = time.clock()
        if use_baum_welch:
            LL, iterations = model.baum_welch(
                    self.order_n_seqs,
                    prior = self.pssm_traits.prior(model),
                    tolerance = self.default_tolerance,
                    callback = self.bw_callback,
                    #max_iterations = 1
            )
        else:
            LL, iterations = model.viterbi_training(
                    self.order_n_seqs,
                    prior = self.pssm_traits.prior(model),
                    tolerance = self.default_tolerance,
                    callback = self.bw_callback
            )
        elapsed = time.clock() - start
        self.logger.info('Trained model: LL=%.3f, iterations=%d, secs/iteration=%.3f' % (LL, iterations, elapsed/iterations))
        return LL, iterations


    def learn_model(self):
        """
        Creates several models, trains them using Baum-Welch and returns the best
        """
        self.logger.info('Learning initial model')

        # work out the parameters for baum-welch and run it
        self.known_bases = sum(self.converter.num_known_bases_order_n(seq) for seq in self.order_n_seqs)
        self.default_tolerance = 1e-6 * self.known_bases

        # create models
        models = [
                hmm.model_states_2_model(self.pssm_traits.new_model())
                for i in xrange(self.num_models)
        ]

        LLs = [ (model,self.train_model(model)[0]) for model in models ]

        # which had the best log likelihood?
        best = max(LLs, key=lambda x:x[1])
        # print best, LLs
        return best


    def improve_model(self, model, LL):
        """
        Generates alternatives as defined by the model's traits and trains them.

        Returns best model, best LL
        """
        best_model, best_LL = model, LL
        for alternative_gen in self.pssm_traits.alternative_generators(model):
            alt_LL = LL
            while True:
                self.logger.info('Trying to improve model')
                alt_model = alternative_gen(model)
                new_LL = self.train_model(alt_model)[0]
                if new_LL < alt_LL + self.default_tolerance: # stop if not improving LL
                    break
                alt_LL = new_LL
            if alt_LL > best_LL:
                best_model, best_LL = alt_model, alt_LL
        return best_model, best_LL


    def analyse_sites(self, model, blank_sites=False):
        """
        Takes a learnt model and finds the sites that make it up.

        If blank_sites is True, then overwrites the sites with missing data
        """

        # analyse the sites that we found
        background_states = self.pssm_traits.background_states
        rev_comps = self.pssm_traits.reverse_complements

        sites = []
        for order_n_seq, order_0_seq in zip(self.order_n_seqs, self.order_0_seqs):
            seq_sites = []
            for site, states, start, is_rev_comp in find_non_rev_comp_sites(order_n_seq, model, background_states, rev_comps):
                seq_sites.append((site, states, start, is_rev_comp)) # copy the site
                if blank_sites:
                    order_n_seq[start:start+len(site)] = model.M # blank the site out in the order n sequence
            sites.append(seq_sites)

        total_sites = sum(len(seq_sites) for seq_sites in sites)
        self.logger.info('Found %d sites' % total_sites)

        return total_sites, sites


    def learn_pssm(self):
        """
        Finds a pssm in the sequences

        Creates several models and trains them using Baum-Welch.
        Selects the best model, finds the sites and blanks them out
        """
        self.LL_history = []
        model, LL = self.learn_model()
        learnt_LL_history = self.LL_history
        self.LL_history = []

        if numpy.isfinite(LL):
            self.LL_history = []
            model, LL = self.improve_model(model, LL)
            self.logger.info('Improved LL to: %f' % LL)
            improved_LL_history = self.LL_history
            self.LL_history = []
        else:
            improved_LL_history = []

        return model, LL, learnt_LL_history, improved_LL_history

class MultiplePssmLearner(object):
    def __init__(self, pssm_learner, out_dir):
        self.pssm_learner = pssm_learner
        self.out_dir = out_dir
        if not os.access(self.out_dir, os.X_OK):
            os.makedirs(self.out_dir)

    def learn_pssms(self, max_to_find = 3, max_attempts = 10):
        '''
        Learns pssms, if interesting blanks the sites out and carries on. If not
        interesting ignores and tries again

        Returns (found, discarded) where found is a list of (model,sites,counts)
        and discarded is a list of models
        '''
        found = []
        discarded = []
        pssm_idx = 0
        while len(found) < max_to_find and len(found)+len(discarded) < max_attempts:

            # look for a pssm
            self.pssm_learner.logger.info('************** Looking for pssm %d' % pssm_idx)
            model, LL, learnt_LLs, improved_LLs = self.pssm_learner.learn_pssm()

            # is it an interesting pssm
            is_interesting = self.pssm_learner.pssm_traits.is_model_interesting(model)
            self.pssm_learner.logger.info('pssm %d %s interesting' % (pssm_idx, is_interesting and 'is' or 'is not'))
            total_sites, sites = self.pssm_learner.analyse_sites(model, blank_sites=is_interesting)
            if total_sites and is_interesting:
                found.append((model, sites))
            else:
                discarded.append(model)

            if 0 == total_sites:
                self.pssm_learner.logger.info("Didn't find any sites")
                continue

            # save the sites and the model
            basename = os.path.join(self.out_dir, 'pssm_%d' % pssm_idx)
            self.save_sites(sites, open(basename+'_sites.txt', 'w'))
            self.save_model(model, learnt_LLs, improved_LLs, basename)
            pssm_idx += 1
        return found, discarded

    def save_model(self, model, learnt_LLs, improved_LLs, basename):
        import pylab
        pickle.dump( model, open(basename+'_model.pickle', 'w') )
        self.pssm_learner.pssm_traits.log_info(model, self.pssm_learner.logger)
        self.pssm_learner.pssm_traits.write_logo(model, basename+'_logo.png')
        pickle.dump(learnt_LLs, open(basename+'_LLs_learnt.pickle', 'w'))
        pickle.dump(improved_LLs, open(basename+'_LLs_improved.pickle', 'w'))
        pylab.figure()
        pylab.plot(learnt_LLs + improved_LLs)
        pylab.savefig(basename+'_LLs.png')
        pylab.close()

    def save_sites(self, sites, f):
        f.write('# site sequence;sequence;start;is reverse complement; state sequence\n')
        for seq_idx, seq_sites in enumerate(sites):
            for site, states, start, is_rev_comp in seq_sites:
                if is_rev_comp:
                    site = hmm.pssm.rev_comp(site)
                f.write('%s;%d;%d;%d;%s\n' % (hmm.pssm.numpy_to_seq(site), seq_idx, start, is_rev_comp, states))





if '__main__' == __name__:
    from hmm.pssm import PssmTraits, create_background_model, seq_to_numpy, random_sequence, information_content
    from random import random

    site = seq_to_numpy('aaactcaa')
    K = len(site)
    rev_comp_site = seq_to_numpy('ttgagttt')
    num_seqs = 60
    seq_length = K + 30
    start = 20
    def gen_sequence():
        seq = random_sequence(seq_length)
        if random() > .5:
            seq[start:start+K] = site
        else:
            seq[start:start+K] = rev_comp_site
        return seq

    p_binding_site = .01
    order = 1
    num_background_states = 2

    # try to learn pssm for these sequences
    seqs = [ gen_sequence() for i in xrange(num_seqs) ]

    def bw_callback(LL):
        print LL

    traits = PssmTraits(K, p_binding_site, order, num_background_states, create_background_model)
    pssm_learner = PssmLearner(seqs, traits)
    #pssm_learner.bw_callback = bw_callback
    multiple_learner = MultiplePssmLearner(pssm_learner, 'test_multiple')
    found, discarded = multiple_learner.learn_pssms(max_to_find = 2, max_attempts = 2)
    print 'Found %d models, discarded %d' % (len(found), len(discarded))
    for model, sites, counts in found:
        print '# sites: %d' % len(sites)
        dist = traits.pssm_dist(model)
        IC = information_content(dist)
        print 'IC: %f' % IC
        #image.show()
