#
# Copyright John Reid 2008
#

"""
Code to implement a MEME-like single gap motif finding algorithm.
"""


import hmm, hmm.pssm, logging, time, math
from hmm.pssm import single_gap
from sequences import *



class SingleGapAlgorithm(object):
    """
    MEME-like single gap motif finding algorithm.
    """

    def __init__(
      self,
      K=12,
      init_K_mer_length=10,
      max_K_mers_to_evaluate=150,
      expected_sites_per_sequence=1.
    ):
        self.K = K
        "Length of the motifs (excluding the gap)."

        self.init_K_mer_length = init_K_mer_length
        "Length of K-mers to look for in sequences to use as starting points."

        self.max_K_mers_to_evaluate = max_K_mers_to_evaluate
        "Maximum number of K-mers to evaluate as starting points."

        # create something to build the gapped pssms
        self.builder = single_gap.SingleGappedPssmBuilder(K=K, gap_position=K/2, markov_order=0, M=4)
        "Used to initialise motifs."

        self.initialisation_pseudo_count = .4
        "How much to add to the initialisation K-mers."

        self.expected_sites_per_sequence = expected_sites_per_sequence
        "Expected number of sites per sequence."


    def model_for_initialisation_K_mer(self, K_mer, p_binding_site):
        """
        Create a model initialised by this K-mer.
        """
        emission_distributions = numpy.ones((self.K, 4)) * self.initialisation_pseudo_count
        for k, base in enumerate(K_mer):
            if 4 == base:
                emission_distributions[k] += .25
            else:
                emission_distributions[k, base] += 1.
        gap_emissions = numpy.ones((4,)) / 4.
        pssm, in_states, out_states = self.builder.create(
          p_gap=.5,
          non_gap_emissions=emission_distributions,
          gap_emissions=gap_emissions
        )
        model = hmm.as_model(
          single_gap.add_to_simple_background_model(
            model=pssm,
            in_states=in_states,
            out_states=out_states,
            p_binding_site=p_binding_site
          )
        )
        return model



    def evaluate_initialisation_K_mer(self, K_mer, sequences):
        """
        Find out how good a starting point this K_mer would be.
        """
        pssms = self.pssms_for_K_mer(K_mer)
        max_scores_per_pssm = numpy.array(
          [
            hmm.max_scores_in_sequences(pssm, sequences)
            for pssm
            in pssms
          ]
        )
        best_scores = max_scores_per_pssm.max(axis=0)
        null_score = (self.K+1)*math.log(.25)
        best_scores[best_scores<null_score] = null_score
        return best_scores.sum()


    def pssms_for_K_mer(self, K_mer):
        """
        Return 4 pssms for gapped/ungapped and complementary/uncomplementary strands.
        @return: pssm, comp_pssm, gap_pssm, gap_comp_pssm
        """
        dist = self.nucleo_dist_from_K_mer(K_mer, include_gap=False)
        gap_dist = self.nucleo_dist_from_K_mer(K_mer, include_gap=True)
        pssm = hmm.calculate_log_scores(dist)
        comp_pssm = hmm.calculate_complementary_scores(pssm)
        gap_pssm = hmm.calculate_log_scores(gap_dist)
        gap_comp_pssm = hmm.calculate_complementary_scores(gap_pssm)
        return pssm, comp_pssm, gap_pssm, gap_comp_pssm


    def nucleo_dist_from_K_mer(self, K_mer, include_gap=False):
        # initialise our distribution to 'N's
        nucleo_dist = numpy.ones((self.K+1, 4)) * self.initialisation_pseudo_count

        # calculate where to put the K-mer in the distribution
        mer_len = len(K_mer)
        start = self.K/2 - mer_len/2
        if start < 0:
            raise RuntimeError('K_mer too large to build PSSM from.')

        # add K_mer to distribution
        for i, c in enumerate(K_mer):
            if include_gap and i >= mer_len/2: # skip the gap if necessary
                pos = start+i+1
            else:
                pos = start+i
            if 4 != c:
                nucleo_dist[pos,c] += 1
            else:
                nucleo_dist[pos] += .25

        # normalise
        for dist in nucleo_dist:
            dist /= dist.sum()

        return nucleo_dist


    def yield_evaluations(self, nmer_counts, preprocessed_sequences):
        from heapq import heapify, heappop

        start = time.time()
        counts = list((-count, i, K_mer) for i, (K_mer, count) in enumerate(nmer_counts.counts()))
        heapify(counts)
        logging.info('Took %f seconds to heapify', time.time()-start)
        #import IPython; IPython.Debugger.Pdb().set_trace()

        left_to_evaluate = self.max_K_mers_to_evaluate
        logging.info('Evaluating up to %d %d-mers', left_to_evaluate, self.init_K_mer_length)
        while counts and left_to_evaluate > 0:
            count, i, K_mer = heappop(counts)
            evaluation = self.evaluate_initialisation_K_mer(K_mer, preprocessed_sequences)
            logging.info('K-mer evaluation: %s - %f', str(K_mer), evaluation)
            yield K_mer, evaluation
            left_to_evaluate -= 1



    def __call__(self, sequences):
        """
        Run the motif finding algorithm.
        """

        preprocessed_sequences = hmm.preprocess_sequences(sequences)

        # how big are the sequences
        num_bases = sum(len(s) for s in sequences)

        # find all K-mers collapsed with their reverse complements
        logging.info('Finding all %d-mers in sequences', self.init_K_mer_length)
        start = time.time()
        nmer_counts = hmm.ReverseComplementCollapsingCounter(self.init_K_mer_length)
        hmm.count_mers(
          sequences,
          n=self.init_K_mer_length,
          callback=nmer_counts
        )
        logging.info('Took %f seconds to find %d-mers', time.time()-start, self.init_K_mer_length)

        p_binding_site = (self.expected_sites_per_sequence * len(sequences)) / num_bases
        logging.info('Found %d %d-mers', nmer_counts.num_counts(), self.init_K_mer_length)
        start = time.time()
        best_starting_point = max(
          self.yield_evaluations(nmer_counts, preprocessed_sequences),
          key=lambda x: x[1]
        )
        logging.info('Evaluation took %f seconds', time.time() - start)
        logging.info('Best starting point: %s: %f' % best_starting_point)

        model = self.model_for_initialisation_K_mer(best_starting_point[0], p_binding_site)
        logging.info('Running Baum-Welch')
        start = time.time()
        LL, num_iterations = model.baum_welch(preprocessed_sequences)
        logging.info('Baum-Welch took %f seconds', time.time()-start)
        logging.info('Achieved LL: %f in %d iterations', LL, num_iterations)

        return model





if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    def synthetic():
        identifier = 'synthetic-sequences-K10-g0.50-N200-L200-seed4-1'
        sequences = [
          hmm.pssm.seq_to_numpy(s)
          for s
          in convert_seqs('synthetic-2/%s.fa' % identifier)
        ]
        return identifier, sequences

    def fragment(identifier = 'T00594'):
        sequences = seqs_for_fragment(identifier)
        return identifier, sequences

    identifier, sequences = synthetic()
    identifier, sequences = fragment()

    algorithm = SingleGapAlgorithm()

    model = algorithm(sequences)

    emissions, gap_probs = algorithm.builder.get_emissions_and_gap_probabilities(model, offset=1)
    import hmm.pssm.logo as logo
    image = logo.pssm_as_image(emissions, transparencies=gap_probs)
    image.save("single-gap-results/%s.png" % identifier, "PNG")
