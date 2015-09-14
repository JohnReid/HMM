#
# Copyright John Reid 2008
#

"""
Code to test simple PSSMs.
"""



if '__main__' == __name__:
    import hmm, numpy as N, numpy.random as R, time, logging

    logging.basicConfig(level=logging.INFO)

    nucleo_dists = N.array(
      [
        [.85, .05, .05, .05],
        [.05, .85, .05, .05],
        [.85, .05, .05, .05],
        [.05, .85, .05, .05],
        [.85, .05, .05, .05],
        [.05, .85, .05, .05],
        [.85, .05, .05, .05],
        [.05, .85, .05, .05],
        [.85, .05, .05, .05],
        [.05, .85, .05, .05],
        [.85, .05, .05, .05],
        [.05, .85, .05, .05],
      ]
    )
    logging.info('Calculating log scores')
    pssm_scores = hmm.calculate_log_scores(nucleo_dists)
    logging.info('Calculating complementary scores')
    comp_scores = hmm.calculate_complementary_scores(pssm_scores)
    logging.info('Preprocessing sequence')
    seq = hmm.preprocess_sequence(N.array([0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,4]))

    logging.info('Scoring sequence')
    logging.info(hmm.score_sequence(pssm_scores, seq))
    logging.info('Scoring sequence with complementary scores')
    logging.info(hmm.score_sequence(comp_scores, seq))

    logging.info(hmm.max_score_in_sequence(pssm_scores, seq))
    logging.info(hmm.max_score_in_sequence(comp_scores, seq))

    long_seqs = [R.random_integers(0,4,size=10000) for i in xrange(100)]
    long_seqs_preprocessed = hmm.preprocess_sequences(long_seqs)

    logging.info('Starting to time max scores on long sequences.')
    start = time.time()
    max_scores = hmm.max_scores_in_sequences(pssm_scores, long_seqs_preprocessed)
    max_comp_scores = hmm.max_scores_in_sequences(comp_scores, long_seqs_preprocessed)
    logging.info('Max scores (and complementary scores) over sequences took %f secs' % (time.time()-start))

    long_seq = R.random_integers(0,4,size=1000000)
    long_seq_preprocessed = hmm.preprocess_sequence(long_seq)

    logging.info('Starting to time max scores on long sequence.')
    start = time.time()
    logging.info(hmm.max_score_in_sequence(pssm_scores, long_seq_preprocessed))
    logging.info(hmm.max_score_in_sequence(comp_scores, long_seq_preprocessed))
    logging.info('Max scores (and complementary scores) for %d bases took %f secs' % (len(long_seq), time.time()-start))
