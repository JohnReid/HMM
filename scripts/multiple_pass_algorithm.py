#
# Copyright John Reid 2008
#

class MultiplePassAlgorithm(object):
    """
    Takes an algorithm that finds binding sites and apply it multiple times.
    """

    def __init__(self, single_pass_algorithm, num_passes, remove_sites_fn):
        """
        @arg single_pass_algorithm: Algorithm that finds binding sites.
        @arg num_passes: Number of passes to make of the algorithm.
        @arg remove_sites_fn: Function that removes sites from sequences given results from the single pass algorithm.
        """

        self.single_pass_algorithm = single_pass_algorithm
        "Algorithm that finds binding sites."

        self.num_passes = num_passes
        "Number of passes to make of the algorithm."

        self.remove_sites_fn = remove_sites_fn
        "Function that removes sites from sequences given results from the single pass algorithm."



    def __call__(self, sequences):
        "Run the single pass algorithm multiple times."

        results = []

        # for each pass
        for i in xrange(num_passes):

            # run the single pass algorithm
            result = self.single_pass_algorithm(sequences)
            results.append(result)

            # blank the bases in sequence that were found
            self.remove_sites_fn(sequences, result)

        return results
