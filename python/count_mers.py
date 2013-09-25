#
# Copyright John Reid 2007-2008
#

"""
Code to analyse sequences for most common n-mers
"""


import hmm, math, numpy

class _MostCommonNMer(object):
    def __init__(self):
        self.best_mer = None
        self.best_count = 0

    def __call__(self, mer, count):
        if 4 not in mer:
            if self.best_count < count:
                self.best_mer = mer
                self.best_count = count

class _AllNMers(object):
    def __init__(self):
        self.n_mers = []

    def __call__(self, mer, count):
        if 4 not in mer:
            self.n_mers.append((mer, count))

def most_common_n_mer(seqs, n):
    callback = _MostCommonNMer()
    hmm.count_mers(seqs, n, callback)
    return callback.best_mer, callback.best_count

def top_k_n_mers(seqs, n, k):
    all = _AllNMers()
    hmm.count_mers(seqs, n, all)
    import heapq
    return heapq.nlargest(k, all.n_mers, key=lambda x: x[1])

class CountingDict(dict):
    'A dictionary that initialises missing values with 0'
    def __missing__(self, k):
        self[k] = 0
        return self[k]



def collapse_rev_comps(n_mers):
    '''
    Takes a sequence of (n_mer, count) tuples and collapses them so that
    reverse complements counts are aggregated
    '''
    if not len(n_mers):
        return []
    n = len(n_mers[0][0])
    import hmm.pssm
    converter = hmm.MarkovOrderConverter(4, n-1)
    counts = CountingDict()
    for mer, count in n_mers:
        mer.dtype = numpy.int32
        positive_idx = converter.convert_to_order_n_observation(mer)
        rev_comp_mer = hmm.pssm.rev_comp(mer)
        negative_idx = converter.convert_to_order_n_observation(rev_comp_mer)
        idx = min(positive_idx, negative_idx)
        counts[idx] += count
    return [
            (converter.convert_from_order_n_observation(idx), count)
for idx, count
in counts.iteritems()
    ]

class _LogFactorial(dict):
    def __missing__(self, k):
        if 1 == k:
            self[k] = 0
        else:
            self[k] = self[k-1] + math.log(k)
        return self[k]

def calculate_n_mer_significances(seqs, n, background=None):
    '''
    Counts all n-mers in the sequences and assesses the significance of each
    count w.r.t. the background_model

    If the background model is not specified, a uniform distribution over
    the bases is assumed
    '''
    from sys import getrecursionlimit
    all = _AllNMers()
    hmm.count_mers(seqs, n, all)
    collapsed = collapse_rev_comps(all.n_mers)
    log_fact = _LogFactorial()
    total_counts = sum(count for mer, count in collapsed)
    for i in xrange(1,total_counts,getrecursionlimit()/2):
        log_fact[i]
    log_fact_total = log_fact[total_counts]
    if None == background:
        background_LL = n * math.log(.25)
        foreground_LL = math.log(1.0 - math.exp(background_LL))
    result = []
    for mer, count in collapsed:
        if None != background:
            background_LL = background.LL(mer)
            foreground_LL = math.log(1.0 - math.exp(background_LL))
        log_bernoulli = (
                log_fact_total
                - log_fact[count]
                - log_fact[total_counts-count]
                + count * background_LL
                + (total_counts-count) * foreground_LL
        )
        result.append((mer, count, log_bernoulli))
    result.sort(cmp=lambda x,y: cmp(x[2], y[2]))
    return result

def most_significan_n_mer(seqs, n, background=None):
    '''
    Finds the most significant n-mer
    '''
    significant_n_mers = calculate_n_mer_significances(seqs, n, background)
    return significant_n_mers[0]

if '__main__' == __name__:
    import numpy
    seqs = [
            numpy.array([0,1,2,3,4]),
            numpy.array([0,1,2,3,4]),
            numpy.array([0,1,2])
    ]
    print most_common_n_mer(seqs, 3)
    print most_common_n_mer(seqs, 5)

    print
    c=hmm.MarkovOrderConverter(4,2)
    a=numpy.array([0,1,2])
    order_n_obs = c.convert_to_order_n_observation(a)
    a_copy = c.convert_from_order_n_observation(order_n_obs)
    print a
    print a_copy
    print a.all() == a_copy.all()

    print
    seqs = [
            numpy.array([3,2,1]),
            numpy.array([0,1,2]),
            numpy.array([0,1,2]),
            numpy.array([1,2,3]),
    ]
    n = 3
    all = _AllNMers()
    hmm.count_mers(seqs, n, all)
    for mer, count in collapse_rev_comps(all.n_mers):
        print mer, count

    print
    significant_mers = calculate_n_mer_significances(seqs, 3, None)
    print significant_mers[:10]
