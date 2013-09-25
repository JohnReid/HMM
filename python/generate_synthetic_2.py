#
# Copyright John Reid 2008
#


"""
Generates synthetic data for tests
"""

import hmm.pssm.single_gap as single_gap
import hmm.pssm.logo as logo
import hmm, hmm.pssm, numpy, numpy.random, os.path
from optparse import OptionParser


parser = OptionParser()
parser.add_option(
  "-K",
  dest="K",
  default=10,
  help="Length of ungapped pssm",
)
parser.add_option(
  "-N",
  dest="N",
  default=20,
  help="Number of sequences",
)
parser.add_option(
  "-L",
  dest="L",
  default=100,
  help="Average length of sequences",
)
parser.add_option(
  "-e",
  "--exp-sites-per-sequence",
  default=2.,
  dest="exp_sites_per_sequence",
  help="Expected number of sites in each sequence",
)
parser.add_option(
  "-p",
  "--p-gap",
  default=.5,
  dest="p_gap",
  help="Probability of a gap in each binding site",
)
parser.add_option(
  "-s",
  "--seed",
  dest="seed",
  default=1,
  help="Seed for random numbers",
)
(options, args) = parser.parse_args()


K = int(options.K)
"Length of ungapped pssm."

p_gap = float(options.p_gap)
"Probability of a gap."

N = int(options.N)
"Number of sequences."

L = int(options.L)
"Average length of sequences."

exp_sites_per_sequence = float(options.exp_sites_per_sequence)
"Expected number of sites in each sequence."

seed = int(options.seed)
"Seed for random numbers."



print 'Going to generate %d sequences of average length %d' % (N, L)
print 'The binding sites are of length %d plus an optional gap with probability %f' % (K, p_gap)
print 'Expect to find %f binding sites per sequence' % exp_sites_per_sequence
print 'Seeding the random number generator with %d' % seed



# seed all the RNGs that we use
hmm.seed_rng(seed)
numpy.random.seed(seed)


# create something to build the gapped pssms
builder = single_gap.SingleGappedPssmBuilder(K=K, gap_position=K/2, markov_order=0, M=4)



# create our emission distributions
dirichlet_prior_strengths = [.01, .1, 1.]
emissions = [
  numpy.array(
    [
      hmm.dirichlet_draw(numpy.ones(builder.M) * strength)
      for k in xrange(builder.K)
    ]
  )
  for strength in dirichlet_prior_strengths
]
gap_emissions = [
  hmm.dirichlet_draw(numpy.ones(builder.M) * strength)
  for strength in dirichlet_prior_strengths
]


# create out single gapped pssms
pssms = [
  builder.create(
    p_gap,
    non_gap,
    gap
  )
  for non_gap, gap
  in zip(emissions, gap_emissions)
]


# create our complete models (by adding a background model)
p_binding_site = exp_sites_per_sequence / L
models = [
  hmm.as_model(single_gap.add_to_simple_background_model(
    model=pssm[0],
    in_states=pssm[1],
    out_states=pssm[2],
    p_binding_site=p_binding_site))
  for pssm
  in pssms
]



# write our logos
# convert to sequences and write fasta
def tag(sample_idx):
    return "K%d-g%.2f-N%d-L%d-seed%d-%d" % (K, p_gap, N, L, seed, sample_idx)
print 'Writing logos'
for i, model in enumerate(models):
    emissions, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)
    image = logo.pssm_as_image(emissions, transparencies=gap_probs)
    image.save("synthetic-2/generating-pssm-logo-%s.png" % tag(i), "PNG")



def yield_sites_in_states(state_sequence, background_states):
    "Yield [start,end) tuples for each site in the state sequence."
    in_site = False
    for i, state in enumerate(state_sequence):
        if not in_site:
            if state not in background_states:
                site_start = i
                in_site = True
        else:
            # we are already in site
            if state in background_states:
                yield site_start, i
                in_site = False
    if in_site:
        yield site_start, i




# sample sequences
def sample_from(model, N, L):
    "Sample N sequences with an expected average length of L."
    from scipy.stats import poisson
    return [
      model.sample(l)
      for l in poisson(L).rvs(size=N)
    ]
print 'Sampling %d sequences with expected average length %d' % (N, L)
samples = [sample_from(model, N, L) for model in models]
print 'Binding bases per model:', [
  sum(sum(s[0] != 0) for s in sample)
  for sample in samples
]
print 'Sites per sample:', [
  sum(len([site for site in yield_sites_in_states(s[0], (0,))]) for s in sample)
  for sample in samples
]


print 'Converting sequences'
sequence_sets = [
  [
    hmm.pssm.numpy_to_seq(s[1])
    for s in sample
  ]
  for sample in samples
]
print 'Writing sequences'
for i, sequences in enumerate(sequence_sets):
    f = open("synthetic-2/synthetic-sequences-%s.fa" % tag(i), 'w')
    for j, s in enumerate(sequences):
        f.write('> sequence %d\n' % j)
        f.write(s)
        f.write('\n')
    f.close()


if False:
    for i, sample in enumerate(samples):
        print '10 largest counts in sequence set: %d' % i
        nmer_counts = hmm.ReverseComplementCollapsingCounter(K)
        hmm.count_mers(
          [sequence[1].astype(int) for sequence in sample],
          n=K,
          callback=nmer_counts
        )
        import heapq
        print "\n".join(
          "%s : %d" % (hmm.pssm.numpy_to_seq(nmer), count)
          for nmer, count
          in heapq.nlargest(10, nmer_counts.counts(), key=lambda count: count[1])
        )
