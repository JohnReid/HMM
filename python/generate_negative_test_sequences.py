#
# Copyright John Reid 2008
#

"""
Code to load a HMM from its pickle file, generate a number of sequences from it and save to a fasta file
"""

import logging, cPickle
from optparse import OptionParser
from sequences import numpy_to_seq



logging.basicConfig(level=logging.DEBUG)

#
# Build, parse and log the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-n",
  "--num-seqs",
  type='int',
  dest="num_seqs",
  help="The number of sequences to generate."
)
option_parser.add_option(
  "-l",
  "--seq_length",
  type='int',
  dest="seq_length",
  help="The length of the sequences to generate."
)
option_parser.add_option(
  "-m",
  "--model-file",
  dest="model_file",
  help="The filename of the pickled HMM."
)
option_parser.add_option(
  "-o",
  "--output",
  dest="output",
  help="The output filename."
)
option_parser.add_option(
  "--fasta-file-line-length",
  dest="fasta_file_line_length",
  type='int',
  default=80,
  help="The line length of the fasta file."
)
# parse them
options, args = option_parser.parse_args()
# log them
for option in option_parser.option_list:
    if option.dest:
        logging.info('%s: %s (%s)', option.dest, str(getattr(options, option.dest)), option.help)


# load HMM
logging.info('Loading HMM from %s', options.model_file)
model = cPickle.load(open(options.model_file))
hmm_info = 'HMM: states %d; markov order %d; filename %s' % (model.N, model.converter.n, options.model_file)

# create sequences
logging.info('Creating %d sequences of length %d', options.num_seqs, options.seq_length)
samples = [model.sample(options.seq_length) for i in xrange(options.num_seqs)]
sequences = map(numpy_to_seq, (seq for states, seq in samples))

# save sequences
chunk_length = options.fasta_file_line_length
logging.info('Saving sequences to %s', options.output)
output = open(options.output, 'w')
for i, seq in enumerate(sequences):
    print >>output, '> Sequence %d from %s' % (i, hmm_info)
    for chunk_start in xrange(0, len(seq), chunk_length):
        print >>output, seq[chunk_start:chunk_start+chunk_length]
