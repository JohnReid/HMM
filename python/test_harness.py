#
# Copyright John Reid 2008
#


"""
A test harness to evaluate gapped pssm models.
"""


import os, os.path, sys, logging, numpy, glob
from optparse import OptionParser
from time import strftime
from infpy.roc import RocCalculator, update_roc, get_new_roc_parameter, plot_roc_points
from sequences import convert_seqs, seq_to_numpy
from fragments import all_fragments
from gapped_pssms.parse_gapped_format import parse_models, build_hmm_from_semi_parsed
from itertools import imap


class DictOfRocs(dict):
    def __missing__(self, k):
        self[k] = RocCalculator()
        return self[k]

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed([])

def ensure_dir_exists(dir):
    if os.access(dir, os.X_OK):
        return
    else:
        # make sure parent exists
        parent_dir = os.path.dirname(dir)
        if parent_dir:
            ensure_dir_exists(parent_dir)
            os.mkdir(dir)





def evaluate_model(model, sequence):
    """
    Evaluates the model against the sequence.

    @return: True if there is at least one hit in the sequence
    """
    hmm, traits = model
    #from IPython.Debugger import Pdb; Pdb().set_trace()
    #ipshell() # this call anywhere in your program will start IPython
    LL, states = hmm.viterbi(sequence)
    # we have a hit if we find at least K/2 states in the state sequence that are not in the
    # background
    return sum(state not in traits.background_states for state in states) > traits.K / 2



def evaluate_models(models, sequence):
    """
    Evaluates the models against the sequence.

    @return: True if there is at least one hit for one of the models in the sequence
    """
    if not len(models):
        return False
    else:
        return reduce(
          bool.__or__,
          (evaluate_model(model, sequence) for model in models)
        )



def generate_roc_data(models, positive_sequences, negative_sequences):
    """
    Yields truth, prediction pairs for the model over the positive and negative sequences
    """
    for seq in positive_sequences:
        yield True, evaluate_models(models, seq)
    for seq in negative_sequences:
        yield False, evaluate_models(models, seq)





def get_roc_for_sequences(p_binding_site, positive_sequences, negative_sequences, pssms):
    models = [
      build_hmm_from_semi_parsed(
        parsed,
        p_binding_site=p_binding_site
      )
      for parsed
      in pssms
    ]
    roc = RocCalculator()
    update_roc(roc, generate_roc_data(models, positive_sequences, negative_sequences))
    return roc



def sequences_from_fasta(fasta):
    """Yields sequences from fasta file."""
    import corebio.seq_io.fasta_io
    return imap(
      lambda s: s.strip('nN'),
      imap(
        str,
        corebio.seq_io.fasta_io.iterseq(
          open(fasta, 'r'),
          corebio.seq.dna_alphabet
        )
      )
    )


class NegativeSeqsFromDir(dict):
    """Dict that maps datasets to generators that return negative sequences."""
    def __init__(self, dir):
        self.dir = dir

    def __missing__(self, dataset):
        self[dataset] = sequences_from_fasta(os.path.join(self.dir, 'negative-%s.fa' % dataset))
        return self[dataset]


class NegativeSeqsFromFasta(dict):
    """Dict that maps datasets to a generator that returns negative sequences from a single fasta file."""
    def __init__(self, fasta):
        self.fasta = fasta

    def __missing__(self, dataset):
        self[dataset] = sequences_from_fasta(self.fasta)
        return self[dataset]


class PositiveNegativeSequences(dict):
    """Dict that maps (dataset, cross_fold_idx) pairs to length matched positive/negative sequences"""
    def __init__(self, negative_seq_generators):
        self.negative_seq_generators = negative_seq_generators

    def __missing__(self, key):
        dataset, cross_fold_idx = key
        positive_seqs = [
          seq_to_numpy(s)
          for s
          in convert_seqs(sequence_files['%s-%d-validate' % (dataset, cross_fold_idx)])
        ]
        # load the negative set and match the lengths of the positive sequences

        negative_seqs = [
          seq_to_numpy(neg[:len(pos)])
          for pos, neg
          in zip(
            positive_seqs,
            self.negative_seq_generators[dataset]
          )
        ]
        if len(negative_seqs) != len(positive_seqs):
            raise RuntimeError('Not enough sequences in negative set to match positive sequences. %d positive, %d negative sequences' % (len(positive_seqs), len(negative_seqs)))

        for i, (pos, neg) in enumerate(zip(positive_seqs, negative_seqs)):
            if len(neg) < len(pos):
                raise RuntimeError(
                  'Not enough bases in negative sequence %d to match length of positive sequence: positive sequence has %d bases and negative sequence has %d bases' % (
                    i,
                    pos.shape[0],
                    neg.shape[0]
                  )
                )

        self[key] = (positive_seqs, negative_seqs)
        return self[key]



def get_roc_for(p_binding_site, pssms_for_method, sequence_dict):
    result = RocCalculator()
    result_per_dataset = DictOfRocs()
    for (dataset, cross_fold_idx), pssms in pssms_for_method.iteritems():

        positive_seqs, negative_seqs = sequence_dict[(dataset, cross_fold_idx)]

        logging.info(
          'Analysing %d positive and negative sequences (%d bases) for %s-%d with p(binding site)=%f',
          len(positive_seqs),
          sum(len(pos) for pos in positive_seqs),
          dataset,
          cross_fold_idx,
          p_binding_site
        )
        roc = get_roc_for_sequences(
          p_binding_site,
          positive_seqs,
          negative_seqs,
          pssms
        )
        logging.info('%s-%d: p(binding site): %f\n%s' % (dataset, cross_fold_idx, p_binding_site, str(roc)))
        result += roc
        result_per_dataset[dataset] += roc
    return result, result_per_dataset



logging.basicConfig(level=logging.INFO)
logging.getLogger('').addHandler(logging.FileHandler('test-harness.log', 'w'))

time_str = strftime('%Y-%m-%d_%H-%M-%S')

#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-m",
  "--model-dir",
  dest="model_dir",
  help="A directory containing the gapped pssm model(s)"
)
option_parser.add_option(
  "-n",
  "--negative-sequences",
  default='vm-r1-back-mouse-NCBIM36-reduced',
  dest="negative_sequences",
  help="The name of the dataset with the negative sequences"
)
option_parser.add_option(
  "--negative-directory",
  default=None,
  dest="negative_dir",
  help="A directory containing a set of negative sequences (1 per dataset). These should be named e.g. negative-T00140.fa, etc.."
)
option_parser.add_option(
  "-r",
  "--roc-output-file",
  default='test-harness-roc',
  dest="roc_output_file",
  help="The filename to write the ROC curve to"
)
option_parser.add_option(
  "-t",
  "--roc-title",
  dest="roc_title",
  help="The title of the ROC plot"
)
option_parser.add_option(
  "-o",
  "--roc-statistics-file",
  default='roc-statistics.txt',
  dest="roc_statistics_file",
  help="Filename to write ROC statistics to"
)
option_parser.add_option(
  "--num-points",
  dest="num_points",
  default="16",
  help="The number of points to generate on the ROC curve"
)
option_parser.add_option(
  "--max-p-binding-site",
  dest="max_p_binding_site",
  default="0.1",
  help="The largest p(binding site) to use in Viterbi algorithm, will generate right-most point on ROC."
)
option_parser.add_option(
  "--min-p-binding-site",
  dest="min_p_binding_site",
  default="0.0001",
  help="The smallest p(binding site) to use in Viterbi algorithm, will generate left-most point on ROC."
)
# sys.argv='test_harness.py -m test-method -n background-mouse'.split()
options, args = option_parser.parse_args()
logging.info('Model directory:                            %s', options.model_dir)
logging.info('ROC output file:                            %s', options.roc_output_file)
if not options.roc_title:
    options.roc_title = options.roc_output_file
logging.info('ROC title:                                  %s', options.roc_title)
logging.info('ROC statistics file:                        %s', options.roc_statistics_file)
if None == options.model_dir:
    raise RuntimeError('No model directory specified in options')
logging.info('# ROC points:                               %s', options.num_points)
logging.info('Max p(binding site):                        %s', options.max_p_binding_site)
logging.info('Min p(binding site):                        %s', options.min_p_binding_site)
if options.negative_dir:
    logging.info('Dataset specific negative sequences dir:    %s', options.negative_dir)
else:
    logging.info('Negative sequence set:                      %s', options.negative_sequences)



#
# get the sequence sets
#
def sequence_name_filename_map(
  test_harness_seq_dir,
  num_fold_cross_validation=5
):
    def make_filename(tag):
        return os.path.join(test_harness_seq_dir, '%s.fa' % tag)
    class MakeFilenameDict(dict):
        def __missing__(self, k):
            filename = make_filename(k)
            if os.access(filename, os.R_OK):
                self[k] = filename
            else:
                raise RuntimeError('No such file: %s' % filename)
            return self[k]
    sequence_files = MakeFilenameDict()
    def fragment_validate_pair(fragment, i):
        return (
          '%s-%d-validate' % (fragment, i),
          make_filename('vm-%s-test-x%d-corrected' % (fragment, i))
        )
    def fragment_training_pair(fragment, i):
        return (
          '%s-%d-training' % (fragment, i),
          make_filename('vm-%s-train-x%d-corrected' % (fragment, i))
        )
    def add_fragment(fragment):
        for i in xrange(1, num_fold_cross_validation+1):
            sequence_files.update((
              fragment_validate_pair(fragment, i),
              fragment_training_pair(fragment, i),
            ))
    for fragment in all_fragments:
        add_fragment(fragment)
    add_fragment('T99999')
    return sequence_files

if os.name == 'nt':
    test_harness_seq_dir = 'C:/Dev/MyProjects/Bio/hmm/python/test-harness-seqs'
else:
    test_harness_seq_dir = '/home/kenneth/vm-intermediate-results'
    if not os.path.exists(test_harness_seq_dir):
        test_harness_seq_dir = '/home/reid/Analysis/GappedPssms/fragments/cross-validate/'
sequence_files = sequence_name_filename_map(test_harness_seq_dir)
class SequenceMap(dict):
    def __missing__(self, k):
        self[k] = [
          seq_to_numpy(s)
          for s
          in convert_seqs(sequence_files[k])
        ]
        return self[k]
sequence_map = SequenceMap()



logging.info('Loading sequences')
if options.negative_dir:
    sequence_dict = PositiveNegativeSequences(NegativeSeqsFromDir(options.negative_dir))
else:
    sequence_dict = PositiveNegativeSequences(NegativeSeqsFromFasta(sequence_files[options.negative_sequences]))


#
# Find the models that the algorithm has generated
#
# for each file called *.pssm in the model directory
logging.info('Looking for PSSMs in: %s' % options.model_dir)
pssms_for_method = dict() # a dict of parsed models indexed by dataset/cross-validation set index tuple
for file in os.listdir(options.model_dir):
    base, ext = os.path.splitext(file)
    if '.pssm' == ext:
        dataset, cross_fold_idx = base.split('-')
        cross_fold_idx = int(cross_fold_idx)
        logging.info(
          'Reading models from: %s for dataset: %s and cross-fold validation set: %d',
          file,
          dataset,
          cross_fold_idx
        )
        pssms_for_method[(dataset, cross_fold_idx)] = list(
          parse_models(
            open(os.path.join(options.model_dir, file))
          )
        )






#
# Generate a ROC point for each p_binding_site for these PSSMs
#
#p_binding_sites = [.1]
num_roc_points = int(options.num_points)
min_p_binding_site = float(options.max_p_binding_site)
max_p_binding_site = float(options.min_p_binding_site)
log_min_p_binding_site = numpy.log10(min_p_binding_site)
log_max_p_binding_site = numpy.log10(max_p_binding_site)
p_binding_sites = numpy.power(
  10,
  numpy.arange(
    log_min_p_binding_site,
    log_max_p_binding_site,
    (log_max_p_binding_site-log_min_p_binding_site)/num_roc_points
  )
)[:num_roc_points]
#p_binding_sites.reverse()
logging.info('Using the following p(binding) values: %s', str(p_binding_sites))
rocs = dict(
  (p_binding_site, get_roc_for(p_binding_site, pssms_for_method, sequence_dict))
  for p_binding_site
  in p_binding_sites
)
logging.info('Generated all ROC points.')
f = open(options.roc_statistics_file, 'w')
f.write('# p(binding);TP;TN;FP;FN\n')
for p_binding_site in p_binding_sites:
    roc = rocs[p_binding_site][0]
    logging.info('Complete ROC for p(binding site) = %f:\n%s' % (p_binding_site, str(roc)))
    f.write('%f;%d;%d;%d;%d\n' % (p_binding_site, roc.tp, roc.tn, roc.fp, roc.fn))
f.close()

def plot_rocs(rocs, filename, fig_title):
    from pylab import figure, close, savefig, plot, title, xlabel, ylabel
    figure()
    plot_roc_points(
      rocs,
      marker='s',
      color='black',
      linestyle='--'
    )
    plot(
      [0,1],
      [0,1],
      color='black',
      linestyle=':'
    )
    title(fig_title)
    xlabel('1 - specificity: 1-TN/(TN+FP)')
    ylabel('sensitivity: TP/(TP+FN)')
    logging.info('Saving ROC to %s' % filename)
    savefig(filename)

# plot the ROC curves
plot_rocs([rocs[p][0] for p in p_binding_sites[::-1]], options.roc_output_file, options.roc_title)
for dataset in rocs[p_binding_sites[0]][1]:
    plot_rocs([rocs[p][1][dataset] for p in p_binding_sites[::-1]], '%s' % dataset, 'ROC for %s' % dataset)
