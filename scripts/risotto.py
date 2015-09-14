#
# Copyright John Reid 2008
#

"""
Methods to call risotto, the suffix trie mismatch motif finder.
"""

import logging, os, os.path



def _write_alphabet_file(filename):
    f = open(filename, 'w')
    f.write(
  """Type:Nucleotides
  A
  C
  G
  T
  """
    )
    f.close()

def _write_fasta_file(filename, sequences):
    from hmm.pssm import numpy_to_seq
    f = open(filename, 'w')
    for i, s in enumerate(sequences):
        f.write('> %i\n')
        f.write(s)
        f.write('\n')
    f.close()


def _parse_output(filename):
    f = open(filename)
    while not f.next().startswith('='):
        True
    results = []
    while True:
        line = f.next().strip()
        if line.startswith('Nb models'):
            tmp, num_models = line.split(': ')
            num_models = int(num_models)
            if num_models != len(results):
                raise RuntimeError('Problem parsing risotto output: %s' % filename)
            break
        results.append(line.split())
    return results

def risotto(
  sequences,
  quorum=40,
  K_min=10,
  K_max=10,
  max_mismatches=2
):
    """
    Calls risotto on the sequences.

    @arg sequences: The sequences.
    @arg K_min: The minimum length of motif to look for.
    @arg K_max: The minimum length of motif to look for.
    @arg quorum: The percentage of sequences we must find the motif in.
    @arg max_mismatches: The maximum number of mismatches.
    @return: A list of (consensus, ?, # sequences, # occurences) tuples, one for each motif found.
    """

    from os import tempnam

    alphabet_file = tempnam()
    logging.debug('Risotto alphabet file: %s', alphabet_file)
    _write_alphabet_file(alphabet_file)

    fasta_file = tempnam()
    logging.debug('Risotto fasta file: %s', fasta_file)
    _write_fasta_file(fasta_file, sequences)

    output_file = tempnam()
    logging.debug('Risotto output file: %s', output_file)

    try:
        cmd = 'max_extensibility.exe %s %s %s %f 1 %d %d %d' % (
          alphabet_file,
          fasta_file,
          output_file,
          quorum,
          K_min,
          K_max,
          max_mismatches
        )
        if os.system(cmd):
            raise RuntimeError('Could not run risotto with args: %s' % cmd)

        return _parse_output(output_file)

    finally:
        if os.access(alphabet_file, os.R_OK):
            os.remove(alphabet_file)
        if os.access(fasta_file, os.R_OK):
            os.remove(fasta_file)
        if os.access(output_file, os.R_OK):
            os.remove(output_file)




if '__main__' == __name__:
    logging.basicConfig(level=logging.DEBUG)

    from sequences import convert_seqs, rev_comp, numpy_to_seq, seq_to_numpy
    test_seq_file = 'synthetic-2/synthetic-sequences-K10-g0.50-N200-L200-seed4-0.fa'
    sequences = convert_seqs(test_seq_file)
    #sequences.extend([numpy_to_seq(rev_comp(seq_to_numpy(s))) for s in sequences])
    #sequences[0] = sequences[0].replace('a','n')
    results = risotto(sequences, max_mismatches=1)
    results.sort(key=lambda x: x[2]) # sort by number of sequences that motif is found in
    print results[-10:]
