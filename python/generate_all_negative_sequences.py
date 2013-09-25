#
# Copyright John Reid 2008
#

"""
Code to generate negative test sequences for those fragments in the test harness.
"""

from gapped_pssms.data import fasta_file_for_fragment, test_set_fragments
import sys

def sequences_from_fasta(fasta):
    """Yields sequences from fasta file."""
    import corebio.seq_io.fasta_io
    from itertools import imap
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



for fragment in test_set_fragments:
    seqs = list(sequences_from_fasta(fasta_file_for_fragment(fragment)))
    seq_length = max(len(s) for s in seqs)
    num_seqs = len(seqs)
    sys.argv = ('generate_negative_test_sequences.py -m ..\..\Python\%s-bg-model.pickle -n %d -l %d -o negative-%s.fa' % (
      fragment,
      num_seqs,
      seq_length,
      fragment,
    )).split()
    execfile('generate_negative_test_sequences.py')
