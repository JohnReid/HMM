#
# Copyright John Reid 2007
#


"""
Code to manipulate sequences
"""

import numpy, math, os, os.path

from hmm.pssm import *

if os.name == 'nt':
    fragments_dir = 'c:/Analysis/GappedPssms/fragments/'
    test_set_dir = os.path.join(fragments_dir, 'test-sets')
else:
    fragments_dir = '/home/reid/Dev/MyProjects/Bio/hmm/python/fragments'
    test_set_dir = os.path.join(fragments_dir, 'test-sets')


def fragment_filename( fragment ):
    return os.path.join(fragments_dir, '%s.fa.masked' % fragment)

def fragment_test_set_filename( fragment, set_idx = 0 ):
    return os.path.join(test_set_dir, fragment, '%scorrectedRM-testset-%d.fasta' % (fragment, set_idx))

def fragment_training_set_filename( fragment, set_idx = 0 ):
    return os.path.join(test_set_dir, fragment, '%scorrectedRM-trainingset-%d.fasta' % (fragment, set_idx))

def convert_seqs( filename ):
    import corebio.seq_io.fasta_io
    return [
            str(s).strip('nN')
            for s
            in corebio.seq_io.fasta_io.iterseq(
                    open( filename, 'r' ),
                    corebio.seq.dna_alphabet
            )
    ]

def seqs_for_fragment(
        fragment,
        verbose = False,
        use_newer_files = True
):
    if verbose: print 'Fragment = %s' % fragment

    if use_newer_files:
        test_seqs = convert_seqs(fragment_test_set_filename(fragment))
        training_seqs = convert_seqs(fragment_training_set_filename(fragment))

        import itertools
        seqs = [
                seq_to_numpy( seq )
                for seq in itertools.chain(test_seqs, training_seqs)
        ]
        test_seqs = training_seqs = None # discard unneeded sequences

    else: # use older fasta files
        bio_seqs = convert_seqs(fragment_filename(fragment))
        seqs = [
                seq_to_numpy( seq )
                for seq in bio_seqs
        ]
        bio_seqs = None # discard unneeded bio object

    if verbose:
        print '%d sequences' % len( seqs )
        print '%d bases' % num_bases( seqs )

    return seqs

def informative_states(
        model,
        information_threshold = 0.5
):
    """
    Identifies those states in the model that have an infomation content higher
    than the threshold

    Returns: (
            informative states,
            bits of information in the model
    )
    """
    informative_states = set()
    B = model.B
    bits = 0.0
    for i in xrange( model.N ):
        #print B[ i,: ]
        ic = math.log( 4, 2 )
        for x in xrange( 4 ):
            if B[ i, x ] > 0:
                ic += B[ i, x ] * math.log( B[ i, x ], 2 )
        #print i, ic
        if ic > information_threshold: informative_states.add( i )
        bits += ic
    return informative_states, bits

def remove_sites_from_seqs(
        model,
        seqs,
        states_to_remove,
        start_states
):
    """
    Removes sites from the sequences by changing the data to unknown

    Returns: (
            number of sites found,
            number bases removed
    )
    """
    num_sites = 0
    removed = 0
    sites = [ ]
    for s, seq in enumerate( seqs ):
        log_P_star, q_star = model.viterbi( seq )
        for i, q in enumerate( q_star ):
            if q in start_states:
                num_sites += 1
                sites.append( ( s, i ) )
            if q in states_to_remove:
                removed += 1
                seq[ i ] = model.M #make this base unknown
    return num_sites, removed


def get_all_sequences_for_set_and_idx(filename_fn, test_set_idx):
    from itertools import chain
    def seqs_by_fragment():
        for fragment in all_fragments:
        #for fragment in [ 'T00594' ]:
            for s in seqs_to_numpy(convert_seqs(filename_fn(fragment, test_set_idx))):
                yield s
    return [ s for s in seqs_by_fragment() ]

def get_all_training_sequences(test_set_idx):
    return get_all_sequences(fragment_training_set_filename, test_set_idx)

def get_all_validation_sequences(test_set_idx):
    return get_all_sequences(fragment_test_set_filename, test_set_idx)

def get_all_sequences():
    from itertools import chain
    print 'Getting all sequences'
    seqs = [
            s for s in
            chain(
                    *(
                            get_all_sequences_for_set_and_idx(fragment_training_set_filename, 0),
                            get_all_sequences_for_set_and_idx(fragment_test_set_filename, 0)
                    )
            )
    ]
    print 'Got them'
    return seqs

def write_seqs_as_fasta(seqs, filename):
    f = open(filename, 'w')
    for i, seq in enumerate(seqs):
        f.write('> sequence %d' % (i+1))
        f.write('\n')
        f.write(numpy_to_seq(seq))
        f.write('\n')
