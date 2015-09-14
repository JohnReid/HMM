#
# Copyright John Reid 2007
#

"""
Data specific to the TRANSFAC ChIP-chip data
"""

all_fragments = [
        'T00594',
        'T00163',
        'T00759',
        'T00368',
        'T03286',
        'T00671',
        'T00140',
        'T09363',
        'T03828',
        'T08969',
        'T00781'
]

fragment_details = '''
T00594:  156324/ 166306 (93%) known bases in  204 sequences
T00163:  168446/ 181045 (93%) known bases in  208 sequences
T00759:  213149/ 236168 (90%) known bases in  335 sequences
T00368:  231204/ 252945 (91%) known bases in  286 sequences
T03286:  252909/ 274166 (92%) known bases in  308 sequences
T00671:  350779/ 526136 (66%) known bases in  584 sequences
T00140:  466566/ 507565 (91%) known bases in  718 sequences
T09363: 1107182/1462020 (75%) known bases in 1031 sequences
T03828: 1738443/1876753 (92%) known bases in 2077 sequences
T08969: 2991217/4183054 (71%) known bases in 2885 sequences
T00781: 8075481/8588852 (94%) known bases in 8475 sequences
'''

if '__main__' == __name__:
    from sequences import *
    from count_mers import *
    from background_models import *
    import hmm.pssm
    for fragment in all_fragments:
        background_model = global_background_model_cache().get_model(order=0, num_mosaics=1, fragment=fragment)
        seqs = seqs_for_fragment(fragment, verbose=False)
        num_bases = hmm.pssm.num_bases(seqs)
        num_known_bases = hmm.pssm.num_known_bases(seqs)
        print '%s: %7d/%7d (%d%%) known bases in %4d sequences' % (fragment, num_known_bases, num_bases, (100*num_known_bases/num_bases), len(seqs))
        significant_mers = calculate_n_mer_significances([seq for seq in seqs if len(seq)], 5, None)
        print '\n'.join(str(mer) for mer in significant_mers[:10])
