#
# Copyright John Reid 2008
#

"""
Code to rank PSSMs by "interesting-ness".

Information content.
Low-order predictability.
Number of sequences with sites.
"""


from parse_gapped_format import parse_models
from itertools import imap
import glob, logging, shutil, os
from gapped_pssms.pssm_score import *
from gapped_pssms.data import test_set_fragments

logging.basicConfig(level=logging.DEBUG)


def calculate_emissions(model):
    emissions = numpy.zeros((model.N, model.M))
    for i in xrange(model.N):
        assert model.emissions[i][0] == i
        emissions[i] = model.emissions[i][1]
    return emissions


M = 4
for fragment in test_set_fragments:
    for cross_fold in xrange(1, 6):
        logging.info('%s %d', fragment, cross_fold)
        results = list()
        for pssm_file in glob.glob('typical-pssms/%s-%d-*.pssm' % (fragment, cross_fold)):
            models = parse_models(open(pssm_file))
            for model in models:
                assert model.M == M
                emissions = calculate_emissions(model)
                first_order_entropy_score = calculate_first_order_entropy_score(emissions)
                information_content_score = calculate_information_content_score(emissions)
                overall_score = geometric_mean((first_order_entropy_score, information_content_score))
                results.append((overall_score, pssm_file))

        results.sort(reverse=True)
        files = []
        for i, (score, file) in enumerate(results):
            src = file.replace('.pssm', '.png')
            dest = 'typical-pssms/rescored/%s-%d-%03d.png' % (fragment, cross_fold, i)
            files.append(dest)
            shutil.copy(src, dest)
        montage_file = 'typical-pssms/rescored/montages/%s-%d.png' % (fragment, cross_fold)
        montage_cmd = 'montage -tile 1x -geometry 1200 %s %s' % (' '.join(files), montage_file)
        os.system(montage_cmd)
