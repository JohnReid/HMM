#
# Copyright John Reid 2007
#

"""
Code to run MEME on sequences and convert output to HMM
"""


import os, numpy
from TAMO.MD.MEME import Meme

def run_system(cmd):
    result = os.system(cmd)
    if result:
        raise RuntimeError('Could not run (%d): %s' % (result, cmd))

def run_cygwin(cmd):
    dos_cmd = 'C:/apps/cygwin/bin/bash -c "%s"' % cmd
    print 'Running cygwin: %s' % cmd
    run_system(dos_cmd)

def run_meme(
        input_file,
        output_file,
        maxiter=None,
        distance=None,
        minw=12,
        maxw=22,
        nmotifs=None,
        mod='zoops',
        maxsize=30000000,
        extra_args=None
):
    cmd = '/home/ubcg52d/bin/meme %s -dna -revcomp' % input_file
    if None != maxiter:
        cmd += ' -maxiter %d' % maxiter
    if None != distance:
        cmd += ' -distance %f' % distance
    if None != minw:
        cmd += ' -minw %d' % minw
    if None != maxw:
        cmd += ' -maxw %d' % maxw
    if None != nmotifs:
        cmd += ' -nmotifs %d' % nmotifs
    if None != mod:
        cmd += ' -mod %s' % mod
    if None != maxsize:
        cmd += ' -maxsize %d' % maxsize
    if None != maxiter:
        cmd += ' -maxiter %d' % maxiter
    if None != extra_args:
        cmd += ' %s' % ' '.join(str(a) for a in extra_args)
    cmd += ' > %s' % output_file
    try:
        run_cygwin(cmd)
    except:
        if os.access(output_file, os.R_OK):
            os.remove(output_file)
        raise

def meme_motif_base_dist(base_dist):
    return [
            base_dist['A'],
            base_dist['C'],
            base_dist['G'],
            base_dist['T'],
    ]

def meme_motif_dist(motif):
    return numpy.array(
            [
              meme_motif_base_dist(bd) for bd in motif.P
            ]
    )

if '__main__' == __name__:
    output_file = 'meme_output.txt'
    meme = Meme(output_file)
    for motif in meme.motifs:
        dist = meme_motif_dist(motif)
        print dist
