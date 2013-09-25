#
# Copyright John Reid 2009
#

"""
Tests HMM Viterbi algorithm using new no-copy data interface.
"""

import _hmm, numpy, os

print '_hmm.__debug__ = %s' % _hmm.__debug__

# create model
M=4
model_by_states = _hmm.ModelByStates(M=M, markov_order=0)

# add background state
bg = model_by_states.add_state()
bg.pi = model_by_states.add_parameter(1.)
uniform_param = model_by_states.add_parameter(.25)
for m in xrange(bg.M):
    bg.b[m] = uniform_param

specific_states = [model_by_states.add_state() for i in xrange(M)]

# connect background to initial binding site states
p_binding_site = 0.1
binding_param = model_by_states.add_parameter(p_binding_site)
not_binding_param = model_by_states.add_parameter(1.-p_binding_site)
bg.add_successor(specific_states[0], binding_param)
bg.add_successor(bg, not_binding_param)
always_one_param = model_by_states.add_parameter(1.)
for i in xrange(M-1):
    specific_states[i].add_successor(specific_states[i+1], always_one_param)
specific_states[M-1].add_successor(bg, always_one_param)

# set up specific emissions
likely_emission_param = model_by_states.add_parameter(1.-(M-1)*.02)
unlikely_emission_param = model_by_states.add_parameter(.02)
for i, specific_state in enumerate(specific_states):
    for b in xrange(M):
        if b == i:
            specific_state.b[b] = likely_emission_param
        else:
            specific_state.b[b] = unlikely_emission_param

model = _hmm.model_states_2_model(model_by_states)
model.normalise()

seq = numpy.array(
    [
        1, 1, 0, 1, 2, 3, 0, 3, 4, 0, 1, 2, 3, 3, 2, 3, 1
    ],
    dtype=numpy.uint16
)

LL1, s1 = model.viterbi_2(seq)
print LL1, s1

#
# Try the types it should work for
#
for t in [
    numpy.int8,
    numpy.int16,
    numpy.int32,
    numpy.int64,
    numpy.uint8,
    numpy.uint16,
    numpy.uint32,
    numpy.uint64,
]:
    seq_othertype = numpy.array(
        [
            1, 1, 0, 1, 2, 3, 0, 3, 4, 0, 1, 2, 3, 3, 2, 3, 1
        ],
        dtype=t
    )
    LL2, s2 = model.viterbi_2(seq_othertype)
    assert (s1 == s2).all()
    assert LL1 == LL2

#
# Try the types it shouldn't work for
#
for t in [
    numpy.bool,
    numpy.float64,
    numpy.complex128,
]:
    try:
        seq_othertype = numpy.array(
            [
                1, 1, 0, 1, 2, 3, 0, 3, 4, 0, 1, 2, 3, 3, 2, 3, 1
            ],
            dtype=t
        )
        LL1, s1 = model.viterbi_2(seq_othertype)
    except:
        continue
    # should not get here - exception should have been raised
    raise RuntimeError('Should not be able to run Viterbi algorithm on numpy array of this type. %s' % t)
