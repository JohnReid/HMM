#
# Copyright John Reid 2007,2008
#

"""
Code to read and write HMMs.
"""

from __future__ import with_statement
import hmm

class ModelBuilder(object):
    def __init__(self):
        self.model = hmm.ModelByStates()

    def __call__(self, l):
        try:
            # ignore empty and commented lines
            if not len(l): return
            if l.startswith('#'): return

            fields=l.split(':')
            if 'MODEL' == fields[0]: self.parse_model_definition(fields[1])
            elif 'TRANSITION' == fields[0]: self.parse_transition(fields[1])
            elif 'EMISSIONS' == fields[0]: self.parse_emissions(fields[1])
            else: RuntimeError('Cannot parse model: ' + fields[1])
        except:
            import sys
            raise RuntimeError('Cannot parse: "%s"\n%s' % (l, str(sys.exc_info()[1])))

    def parse_model_definition(self, match):
        # expect model definition line to be "MODEL:N,M"
        fields = match.split(',')
        if 2 != len(fields): raise RuntimeError( 'Expecting line like "MODEL:N,M": ' + match )
        N = int(fields[0])
        M = int(fields[1])
        for i in xrange(N): self.model.add_state(M,self.model.add_parameter(1.0/N))

    def parse_transition(self, match):
        if not self.model.N: raise RuntimeError('Must have model definition line before transitions')
        i, j, p = match.split(';')
        from_ = self.model.states[int(i)]
        to_ = self.model.states[int(i)]
        from_.add_successor(to_, self.model.add_parameter(float(p)))

    def parse_emissions(self, match):
        if not self.model.N: raise RuntimeError('Must have model definition line before emissions')
        fields = match.split(';')
        state = self.model.states[int(fields[0])]
        emissions = fields[1].split(',')
        if len(emissions) != self.model.M: raise RuntimeError( 'Wrong number of emissions: ' + emissions )
        for i, b in enumerate(emissions):
            state.b[i] = self.model.add_parameter(float(b))

def build_model(lines):
    "Takes a sequence of lines and builds a model from them"
    builder = ModelBuilder()
    for l in lines: builder(l)
    return builder.model

def write_model(model, stream, write_state = None):
    """
    Writes a model to the stream

    if None != write_state then calls write_state(i) for each state to
    determine if that state is output
    """

    # get a map from states we will output to indexes
    state_map = dict();     j = 0
    for state in model.states:
        if not write_state or write_state(state): state_map[state] = j; j += 1

    stream.write('MODEL:%d,%d\n\n' % (j,model.M))

    for state in model.states:
        if state in state_map:
            i = state_map[state]
            for successor in state.successors:
                if successor.state in state_map:
                    j = state_map[successor.state]
                    p = model.parameters[successor.a.idx]
                    stream.write('TRANSITION:%d;%d;%f\n' % (i,j,p))
    stream.write('\n')

    for state in model.states:
        if state in state_map:
            i = state_map[state]
            stream.write('EMISSIONS:%d;%s\n' % (i, ','.join([str(model.parameters[b.idx]) for b in state.b])))
    stream.write('\n')


if '__main__' == __name__:
    to_parse = """
# MODEL:N,M - where N is the number of states and M the size of the output alphabet
MODEL:9,4

# TRANSITION:i;j,p - where p(transition i to j) = p
TRANSITION:0;1;1
TRANSITION:1;2;.5
TRANSITION:1;3;.5
TRANSITION:2;3;1
TRANSITION:3;4;.1
TRANSITION:3;5;.9
TRANSITION:4;5;1
TRANSITION:5;6;1
TRANSITION:6;7;1
TRANSITION:7;8;1

# EMISSIONS:i;b1,b2,b3,b4 - where p(state i emits base i) = bi
EMISSIONS:0;.5,.5,0,0
EMISSIONS:1;0,0,.5,.5
EMISSIONS:2;1,0,0,0
EMISSIONS:3;0,0,0,1
EMISSIONS:4;0,0,1,0
EMISSIONS:5;0,1,0,0
EMISSIONS:6;1,0,0,0
EMISSIONS:7;1,0,0,0
EMISSIONS:8;0,0,1,0
"""
    model = build_model( to_parse.split('\n') )
    m = hmm.model_states_2_model(model)
    if False:
        hmm.graph_as_svg(
                m,
                'test',
                'model_io_test',
                graphing_keywords = {
                        'show_dists' : lambda l: True,
                        'state_labels' : None
                },
                neato_properties = { '-Elen' : '2' }
        )
    with open('model.mdl', 'w') as f: write_model(model, f)
