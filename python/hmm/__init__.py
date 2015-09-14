#
# Copyright John Reid 2007, 2008, 2015
#

"""
Code that wraps C++ HMM python module L{hmm._hmm}.
Also contains code useful in its own right.
"""


from ._build import *
import os

if os.name == 'nt':
    try:
        import _fpu as fpu
    except:
        from warnings import warn
        warn('Could not import FPU module: Only used on Windows')

def inf():
    return float("1e9999")

def nan():
    return inf()-inf()

def module_info():
    "@return: A string giving information about the module."
    return 'HMM python module: %s (%s)' % (
            _hmm.__version__,
            _hmm.__debug__ and 'debug' or 'release'
    )

def as_state_model(model):
    """
    If model is of type ModelByStates, returns it, otherwise returns
    model_2_model_states(model)
    """
    if isinstance(model, ModelByStates):
        return model
    else:
        return model_2_model_states(model)


def as_model(model):
    """
    If model is of type Model, returns it, otherwise returns
    model_states_2_model(model)
    """
    if isinstance(model, Model):
        return model
    else:
        return model_states_2_model(model)

def graph_model(
        model,
        include_emissions = False,
        state_labels = None, # function mapping state indices to labels
        show_dists = False,
        label_edges = True
):
    """
    Returns a boost.graph.Digraph of the model

    if include_emissions is true also includes nodes for emissions

    if state_labels is a function it is called for each node to return a label
            nodes can be groupd as one by returning the same label

    if show_dists is a function it is called for each label to see if
            we should plot the distribution for that labelled node
    """
    import boost.graph as bgl
    g = bgl.Digraph()
    label = g.add_vertex_property( 'label', 'string' )
    if label_edges:
        edge_label = g.add_edge_property( 'label', 'string' )
    shape = g.add_vertex_property( 'shape', 'string' )
    if show_dists: shapefile = g.add_vertex_property( 'shapefile', 'string' )

    # if we have not been supplied with labels for the states use the indices
    if None == state_labels:
        def state_labels( i ):
            return str(i)

    # add each label as a vertex
    vertices = dict(
            [
                    (l,None)
                    for l
                    in set(
                            [
                                    state_labels(i)
                                    for i in xrange( model.N )
                            ]
                    )
            ]
    )
    for l in vertices.keys():
        v = g.add_vertex()
        vertices[ l ] = v
        label[ v ] = l

    # add the transitions
    for i in xrange( model.N ):
        for j in xrange( model.N ):
            params = model.get_transition_parameterisation( i, j )
            if None != params:
                p = model.a( i, j )
                if p > 0.0:
                    v1 = vertices[ state_labels( i ) ]
                    v2 = vertices[ state_labels( j ) ]
                    e = g.edge( v1, v2 )
                    if not e:
                        e = g.add_edge( v1, v2 )
                        if label_edges:
                            edge_label[ e ] = '0.0'
                        if label_edges:
                            edge_label[ e ] = '%.3f' % (p + float( edge_label[ e ] ) )


    # add the emissions if asked to
    if include_emissions:
        for i in xrange( model.N ):
            for k in xrange( model.M ):
                #only include if '+'ve probability
                p = model.b( i, k )
                if p > 0.0:
                    v = g.add_vertex()
                    label[ v ] = '%d' % k
                    shape[ v ] = 'diamond'

                    e = g.add_edge( vertices[ state_labels( i ) ], v )
                    if label_edges:
                        edge_label[ e ] = '%.3f' % p

    # show distributions if asked
    if show_dists:
        from hmm.weblogo import format_weblogo_from_dist
        import numpy
        B = model.B
        for l, v in vertices.iteritems():
            if show_dists( l ):
                dist = numpy.array(
                        [
                                B[i]
                                for i in xrange( model.N )
                                if state_labels(i) == l
                        ]
                )
                basename = 'emission_%s' % l
                format_weblogo_from_dist(
                        dist,
                        basename,
                        'png',
                        convert_args = '-scale 30x90' )
                shapefile[ v ] = basename + '.png'
                label[ v ] = ''

    return g

def graph_as_svg(
        model,
        basename,
        directory = None,
        graphing_keywords = { },
        neato_properties = { },
        output_format = 'svg'
):
    """
    Writes the model as a svg file
    """
    import os
    current_dir = os.getcwd()
    if directory:
        if not os.access(directory, os.X_OK): os.makedirs(directory)
        os.chdir(directory)

    try:
        g = graph_model( model, **graphing_keywords )
        dot_filename = '%s.dot' % basename
        svg_filename = '%s.%s' % (basename, output_format)
        g.write_graphviz( dot_filename )
        cmd = 'neato %s -T%s -o %s -Gscale=overlap %s' % (
                dot_filename,
                output_format,
                svg_filename,
                " ".join( [ '%s=%s' % (k,v) for k, v in neato_properties.iteritems() ] )
        )
        #print 'Running: %s' % cmd
        os.system( cmd )
    finally:
        if directory: os.chdir( current_dir )

def emission_entropy( model ):
    "Entropy of each state's emission distribution"
    def entropy( row ):
        import math
        sum = 0.0
        for x in row:
            try:
                sum -= x * math.log( x )
            except OverflowError:
                pass # this is ok - caused by tiny x
        return sum
    import numpy
    return numpy.array( [ entropy( row ) for row in model.B ] )


class ReverseComplementCollapsingCounter(object):
    """
    Counts n-mers but collapses reverse complement counts together.
    """
    def __init__(self, n):
        self._n = n
        self._counts = dict()
        self._converter = MarkovOrderConverter(4, n-1)

    def __call__(self, nmer, count):
        assert len(nmer) == self._n
        import hmm.pssm
        rev_comp_nmer = hmm.pssm.rev_comp(nmer)
        rev_comp_idx = self._converter.convert_to_order_n_observation(rev_comp_nmer.astype(int))
        if rev_comp_idx in self._counts:
            self._counts[rev_comp_idx] += count
        else:
            idx = self._converter.convert_to_order_n_observation(nmer.astype(int))
            self._counts[idx] = count

    def num_counts(self):
        return len(self._counts)

    def counts(self):
        for idx, count in self._counts.iteritems():
            yield self._converter.convert_from_order_n_observation(idx), count
