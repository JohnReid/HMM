#
# Copyright John Reid 2008
#

"""
Code to test building PSSMs.
"""



import hmm, numpy, os

if '__main__' == __name__:

    model = hmm.build_pssm_model(
            0.04,
            numpy.array(
                    [
                            [ 0.97, 0.01, 0.01, 0.01 ],
                            [ 0.01, 0.01, 0.01, 0.97 ],
                            [ 0.01, 0.01, 0.97, 0.01 ],
                    ]
            )
    )


    O = numpy.array(
            [
                    0, 0, 0, 3, 2, 0, 0, 1, 0, 3, 0, 0
            ]
    )



    def test_graphing():
        print '***** Graphing *****'
        g = hmm.graph_model( model, include_emissions = True )
        g.write_graphviz( 'hmm.dot' )
        os.system( 'neato hmm.dot -Tsvg -o hmm.svg -Gscale=overlap -Elen=1.5' )


    def test_viterbi():
        print '***** Viterbi *****'
        log_P_star, q_star = model.viterbi( O )
        expected_result = numpy.array( [0, 0, 1, 2, 3, 0, 0, 6, 5, 4, 0, 0] )
        print log_P_star
        print q_star
        assert len( expected_result ) == len( q_star )
        for i, v in enumerate( q_star ):
            assert expected_result[ i ] == v



    def test_forward():
        print '***** Forward *****'
        LL = []
        for scaling_threshold in [ 1.0, 1e-5, 1e-10 ]:
            print 'Scaling threshold: %f' % scaling_threshold
            log_p_O, alpha, c = model.forward( O, scaling_threshold )
            LL.append( log_p_O )
            print log_p_O
            # print c
        for l1 in LL:
            for l2 in LL:
                assert l1 - l2 / l1 < 1e-8, 'log p(O) should be independent of scaling factor'



    def test_baum_welch():
        print '***** Baum-Welch *****'
        print model.baum_welch( [ O ] )
        for i in xrange( 10 ):
            LL = model.baum_welch_iteration( [ O ] )
            print LL
        g_after_bw = hmm.graph_model( model, include_emissions = True )
        g_after_bw.write_graphviz( 'hmm_after_baum_welch.dot' )
        os.system( 'neato hmm_after_baum_welch.dot -Tsvg -o hmm_after_baum_welch.svg -Gscale=overlap -Elen=1.5' )



    def test_sample():
        print '***** Sample *****'
        length = 30
        states, output = model.sample( length )
        assert length == len( states )
        assert length == len( output )
        print states
        print output



    test_graphing()
    test_viterbi()
    test_forward()
    test_baum_welch()
    test_viterbi() # check this after baum-welch
    test_forward() # check this after baum-welch
    test_sample()
