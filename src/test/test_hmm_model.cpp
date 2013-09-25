
#ifdef _MSC_VER
# pragma warning( disable : 4172 )
#endif //_MSC_VER

#include "test.h"

#include <hmm/model.h>
#include <hmm/forward_backward.h>
#include <hmm/viterbi.h>
#include <hmm/baum_welch.h>
#include <hmm/sample.h>
using namespace hmm;

#include <myrrh/math.h>

namespace boost {
    void
    assertion_failed_msg(
        char const * expr,
        char const * msg,
        char const * function,
        char const * file,
        long line
    ) {
        std::cerr
            << file << ":" << line
            << ", " << function << "()" 
            << "BOOST_ASSERT_MSG(" << expr << ") failed, \""
            << msg << "\"\n";
    }
}




void check_empty_sequences()
{
	cout << "******* check_empty_sequences() " << endl;

	using namespace hmm;

	model::ptr m = create_random_fully_connected_model(3, 4, 1);
	const output_seq O; //empty output sequence

	//forward-backward
	double_array alpha;
	double_array beta;
	double_vec c;
	forward( *m, O, alpha, c );
	forward_LL( *m, O );
	backward( *m, O, c, beta );
	LL_from_backward_probabilities( *m, O, c, beta );

	//viterbi
	unsigned_vec q_star;
	viterbi( *m, O, q_star );

	//baum-welch
	std::vector< output_seq > output_seqs;
	output_seqs.push_back( O );
	model::prior prior( m->N, m->M );
	baum_welch_multi_sequence_iteration( *m, prior, output_seqs );
}

void check_sample()
{
	cout << "******* check_sample() " << endl;

	using namespace hmm;

	unsigned_vec states;
	output_seq sampled_sequence;
	model::ptr m = create_random_fully_connected_model(3, 64, 2);
	sample( *m, 1000, states, sampled_sequence );
	BOOST_CHECK( m->converter.is_valid_order_0_sequence( sampled_sequence.begin(), sampled_sequence.end() ) );
}

test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite * test = BOOST_TEST_SUITE( "HMM test suite" );

	try
	{
		test->add( BOOST_TEST_CASE( &check_sample ), 0 );
		test->add( BOOST_TEST_CASE( &check_empty_sequences ), 0 );
	}
	catch (const std::exception & e)
	{
		cerr << "Registering tests - exception: " << e.what() << endl;
	}

    return test; 
}

