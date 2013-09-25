
#include <hmm/defs.h>
#include "test.h"

#include <hmm/markov_order.h>
using namespace hmm;

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



typedef std::vector< unsigned > unsigned_vec;
typedef std::vector< char > char_vec;

#define ALPHABET_SIZE 5

const unsigned_vec markov_orders = list_of
	(   0 )
	(   1 )
	(   2 )
	(   3 )
	(   4 )
	(   5 )
	(  10 )
	//(  15 ) too big for unsigned
	//( 100 )
	//( 999 )
	;

const unsigned_vec test_seq_order_0 = list_of
	( 0 )
	( 1 )
	( 2 )
	( 3 )
	( 4 )
	;

template< unsigned AlphabetSize, typename SeqT >
void
check_conversion_to_and_from( const SeqT & input_seq, unsigned n )
{
	typedef typename SeqT::value_type char_t;
	typedef markov_order_n_converter< char_t > converter_t;
	converter_t converter( AlphabetSize, n );

	unsigned_vec seq_order_n;
	converter.make_markov_order_n( input_seq.begin(), input_seq.end(), back_inserter( seq_order_n ) );
	BOOST_CHECK( converter.is_valid_order_n_sequence( seq_order_n.begin(), seq_order_n.end() ) );
	BOOST_CHECK_EQUAL( seq_order_n.size(), input_seq.size() );
	
	unsigned_vec seq_order_0;
	converter.make_markov_order_0( seq_order_n.begin(), seq_order_n.end(), back_inserter( seq_order_0 ) );
	BOOST_CHECK_EQUAL( seq_order_0.size(), seq_order_n.size() );

	for( unsigned i = 0; input_seq.size() != i; ++i )
	{
		if( converter.order_n_size == seq_order_n[i] ) BOOST_CHECK_EQUAL( converter.order_0_size, seq_order_0[i] );
		else BOOST_CHECK_EQUAL( seq_order_0[i], input_seq[i] );
	}
}

void
check_markov_order( unsigned n )
{
	cout << "******* check_markov_order(): " << n << endl;
	
	check_conversion_to_and_from< ALPHABET_SIZE >( test_seq_order_0, n );
}



test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite * test = BOOST_TEST_SUITE( "Markov order test suite" );

	try
	{
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_markov_order, 
				markov_orders.begin(),
				markov_orders.end() ),
			0 );
	}
	catch (const std::exception & e)
	{
		cerr << "Registering tests - exception: " << e.what() << endl;
	}

    return test; 
}

