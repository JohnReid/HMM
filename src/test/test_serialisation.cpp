
#include <hmm/defs.h>
#include "test.h"

#include <hmm/model.h>
using namespace hmm;

#include <boost/filesystem/convenience.hpp>

#include <fstream>

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


#define ALPHABET_SIZE 3
#define NUM_STATES 2

const std::vector< unsigned > markov_orders = list_of
	(   0 )
	(   1 )
	(   2 )
	(   3 )
	(   4 )
	(   5 )
	//(  15 ) too big for unsigned
	//( 100 )
	//( 999 )
	;

template< typename T >
void check_serialisation_equals( T & t ) {
	const std::string filename = "test_serialisation.test";

	// save data to archive
    {
		std::ofstream ofs( filename.c_str() );
        boost::archive::text_oarchive oa( ofs );
        // write class instance to archive
        oa << t;
    	// archive and stream closed when destructors are called
    }

    // ... some time later restore the class instance to its orginal state
    T t_copy;
    {
        // create and open an archive for input
        std::ifstream ifs( filename.c_str() );
        boost::archive::text_iarchive ia( ifs );
        // read class state from archive
        ia >> t_copy;
        // archive and stream closed when destructors are called
    }
	boost::filesystem::remove( filename );

	BOOST_CHECK_EQUAL( t, t_copy );
}

typedef boost::multi_array< double, 2 > array_t;
typedef boost::optional< int > optional_t;
BOOST_TEST_DONT_PRINT_LOG_VALUE( optional_t );
BOOST_TEST_DONT_PRINT_LOG_VALUE( array_t );

void check_optional_serialisation() {
	cout << "******* check_optional_serialisation()" << endl;

	optional_t empty;
	optional_t not_empty( 1 );
	check_serialisation_equals( empty );
	check_serialisation_equals( not_empty );
}


void check_multi_array_serialisation() {
	cout << "******* check_multi_array_serialisation()" << endl;

	array_t a( boost::extents[ 2 ][ 2 ] );
	a[0][0] = 1.;
	a[1][0] = 2.;
	a[0][1] = 3.;
	a[1][1] = 4.;
	check_serialisation_equals( a );
}


BOOST_TEST_DONT_PRINT_LOG_VALUE( hmm::markov_order_n_converter<> );
BOOST_TEST_DONT_PRINT_LOG_VALUE( hmm::double_vec );
//BOOST_TEST_DONT_PRINT_LOG_VALUE( hmm::double_array );
BOOST_TEST_DONT_PRINT_LOG_VALUE( hmm::param_idx_vec );
BOOST_TEST_DONT_PRINT_LOG_VALUE( hmm::param_idx_array );

void
check_serialisation( unsigned markov_order )
{
	cout << "******* check_serialisation(): " << markov_order << endl;
	
	const std::string filename = "test_serialisation.hmm"; 
	model::ptr m = create_random_fully_connected_model( NUM_STATES, myrrh::int_power( ALPHABET_SIZE, markov_order + 1 ), markov_order );

	// save data to archive
    {
		std::ofstream ofs( filename.c_str() );
        boost::archive::text_oarchive oa( ofs );
        // write class instance to archive
        oa << m;
    	// archive and stream closed when destructors are called
    }

    // ... some time later restore the class instance to its orginal state
    model::ptr m_copy;
    {
        // create and open an archive for input
        std::ifstream ifs( filename.c_str() );
        boost::archive::text_iarchive ia( ifs );
        // read class state from archive
        ia >> m_copy;
        // archive and stream closed when destructors are called
    }
	boost::filesystem::remove( filename );

	BOOST_CHECK_EQUAL( m->converter, m_copy->converter );
	BOOST_CHECK_EQUAL( m->N, m_copy->N );
	BOOST_CHECK_EQUAL( m->M, m_copy->M );
	BOOST_CHECK_EQUAL( m->P, m_copy->P );
	//BOOST_CHECK_EQUAL( m->theta, m_copy->theta );
	BOOST_CHECK_EQUAL( m->pi_parameterisation, m_copy->pi_parameterisation );
	//BOOST_CHECK_EQUAL( m->pi_normaliser, m_copy->pi_normaliser );
	BOOST_CHECK_EQUAL( m->a_parameterisation, m_copy->a_parameterisation );
	//BOOST_CHECK_EQUAL( m->a_normaliser, m_copy->a_normaliser );
	BOOST_CHECK_EQUAL( m->b_parameterisation, m_copy->b_parameterisation );
	//BOOST_CHECK_EQUAL( m->b_normaliser, m_copy->b_normaliser );
}



test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite * test = BOOST_TEST_SUITE( "Serialisation test suite" );

	try
	{
		test->add( BOOST_TEST_CASE( &check_optional_serialisation ), 0 );
		test->add( BOOST_TEST_CASE( &check_multi_array_serialisation ), 0 );
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_serialisation, 
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

