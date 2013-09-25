/**
@file

Copyright John Reid 2007, 2008, 2012

*/

#ifdef _MSC_VER
# pragma warning( disable : 4172 )
#endif //_MSC_VER

#include "test.h"

#include <hmm/model.h>
#include <hmm/baum_welch.h>
#include <hmm/pssm_model.h>
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



const unsigned_vec num_states = list_of
	(  10 )
	(   1 )
	(   2 )
	(   3 )
	(   5 )
//	(  25 )
	;

const unsigned_vec pssm_lengths = list_of
	(  10 )
	(   1 )
	(   2 )
	(   3 )
	(   5 )
//	(  30 )
//	(  50 )
	;

#define ALPHABET_SIZE 4
#define NUM_SEQUENCES 3
#define SEQ_LENGTH 100
#define MIN_MARKOV_ORDER 2
#define MAX_MARKOV_ORDER 2
#define MAX_ITERATIONS 100
#define TOLERANCE 1e-6


void generate_sequence( unsigned length, unsigned_vec & seq, unsigned markov_order = 0 )
{
	seq.clear();
	if( 0 == markov_order )
	{
		for( unsigned i = 0; length != i; ++i )
		{
			seq.push_back( gsl_rng_uniform_int( myrrh::gsl::get_rng(), ALPHABET_SIZE ) );
		}
	}
	else
	{
		unsigned_vec order_0_seq;
		generate_sequence( length, order_0_seq, 0 );
		markov_order_n_converter< unsigned >( 
			ALPHABET_SIZE, 
			markov_order ).make_markov_order_n( 
				order_0_seq.begin(), 
				order_0_seq.end(), 
				std::back_inserter( seq ) );
	}
}

struct check_ll
{
	boost::optional< double > last_LL;
	bool operator()( double LL )
	{
		static const double tol = 1e-10;

		BOOST_CHECK( MYRRH_ISFINITE( LL ) );

		if( last_LL && *last_LL - LL >= tol ) 
			std::cout << *last_LL << "," << LL << "," << *last_LL - LL << "\n";

		BOOST_CHECK( 
			( ! last_LL ) 
			|| 
			( *last_LL - LL < tol ) );

		last_LL = LL;

		return true;
	}
};

void generate_sequences( std::vector< unsigned_vec > & seqs, unsigned markov_order = 0 )
{
	seqs.resize( NUM_SEQUENCES );
	BOOST_FOREACH( unsigned_vec & seq, seqs ) generate_sequence( SEQ_LENGTH, seq );
}

void check_model( model::ptr m )
{
	//generate sequences
	std::vector< unsigned_vec > seqs;
	generate_sequences( seqs, m->converter.order );

	model::prior prior( m->N, m->M );

	baum_welch_multi_sequence(
		*m,
		prior,
		seqs,
		MAX_ITERATIONS,
		TOLERANCE,
		check_ll() );
}

void
check_hmm_log_likelihood_increases_with_pssm_model( unsigned K )
{
	cout << "******* check_hmm_log_likelihood_increases_with_pssm_model(): " << K << endl;

	//seed the rng
	gsl_rng_set( myrrh::gsl::get_rng(), K );

	check_model( build_random_pssm_model( 0.01, K ) );
}

template< unsigned markov_order >
void
check_hmm_log_likelihood_increases( unsigned n )
{
	cout << "******* check_hmm_log_likelihood_increases(): states=" << n << ", order=" << markov_order << endl;

	//seed the rng
	gsl_rng_set( myrrh::gsl::get_rng(), n * (markov_order + 1) );

	check_model( create_random_fully_connected_model( n, myrrh::int_power( ALPHABET_SIZE, markov_order + 1 ), markov_order ) );
}





test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite * test = BOOST_TEST_SUITE( "HMM LL test suite" );

	try
	{
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_hmm_log_likelihood_increases< 0 >, 
				num_states.begin(),
				num_states.end() ),
			0 );
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_hmm_log_likelihood_increases< 1 >, 
				num_states.begin(),
				num_states.end() ),
			0 );
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_hmm_log_likelihood_increases< 2 >, 
				num_states.begin(),
				num_states.end() ),
			0 );
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_hmm_log_likelihood_increases< 3 >, 
				num_states.begin(),
				num_states.end() ),
			0 );

		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_hmm_log_likelihood_increases_with_pssm_model, 
				pssm_lengths.begin(),
				pssm_lengths.end() ),
			0 );
	}
	catch (const std::exception & e)
	{
		cerr << "Registering tests - exception: " << e.what() << endl;
	}

    return test; 
}

