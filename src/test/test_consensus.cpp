

#include <hmm/consensus.h>
#include "test.h"
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


using namespace hmm;
using namespace hmm::consensus;

std::vector< unsigned > alphabet_sizes = list_of
	( 3 )
	( 4 )
	( 7 )
	( 9 )
	;

void
check_consensus( unsigned alphabet_size )
{
	cout << "******* check_consensus(): " << alphabet_size << endl;

	codes<> consensus_codes( alphabet_size );

	//check the sizes have been calculated correctly
	BOOST_CHECK_EQUAL( consensus_codes.alphabet_size(), alphabet_size );
	BOOST_CHECK_EQUAL( consensus_codes.num_codes(), codes<>::code_t( myrrh::int_power(2, alphabet_size) ) );

	//check we only have the one char encoded in simple codes
	for( unsigned i = 0; alphabet_size != i; ++i )
	{
		const unsigned code = consensus_codes.code_for_char(i);
		for( unsigned j = 0; alphabet_size != j; ++j )
		{
			const bool matches = consensus_codes.matches( j, code );
			BOOST_CHECK_EQUAL( i == j, matches );
		}
		//check the original char is in all the possible codes for this character
		std::vector< unsigned > calculated_codes;
		consensus_codes.codes_for_char( i, std::back_inserter(calculated_codes) );
		BOOST_CHECK( calculated_codes.end() != std::find( calculated_codes.begin(), calculated_codes.end(), code ) );
	}

	//check for more complicated codes that we have the correct chars encoded
	std::set< unsigned > chars = list_of(1)(0);
	const unsigned code = consensus_codes.code_for_chars( chars );
	BOOST_CHECK_EQUAL( consensus_codes.num_chars_in_code( code ), chars.size() );
	std::vector< unsigned > calculated_chars;
	consensus_codes.chars_for_code( code, std::back_inserter(calculated_chars) );
	for( unsigned j = 0; alphabet_size != j; ++j )
	{
		const bool in_chars = std::find( chars.begin(), chars.end(), j ) != chars.end();
		const bool in_calc = std::find( calculated_chars.begin(), calculated_chars.end(), j ) != calculated_chars.end();
		BOOST_CHECK_EQUAL( in_chars, in_calc );
	}
}

void
check_consensus_word_counting()
{
	cout << "******* check_consensus_word_counting()" << endl;
	codes<> consensus_codes( 4 );
	const std::vector<unsigned> seq = list_of(0)(0)(0)(0)(0)(0)(0)(0)(0)(0)(0)(0)(0)(0)(0);
	boost::multi_array<unsigned, 4> counts;
	consensus_codes.count_consensus_L_words( counts, seq.begin(), seq.end() );
}


test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	const char * test_suite_name = "consensus test suite";
	test_suite * test = BOOST_TEST_SUITE( test_suite_name );

	try
	{
		test->add( 
			BOOST_PARAM_TEST_CASE( 
				&check_consensus, 
				alphabet_sizes.begin(),
				alphabet_sizes.end() ),
			0 );
		test->add(BOOST_TEST_CASE( &check_consensus_word_counting ), 0 );
	}
	catch (const std::exception & e)
	{
		cerr << "Registering " << test_suite_name << " - exception: " << e.what() << endl;
	}

    return test; 
}

