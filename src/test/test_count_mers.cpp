#ifdef _MSC_VER
# pragma warning( disable : 4172 )
#endif //_MSC_VER

#include <hmm/count_mers.h>
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

struct output
{
	std::vector< output_seq > prefixes;
	std::vector< unsigned > counts;

	template< typename Prefix >
	void operator()( const Prefix & prefix, const mer_count & count )
	{
		prefixes.push_back( prefix );
		counts.push_back( count.total );
	}
};

void
check_count_mers()
{
	cout << "******* check_count_mers()" << endl;

	typedef std::vector< output_seq > SeqRange;
	const SeqRange seqs = list_of<output_seq>
		(
			list_of(0)(1)(2)(3)(4)
		)
		(
			list_of(0)(1)(2)(3)(4)
		)
		;

	for( unsigned i = 3; 6 != i; ++i )
	{
		output o;
		count_mers< mer_count >( seqs, i, o );
		BOOST_CHECK_EQUAL( o.prefixes.size(), 6-i );
		BOOST_FOREACH( unsigned count, o.counts )
		{
			BOOST_CHECK_EQUAL( count, unsigned( 2 ) );
		}
	}
}





test_suite *
init_unit_test_suite( int argc, char * argv [] )
{
	test_suite * test = BOOST_TEST_SUITE( "HMM count_mers test suite" );

	try
	{
		test->add( BOOST_TEST_CASE( &check_count_mers ), 0 );
	}
	catch (const std::exception & e)
	{
		cerr << "Registering tests - exception: " << e.what() << endl;
	}

    return test; 
}

