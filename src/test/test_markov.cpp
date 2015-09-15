/**
@file

Copyright John Reid 2015

*/

#include <hmm/defs.h>
#include <hmm/markov.h>
#include "test.h"

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/timer.hpp>

using namespace hmm;

typedef seqan::StringSet< seqan::CharString > id_set_t;
typedef seqan::String< seqan::Dna5 > string_t;
typedef seqan::StringSet< string_t > string_set_t;


BOOST_AUTO_TEST_CASE( markov_model_init_test )
{
    typedef complete_markov_model< 3, 4, double > model_t;

    {
        model_t model;
        BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][0], 0. );

        // initialise counts with pseudo-count of 6.
        add_pseudo_counts( model, double( 6 ) );
        BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][0], 6. );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][0], 6. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][0], 6. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[0], 6. );

        // add pseudo-counts of 3.
        add_pseudo_counts( model, double( 3 ) );
        BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][0], 9. );
    }

    {
        model_t model;
        // add counts from a sequence
        vector< size_t > seq = list_of(0)(1)(2)(3)(0);
        add_counts_to_complete( model, seq.begin(), seq.end() );
        BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][3], 1. );
        BOOST_CHECK_EQUAL( model.mm.storage[3][1][2][3], 0. );
        BOOST_CHECK_EQUAL( model.mm.storage[0][0][2][3], 0. );
        BOOST_CHECK_EQUAL( model.mm.storage[0][1][1][3], 0. );
        BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][2], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][3], 1. );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[3][2][3], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][3][3], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][0], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][3], 1. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[1][3], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][0], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[0], 2. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[1], 1. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[2], 1. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.lower_orders.mm.storage[3], 1. );

        // normalise counts
        normalise_counts( model );
        BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][3], 1. );
        BOOST_CHECK_MESSAGE( std::isnan( model.mm.storage[3][1][2][3] ), "Should be NaN" );
        BOOST_CHECK_MESSAGE( std::isnan( model.mm.storage[0][0][2][3] ), "Should be NaN" );
        BOOST_CHECK_MESSAGE( std::isnan( model.mm.storage[0][1][1][3] ), "Should be NaN" );
        BOOST_CHECK_EQUAL( model.mm.storage[0][1][2][2], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][3], 1. );
        BOOST_CHECK_MESSAGE( std::isnan( model.lower_orders.mm.storage[3][2][3] ), "Should be NaN" );
        BOOST_CHECK_MESSAGE( std::isnan( model.lower_orders.mm.storage[1][3][3] ), "Should be NaN" );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][0], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][3], 1. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[1][3], 0. );
        BOOST_CHECK_EQUAL( model.lower_orders.lower_orders.mm.storage[2][0], 0. );
        BOOST_CHECK_CLOSE( model.lower_orders.lower_orders.lower_orders.mm.storage[0], 2./5., 0.0001 );

        // convert to scale
        convert_to_scale( model );
        BOOST_CHECK_EQUAL( model.lower_orders.mm.storage[1][2][3], 0. );
        BOOST_CHECK_CLOSE( model.lower_orders.lower_orders.lower_orders.mm.storage[0], std::log( 2./5. ), 0.0001 );
    }
}



BOOST_AUTO_TEST_CASE( markov_model_evaluation_test )
{
    // create and initialise model.
    typedef complete_markov_model< 3, 4, double > model_t;
    model_t model;
    vector< size_t > seq = list_of(0)(1)(2)(3)(0);
    initialise_from_sequence( model, 1., seq.begin(), seq.end() );

    {
        vector< size_t > eval_seq = list_of(0);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
    }

    {
        vector< size_t > eval_seq = list_of(0)(1);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
    }

    {
        vector< size_t > eval_seq = list_of(0)(1)(2);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
    }

    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
    }

    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
    }

    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 1./4. ), 0.0001 );
    }

    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate( eval_seq.begin(), eval_seq.end(), boost::make_function_output_iterator( make_push_back_cumulative( probs ) ) );
        BOOST_CHECK_EQUAL( eval_seq.size(), probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 3./9. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[6] - probs[5], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[7] - probs[6], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is not large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 1,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 1, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[6] - probs[5], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is not large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 2,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 2, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[5] - probs[4], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is just right size
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 3,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 3, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[4] - probs[3], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 4,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 4, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 2./5. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[2] - probs[1], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[3] - probs[2], std::log( 2./5. ), 0.0001 );
    }

    // test evaluation within context, when context is large enough
    {
        vector< size_t > eval_seq = list_of(0)(1)(2)(3)(0)(1)(2)(3);
        vector< double > probs;
        model.evaluate_in_context(
            eval_seq.begin(),
            eval_seq.begin() + 6,
            eval_seq.end(),
            boost::make_function_output_iterator( make_push_back_cumulative( probs ) )
        );
        BOOST_CHECK_EQUAL( eval_seq.size() - 6, probs.size() );
        BOOST_CHECK_CLOSE( probs[0], std::log( 1./4. ), 0.0001 );
        BOOST_CHECK_CLOSE( probs[1] - probs[0], std::log( 2./5. ), 0.0001 );
    }
}


/**
 * Read a FASTA file into a string set.
 */
void
read_fasta(
    const char * seqPath,
    id_set_t & ids,
    string_set_t & X
) {
    // Open the sequence file
    seqan::SeqFileIn seqFile;
    if( ! seqan::open(seqFile, seqPath) ) {
        throw std::logic_error("Could not open file.");
    }
    // Read the sequences from a file into X.
    std::cout << "Reading sequences from " << seqPath << " ...\n";
    seqan::readRecords(ids, X, seqFile);
}

void
print_seq_info(
    const char * filename,
    string_set_t & sequences
) {
    using namespace seqan;
    Size< string_t >::Type max_length = 0;
    Size< string_t >::Type total_bases = 0;
    for( Size< string_set_t >::Type i = 0; length( sequences ) != i; ++i ) {
        Size< string_t >::Type len = length( value( sequences, i ) );
        total_bases += len;
        max_length = std::max( max_length, len );
    }
    std::cout << "Read " << total_bases << " base pairs from " << length( sequences ) << " sequences from " << filename << "\n";
    std::cout << "Longest sequence has " << max_length << " base pairs\n";
}

const size_t ALPHABET_SIZE = 5;

template< size_t Order, size_t AlphabetSize = ALPHABET_SIZE >
void
check_markov()
{
    std::cout << "**************** Order: " << Order << " *****************\n";
    typedef complete_markov_model< Order, AlphabetSize > mm_t;

    // Read FASTA
    id_set_t ids;
    string_set_t X;
    const char * fasta = "etc/fasta/dm01r.fasta";
    read_fasta(fasta, ids, X);
    print_seq_info(fasta, X);

    // Build Markov model
    using namespace seqan;
    typedef Iterator< string_set_t, Standard >::Type TIterator;
    mm_t mm;
    add_pseudo_counts( mm, 1.0 );
    for (TIterator it = begin(X, Standard()); it != end(X, Standard()); ++it) {
        add_counts_to_complete( mm, begin( *it ), end( *it ) );
    }
    normalise_counts( mm );
    convert_to_scale( mm );

    // Evaluate the model on the same FASTA file
    std::vector< double > lp;
    mm.evaluate( begin( X[0] ), end( X[0] ), std::back_inserter( lp ) );
    for( unsigned i = 0; 6 != i; ++i ) {
        std::cout << lp[i] << "\n";
    }

    // save data to archive
    const char * filename = "serialization.test";
    {
        std::ofstream ofs( filename );
        boost::archive::text_oarchive( ofs ) << mm;
    }

    // ... some time later restore the class instance to its orginal state
    mm_t newmm;
    {
        std::ifstream ifs( filename );
        boost::archive::text_iarchive( ifs ) >> newmm;
    }

    // Evaluate the new model on the same FASTA file
    std::vector< double > newlp;
    newmm.evaluate( begin( X[0] ), end( X[0] ), std::back_inserter( newlp ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( lp.begin(), lp.end(),
                                   newlp.begin(), newlp.end() );
}


bool
init_unit_test()
{
    test_suite * test = &boost::unit_test::framework::master_test_suite();
    test->p_name.value = "Markov model test suite";

    try
    {
#define MAX_ORDER 4
#define ADD_CHECK(z, n, check) test->add(BOOST_TEST_CASE((&check< n, ALPHABET_SIZE >)), 0);
        BOOST_PP_REPEAT(5, ADD_CHECK, check_markov)
    }
    catch (const std::exception & e)
    {
        cerr << "Registering tests - exception: " << e.what() << endl;
        throw  boost::unit_test::framework::setup_error( e.what() );
    }
    return test;
}

int main(int argc, char* argv[])
{
    return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}

