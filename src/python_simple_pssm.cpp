/**
@file

Copyright John Reid 2007, 2012

*/

#ifdef _MSC_VER
#pragma warning( disable : 4172 )
#endif //_MSVC



#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"

#include <boost/python/detail/prefix.hpp>
#include "hmm/model.h"
#include "hmm/simple_pssm.h"
#include "hmm/python.h"

#include <myrrh/python/numpy.h>

//#include  <boost/iterator/pointee.hpp>


namespace hmm {
namespace pssm {
namespace impl {

using boost::python::object;
using boost::python::extract;
using boost::python::tuple;
using myrrh::python::convert_from_python;
using hmm::impl::double_vec_ptr;
using hmm::impl::output_seq_shared_ptr;
using hmm::impl::sequence_vec_shared_ptr;

typedef boost::shared_ptr< double_array > double_array_ptr;



template< typename CPtrType >
CPtrType
extract_or_convert( object py_obj )
{
	extract< CPtrType > e( py_obj );
	if( e.check() ) {
		return e();
	} else {
		CPtrType c_obj( new typename boost::pointee< CPtrType >::type );
		convert_from_python( py_obj, *c_obj );
		return c_obj;
	}
}



double_array_ptr
calculate_log_scores( object py_nucleo_dist )
{
	double_array_ptr nucleo_dist( extract_or_convert< double_array_ptr >( py_nucleo_dist ) );

	double_array_ptr log_scores( new double_array( boost::extents[boost::size(*nucleo_dist)][5] ) );
	calculate_log_scores_for_nmer( *nucleo_dist, *log_scores );

	return log_scores;
}


double_array_ptr
calculate_complementary_scores( object py_log_scores )
{
	double_array_ptr log_scores( extract_or_convert< double_array_ptr >( py_log_scores ) );

	double_array_ptr complementary_scores( new double_array( boost::extents[boost::size(*log_scores)][5] ) );
	hmm::pssm::calculate_complementary_scores( *log_scores, *complementary_scores );

	return complementary_scores;
}


double_vec_ptr
score_sequence( object py_pssm_scores, object py_seq )
{
	double_array_ptr pssm_scores( extract_or_convert< double_array_ptr >( py_pssm_scores ) );
	output_seq_shared_ptr seq( extract_or_convert< output_seq_shared_ptr >( py_seq ) );

	double_vec_ptr result( new double_vec() );

//#define TEST_HOW_QUICK_VEC_VEC_IS
#ifdef TEST_HOW_QUICK_VEC_VEC_IS

	std::vector< std::vector< double > > pssm_in_vec( boost::size( *pssm_scores ) );
	for( int i = 0; boost::size( *pssm_scores ) != i; ++i ) {
		pssm_in_vec[i].resize(5);
		for( int j = 0; 5 != j; ++j ) {
			pssm_in_vec[i][j] = (*pssm_scores)[i][j];
		}
	}

	apply_log_scores_to_sequence( pssm_in_vec, *seq, std::back_inserter( *result ), double_vec_ptr() );

#else //TEST_HOW_QUICK_VEC_VEC_IS

	apply_log_scores_to_sequence( *pssm_scores, *seq, std::back_inserter( *result ), double_vec_ptr() );

#endif //TEST_HOW_QUICK_VEC_VEC_IS

	return result;
}

struct find_max {
	unsigned position;
	unsigned max_position;
	double max_score;
	find_max() : position(0), max_position(-1), max_score( -std::numeric_limits< double >::max() ) { }
	void operator()( double score ) {
		if( score > max_score ) {
			max_position = position;
			max_score = score;
		}
		++position;
	}
};

template< typename Fn >
struct indirect_function_object {
	Fn & fn;
	indirect_function_object( Fn & fn ) : fn( fn ) { }
	template< typename T >
	void operator()( T t ) { return fn( t ); }
};

template< typename Fn >
indirect_function_object< Fn >
make_indirect_function_object( Fn & fn ) {
	return indirect_function_object< Fn >( fn );
}


boost::tuple< double, unsigned > //max score, position
max_score_in_sequence( object py_pssm_scores, object py_seq, object py_background_scores )
{
	double_array_ptr pssm_scores( extract_or_convert< double_array_ptr >( py_pssm_scores ) );
	output_seq_shared_ptr seq( extract_or_convert< output_seq_shared_ptr >( py_seq ) );
	double_vec_ptr background_scores( extract_or_convert< double_vec_ptr >( py_background_scores ) );

	double_vec_ptr result( new double_vec() );

	find_max max_finder;
	apply_log_scores_to_sequence( 
		*pssm_scores, 
		*seq, 
		boost::make_function_output_iterator( make_indirect_function_object( max_finder ) ),
		background_scores );

	return boost::make_tuple( max_finder.max_score, max_finder.max_position );
}



boost::python::list
max_scores_in_sequences( object py_pssm_scores, object py_seqs, object py_background_scores )
{
	double_array_ptr pssm_scores( extract_or_convert< double_array_ptr >( py_pssm_scores ) );
	sequence_vec_shared_ptr seqs( extract_or_convert< sequence_vec_shared_ptr >( py_seqs ) );

	boost::python::list result;

	unsigned i = 0;
	BOOST_FOREACH( const output_seq & seq, *seqs ) {
		double_vec_ptr background_scores;
		if( py_background_scores ) {
			background_scores = extract_or_convert< double_vec_ptr >( py_background_scores[i] );
		}
		find_max max_finder;
		apply_log_scores_to_sequence( 
			*pssm_scores, 
			seq, 
			boost::make_function_output_iterator( make_indirect_function_object( max_finder ) ),
			background_scores );
		result.append( max_finder.max_score );
		++i;
	}

	return result;
}



} //namespace impl
} //namespace pssm


HMM_EXPORT_FN_SPEC
void 
export_simple_pssm()
{
	using namespace boost::python;
	using namespace indexing;
	using namespace myrrh::python;
	using pssm::double_array;
	using pssm::impl::double_array_ptr;

	def( 
		"calculate_log_scores",
		pssm::impl::calculate_log_scores, 
		"Take an array of nucleotide distributions and return the log scores." 
	);

	def( 
		"calculate_complementary_scores",
		pssm::impl::calculate_complementary_scores, 
		"Take an array of log scores and return their reverse complement." 
	);

	def( 
		"max_score_in_sequence",
		pssm::impl::max_score_in_sequence, 
		"Takes a sequence and an array of PSSM log scores and returns a tuple of the best score and its location." 
	);

	def( 
		"max_scores_in_sequences",
		pssm::impl::max_scores_in_sequences, 
		"Takes a sequence of sequences and an array of PSSM log scores and returns a list of the best scores in each sequence." 
	);

	def( 
		"score_sequence",
		pssm::impl::score_sequence, 
		"Takes a sequence and an array of PSSM log scores and returns a sequence of the scores applied to the sequence." 
	);


	class_< 
		double_array, 
		double_array_ptr
	>( 
		"Array",
		"Holds an array of values.",
		no_init
	)
		.def(
			"as_array",
			convert_to_python< double_array >, 
			"Returns a numpy array with a copy of this array's data." 
		)
		;
}


} //namespace hmm
