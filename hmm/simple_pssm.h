/**
@file

Copyright John Reid 2007, 2013

*/

#ifndef HMM_SIMPLE_PSSM_H_
#define HMM_SIMPLE_PSSM_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "hmm/defs.h"

namespace hmm {
namespace pssm {

using myrrh::double_array;


/**
Applies the scores to the sequence beginning at seq_begin.
*/
template< typename SeqIt, typename ScoresRange >
double
apply_log_scores_to_nmer( const ScoresRange & scores, SeqIt seq_begin )
{
	using namespace boost;

	double result = 0.0;
	typedef typename range_value< const ScoresRange >::type value_t;
	BOOST_FOREACH( const value_t & s, scores ) {
		result += s[ *seq_begin ];
		++seq_begin;
	}
	return result;
}

/**
Applies the scores to the sequence.
*/
template< typename ScoresRange, typename SeqRange, typename OutputIt, typename BgScoresPtr >
void
apply_log_scores_to_sequence( 
	const ScoresRange & scores, 
	const SeqRange & sequence, 
	OutputIt output_it,
	BgScoresPtr bg_scores )
{
	using namespace boost;
	typedef typename boost::pointee< BgScoresPtr >::type bg_scores_container;
	typedef typename bg_scores_container::iterator bg_scores_it;
	typedef typename range_iterator< const SeqRange >::type seq_it;
	typename range_iterator< const SeqRange >::type seq_begin = boost::begin( sequence );
	typename range_iterator< const SeqRange >::type seq_end = boost::end( sequence ) - size( scores ) + 1;
	if( seq_end <= seq_begin ) return; //nothing to do, PSSM is too long
	bg_scores_it bg_score;
	if( bg_scores )
		bg_score = bg_scores->begin();
	while( seq_begin != seq_end ) {
		double score = apply_log_scores_to_nmer( scores, seq_begin );
		if( bg_scores )
			score -= *bg_score;
		*output_it = score;
		++output_it;
		++seq_begin;
		if( bg_scores )
			++bg_score;
	}
}



/**
Takes a nucleotide distribution over an n-mer and creates log scores.
*/
template< typename DistRange, typename Scores >
void
calculate_log_scores_for_nmer( const DistRange & dist, Scores & scores )
{
	using namespace boost;
	if( size(scores) != size(dist) ) throw std::runtime_error( "dist and scores must have the same length" );
	for( int i = 0; size(scores) != size_t( i ); ++i ) {
		if( 4 != size(dist[i]) ) throw std::runtime_error( "dist must have elements with length 4" );
		if( 5 != size(scores[i]) ) throw std::runtime_error( "scores must have elements with length 5" );
		for( unsigned j = 0; 4 != j; ++j ) scores[i][j] = log( dist[i][j] );
		scores[i][4] = log( .25 );
	}
}

/**
Takes log scores over an n-mer and returns scores for the complementary strand.
*/
template< typename Scores >
void
calculate_complementary_scores( const Scores & scores, Scores & complementary_scores )
{
	using namespace boost;
	using boost::size;

	if( size( scores ) != size( complementary_scores ) ) throw std::runtime_error( "Scores should have entries of length 5." );

	const unsigned len = size( scores );
	for( unsigned i = 0; len != i; ++i ) 
	{
		const unsigned i_dash = len - i - 1;
		if( size( scores[i] ) != 5 ) throw std::runtime_error( "Scores should have entries of length 5." );
		if( size( complementary_scores[i_dash] ) != 5 ) throw std::runtime_error( "Complementary scores should have entries of length 5." );

		for( unsigned j = 0; 4 != j; ++j ) 
		{
			const unsigned j_dash = 3 - j;
			complementary_scores[i_dash][j_dash] = scores[i][j];
		}
		complementary_scores[i_dash][4] = scores[i][4];
	}
}



} //namespace pssm
} //namespace hmm


#endif //HMM_SIMPLE_PSSM_H_

