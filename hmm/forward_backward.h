/**
@file

Copyright John Reid 2007, 2008

*/

#ifndef HMM_FORWARD_BACKWARD_H_
#define HMM_FORWARD_BACKWARD_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/defs.h"
#include "hmm/model.h"

#include <myrrh/math.h>

#ifdef _MSC_VER
#pragma float_control(push)
//#pragma float_control(except, on)
//#pragma float_control(precise, on)
//#pragma fp_contract(on)
#pragma fenv_access(on)
#endif //_MSC_VER


namespace hmm {



template< typename Range >
double
scale_if_required(
	Range & x,
	double scaling_threshold )
{
	using boost::begin;
	using boost::end;
	using boost::range_value;

	const double sum = std::accumulate( begin( x ), end( x ), 0.0 );

	//do we want to scale?
	if( sum < scaling_threshold )
	{
		//yes so scale all alphas
		const double scale_factor = 1.0 / sum;
		BOOST_FOREACH( typename range_value< Range >::type & v, x ) v *= scale_factor;
		return scale_factor;
	}
	else
	{
		return 1.0;
	}
}

inline
double
log_scale_factors( const double_vec & c )
{
	double result = 0.0;
	BOOST_FOREACH( double x, c ) result += myrrh::safe_log( x );
	return result;
}


/** 
Forward procedure for HMMs implemented with scaling factors, c. 
Returns log likelihood. 
*/
template< typename OutputSeq >
double
forward( 
	const model & m,
	const OutputSeq & O,
	double_array & alpha,
	double_vec & c,
	double scaling_threshold = 1e-8 )
{
	const unsigned N = m.N;
	const unsigned T = boost::size( O );

	alpha.resize( boost::extents[ T ][ N ] );
	c.resize( T );

	if( T == 0 ) return 0.0;

        //more efficient if we calculate A and B ahead of time...
        double_array A;
        calculate_model_transition_matrix(m, A);
        double_array B;
        calculate_model_emission_matrix_with_unknown(m, B);

	//initialisation
	for( unsigned i = 0; N != i; ++i )
	{
		alpha[ 0 ][ i ] = m.pi( i ) * B[ i ][ O[ 0 ] ];
	}
	c[ 0 ] = 1.0;

	//more efficient if we just examine those states that can preceed any given state so calculate them here
	unsigned_vec_vec predecessors;
	calculate_predecessors( m, predecessors );	

	//induction
	for( unsigned t = 0; T - 1 != t; ++t )
	{
		for( unsigned j = 0; N != j; ++j )
		{
			double tmp = 0.0;
			BOOST_FOREACH( unsigned i, predecessors[ j ] )
			//for( unsigned i = 0; N != i; ++i ) 
			{
				tmp += alpha[ t ][ i ] * A[ i ][ j ];
			}
			alpha[ t + 1 ][ j ] = tmp * B[ j ][ O[ t + 1 ] ];
		}

		double_array::reference alpha_row = alpha[ t + 1 ];
		c[ t + 1 ] = 
			( T - 2 != t ) 
				? scale_if_required( alpha_row, scaling_threshold )
				: myrrh::scale( alpha_row ); //always scale the last one
	}

	//termination
	double result = 0.0;
	for( unsigned i = 0; N != i; ++i ) result += alpha[ T - 1 ][ i ];
	return myrrh::safe_log( result ) - log_scale_factors( c );
}


/**
Convenience function to get the log likelihood.
*/
template< typename OutputSeq >
double
forward_LL(
	const model & m,
	const OutputSeq & O )
{
	double_array alpha;
	double_vec c;
	return forward( m, O, alpha, c );
}


/**
The combined log likelihood of a range of sequences.
*/
template< typename SeqRange >
double
LL_multiple_sequences(
	const model & m,
	const SeqRange & sequences )
{
	using namespace boost;

	double LL = 0.0;
	double_array alpha;
	double_vec c;
	BOOST_FOREACH( const typename range_value< SeqRange >::type & seq, sequences ) LL += forward( m, seq, alpha, c );
	return LL;
}



/**
Backward procedure for HMMs implemented with scaling factors, c.
*/
template< typename OutputSeq >
void
backward(
	const model & m,
	const OutputSeq & O,
	const double_vec & c,
	double_array & beta )
{
	const unsigned N = m.N;
	const unsigned T = boost::size( O );

	beta.resize( boost::extents[ T ][ N ] );

	if( 0 == T ) return;

	//initialisation
	for( unsigned i = 0; N != i; ++i )
	{
		beta[ T - 1 ][ i ] = c[ T - 1 ];
	}

	//more efficient if we calculate A and B ahead of time...
	double_array A;
	calculate_model_transition_matrix(m, A);
	double_array B;
	calculate_model_emission_matrix_with_unknown(m, B);

	//more efficient if we just examine those states that can succeed any given state so calculate them here
	unsigned_vec_vec successors;
	calculate_successors( m, successors );	

	//induction
	for( unsigned t = T - 1; 0 != t; --t )
	{
		for( unsigned i = 0; N != i; ++i )
		{
			double tmp = 0.0;
			BOOST_FOREACH( unsigned j, successors[ i ] )
			//for( unsigned j = 0; N != j; ++j ) 
			{
				tmp += A[ i ][ j ] * B[ j ][ O[ t ] ] * beta[ t ][ j ];
			}
			beta[ t - 1 ][ i ] = tmp * c[ t - 1 ];
		}
	}
}





/**
Calculate the LL from the backward probabilities.

Can use this to check forward and backward algorithms agree.
*/
template< typename OutputSeq >
double
LL_from_backward_probabilities(
	const model & m,
	const OutputSeq & O,
	const double_vec & c,
	double_array & beta )
{
	if( 0 == c.size() ) return 0.0;

	const unsigned N = m.N;
	double result = 0.0;
	for( unsigned i = 0; N != i; ++i ) result += m.b( i, O[ 0 ] ) * beta[ 0 ][ i ] * m.pi( i );
	return myrrh::safe_log( result ) - log_scale_factors( c );
}


/**
Check the forward and backward probabilities for consistency.
*/
inline
void
check_forward_backward_consistency(
	const double_array & alpha,
	const double_array & beta,
	const double_vec & c )
{
	HMM_VERIFY( alpha.shape()[ 0 ] == beta.shape()[ 0 ] );
	HMM_VERIFY( alpha.shape()[ 1 ] == beta.shape()[ 1 ] );

	const unsigned T = alpha.shape()[ 0 ];
	const unsigned N = alpha.shape()[ 1 ];

	//use alpha and beta at each time point to calculate p(observed seq) i.e likelihood
	double_vec p( T, 0.0 );
	for( unsigned t = 0; T != t; ++t )
	{
		for( unsigned i = 0; N != i; ++i )
		{
			p[ t ] += alpha[ t ][ i ] * beta[ t ][ i ] / c[ t ];
		}

		//check all calculations are very similar
		if( 0 != t && ! myrrh::is_close( p[ t ], p[ t - 1 ], 1e-6 ) ) 
			throw std::logic_error( "p(obs) not calculated correctly by forward-backward procedures" );	
	}
}

} //namespace hmm


#ifdef _MSC_VER
#pragma float_control(pop)
#endif //_MSC_VER


#endif //HMM_FORWARD_BACKWARD_H_

