/**
@file

Copyright John Reid 2007, 2008, 2012, 2013

*/

#ifndef HMM_TRAINING_H_
#define HMM_TRAINING_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/defs.h"
#include "hmm/model.h"

#include <myrrh/math.h>

#ifdef _DEBUG
# define HMM_CHECK_IS_A_NUMBER( x ) if( ! MYRRH_ISFINITE( x ) ) throw std::logic_error( "not finite!" )
#else //_DEBUG
# define HMM_CHECK_IS_A_NUMBER( x ) 
#endif //_DEBUG


namespace hmm {

typedef boost::tuple< double, unsigned > training_result_tuple;

inline
void
calculate_gamma(
	const double_array & alpha,
	const double_array & beta,
	double_array & gamma )
{
	const unsigned T = alpha.shape()[ 0 ];
	HMM_VERIFY( T == beta.shape()[ 0 ] );

	const unsigned N = alpha.shape()[ 1 ];
	HMM_VERIFY( N == beta.shape()[ 1 ] );

	gamma.resize( boost::extents[ T ][ N ] );

	for( unsigned t = 0; T != t; ++t )
	{
		double normalisation_factor = 0.0;
		for( unsigned i = 0; N != i; ++i ) normalisation_factor += alpha[ t ][ i ] * beta[ t ][ i ];
		for( unsigned i = 0; N != i; ++i ) gamma[ t ][ i ] = alpha[ t ][ i ] * beta[ t ][ i ] / normalisation_factor;
	}
}




namespace impl {

struct null_callback { bool operator()( double LL ) const { return true; } };

template< typename OutputSeq >
void
update_estimators(
	const model & m,
	const OutputSeq & O,
	const double_array & alpha,
	const double_array & beta,
	const double_vec & c,
	double_vec & pi_estimator,
	double_array & a_estimator,
	double_array & b_estimator )
{
	const unsigned N = m.N;
	const unsigned M = m.M;

	const unsigned T = boost::size( O );

	//update the pi estimator
	if( T > 0 )
	{
		double_vec pi_tmp( N );
		for( unsigned i = 0; N != i; ++i ) pi_tmp[ i ] = alpha[ 0 ][ i ] * beta[ 0 ][ i ];
		const double pi_sum = std::accumulate( pi_tmp.begin(), pi_tmp.end(), 0.0 );
		for( unsigned i = 0; N != i; ++i ) pi_estimator[ i ] += pi_tmp[ i ] / pi_sum;
	}

	//more efficient if we just examine those states that can succeed any given state so calculate them here
	unsigned_vec_vec successors;
	calculate_successors( m, successors );	

	//more efficient if we calculate A and B ahead of time...
	double_array A;
	calculate_model_transition_matrix(m, A);
	double_array B;
	calculate_model_emission_matrix_with_unknown(m, B);

	//for each time point
	for( unsigned t = 0; T != t; ++t )
	{
		for( unsigned i = 0; N != i; ++i )
		{
			const double p_forward_and_backward = alpha[ t ][ i ] * beta[ t ][ i ] / c[ t ];
			if( M != O[ t ] ) b_estimator[ i ][ O[ t ] ] += p_forward_and_backward;

			//can't update 'a' for last time point
			if( T - 1 != t )
			{
				//for( unsigned j = 0; N != j; ++j )
				//const unsigned_vec & next_states = successors[i];
				//const int num_successors = next_states.size();
				//for( int s = 0; num_successors != s; ++s )
				BOOST_FOREACH( unsigned j, successors[ i ] )
				{
					//const unsigned j = next_states[ s ];
					const double expected_transition_i_to_j =
						alpha[ t ][ i ]
						* A[ i ][ j ]
						* B[ j ][ O[ t + 1 ] ]
						* beta[ t + 1 ][ j ];
					a_estimator[ i ][ j ] += expected_transition_i_to_j;
				}
			}
		}
	}
}



template<
	typename ParamIdxRange,
	typename EstimateRange
>
void
update_theta_estimates(
	const ParamIdxRange & parameterisation,
	const EstimateRange & estimator,
	double_vec & theta_numerator,
	double_vec & theta_denominator )
{
	using boost::begin;
	using boost::end;
    using boost::range_iterator;

	const double sum = std::accumulate( boost::begin( estimator ), boost::end( estimator ), 0.0 );
    typename range_iterator< const ParamIdxRange >::type p = boost::begin( parameterisation );
    typename range_iterator< const EstimateRange >::type e = boost::begin( estimator );
    for( ; boost::end( parameterisation ) != p; ++p, ++e )
	{
		if( ! *p ) continue; //no updating if nothing to update...
		HMM_CHECK_IS_A_NUMBER( *e );
		theta_numerator[ **p ] += *e;
		theta_denominator[ **p ] += sum;
	}
}

inline
void
update_parameter( double & parameter, double numerator, double denominator )
{
	if( denominator != 0.0 ) parameter = numerator / denominator;
}

inline
void
update_model(
	model & m,
	const double_vec & pi_estimator,
	const double_array & a_estimator,
	const double_array & b_estimator )
{
	double_vec theta_numerator( m.theta.size(), 0.0 );
	double_vec theta_denominator( m.theta.size(), 0.0 );

	for( unsigned i = 0; m.N != i; ++i )
	{
		update_theta_estimates(
			m.a_parameterisation[i],
			a_estimator[i],
			theta_numerator,
			theta_denominator );

		typedef param_idx_array::value_type::const_iterator it;
		it begin = m.b_parameterisation[i].begin();
		it end = begin + m.converter.order_0_size;

		typedef double_array::value_type::const_iterator b_it;
		b_it b_begin = b_estimator[i].begin();
		b_it b_end = b_begin + m.converter.order_0_size;

		while( begin != m.b_parameterisation[i].end() )
		{
			boost::iterator_range< it > range( begin, end );
			boost::iterator_range< b_it > b_range( b_begin, b_end );
			update_theta_estimates(
				range,
				b_range,
				theta_numerator,
				theta_denominator );
			begin = end;
			end += m.converter.order_0_size;

			b_begin = b_end;
			b_end += m.converter.order_0_size;
		}
	}

	update_theta_estimates(
		m.pi_parameterisation,
		pi_estimator,
		theta_numerator,
		theta_denominator );

	for( unsigned p = 0; m.P != p; ++p ) update_parameter( m.theta[p], theta_numerator[p], theta_denominator[p] );

	normalise_model_parameters( m );

#ifdef _DEBUG
	check_consistent( m );
#endif //_DEBUG
}



/** 
Takes multiple sequences and trains using supplied algorithm.
Returns the log likelihood (last but one LL - re-run forward() to get the true last LL) and the number of iterations
*/
template< 
	typename SeqRange,
	typename Algorithm,
	typename Callback
>
training_result_tuple
train_multi_sequence(
	model & m,
	const model::prior & prior,
	const SeqRange & sequences,
	Algorithm algorithm,
	unsigned max_iterations = 0,
	double tolerance = 1e-5,
	Callback callback = impl::null_callback(),
	bool check_LL_improves = false )
{
	double LL = - std::numeric_limits< double >::max();
	double LL_old, LL_improvement;
	unsigned iteration = 0;
	do
	{
		LL_old = LL;
		LL = algorithm( m, prior, sequences );
		HMM_CHECK_IS_A_NUMBER( LL );
		if( ++iteration == max_iterations ) break;
		if( ! callback( LL ) ) break;
		LL_improvement = LL - LL_old;
#ifdef HMM_VERIFY_ON
		if( check_LL_improves ) HMM_VERIFY( LL_improvement >= 0 || fabs( LL_improvement ) < 1e-8 );
#endif //HMM_VERIFY_ON
	}
	while( fabs( LL_improvement ) > tolerance );

	return boost::make_tuple( LL, iteration );
}



} //namespace impl




} //namespace hmm

#endif //HMM_TRAINING_H_

