/**
@file

Copyright John Reid 2007, 2008, 2012

*/

#ifndef HMM_BAUM_WELCH_H_
#define HMM_BAUM_WELCH_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/training.h"
#include "hmm/forward_backward.h"

#ifdef _MSC_VER
#pragma float_control(push)
//#pragma float_control(except, on)
//#pragma float_control(precise, on)
//#pragma fp_contract(on)
#pragma fenv_access(on)
#endif //_MSC_VER




namespace hmm {





/**
Updates the model given the sequence and the forward and backward probabilities in alpha and beta.
*/
template< typename OutputSeq >
void
baum_welch_updater(
	model & m,
	const OutputSeq & O,
	const double_array & alpha,
	const double_array & beta,
	const double_vec & c )
{
#ifdef HMM_VERIFY_ON
	const unsigned T = O.size();
#endif //HMM_VERIFY_ON
	const unsigned N = m.N;
	const unsigned M = m.M;

	HMM_VERIFY( c.size() == T );

	HMM_VERIFY( alpha.shape()[ 0 ] == T );
	HMM_VERIFY( alpha.shape()[ 1 ] == N );

	HMM_VERIFY( beta.shape()[ 0 ] == T );
	HMM_VERIFY( beta.shape()[ 1 ] == N );

	double_array a_estimator( boost::extents[ N ][ N ] );
	fill_multi_array( a_estimator, 0.0 );
	double_array b_estimator( boost::extents[ N ][ M ] );
	fill_multi_array( b_estimator, 0.0 );
	double_vec pi_estimator( N, 0.0 );

	impl::update_estimators(
		m,
		O,
		alpha,
		beta,
		c,
		pi_estimator,
		a_estimator,
		b_estimator );

	impl::update_model(
		m,
		pi_estimator,
		a_estimator,
		b_estimator );
}



/** 
Takes multiple sequences and performs one iteration of Baum-Welch.
Returns the log likelihood (before modification - re-run forward() to get new LL)
*/
template< typename SeqRange >
double
baum_welch_multi_sequence_iteration(
	model & m,
	const model::prior & prior,
	const SeqRange & sequences )
{
	using namespace boost;

	prior.check_matches_model(m);

	double_array a_estimator( prior.a );
	double_array b_estimator( prior.b );
	double_vec pi_estimator( prior.pi );

	double LL = 0.0;
	BOOST_FOREACH( const typename range_value< SeqRange >::type & seq, sequences )
	{
		double_array alpha;
		double_vec c;
		const double LL_for_seq = forward( 
			m,
			seq,
			alpha,
			c );
		HMM_CHECK_IS_A_NUMBER( LL_for_seq );
		LL += LL_for_seq;

		double_array beta;
		backward( 
			m,
			seq,
			c,
			beta );

#ifdef _DEBUG //verify calculations are being done correctly
		{
			const double LL_from_beta = LL_from_backward_probabilities( m, seq, c, beta );
			HMM_CHECK_IS_A_NUMBER( LL_from_beta );
			if( ! myrrh::is_close( LL_for_seq, LL_from_beta, 1e-6 ) )
			{
				throw std::logic_error( HMM_MAKE_STRING( LL_for_seq << " != " << LL_from_beta ) );
			}
			check_forward_backward_consistency( alpha, beta, c );
		}
#endif //_DEBUG

		impl::update_estimators(
			m,
			seq,
			alpha,
			beta,
			c,
			pi_estimator,
			a_estimator,
			b_estimator );
	}

	impl::update_model(
		m,
		pi_estimator,
		a_estimator,
		b_estimator );

	return LL;
}


/** 
Takes multiple sequences and performs Baum-Welch.
Returns the log likelihood (last but one LL - re-run forward() to get the true last LL) and the number of iterations
*/
template< 
	typename SeqRange,
	typename Callback
>
training_result_tuple
baum_welch_multi_sequence(
	model & m,
	const model::prior & prior,
	const SeqRange & sequences,
	unsigned max_iterations = 0,
	double tolerance = 1e-5,
	Callback callback = impl::null_callback() )
{
	return impl::train_multi_sequence(
		m,
		prior,
		sequences,
		baum_welch_multi_sequence_iteration< SeqRange >,
		max_iterations,
		tolerance,
		callback );
}




} //namespace hmm


#ifdef _MSC_VER
#pragma float_control(pop)
#endif //_MSC_VER

#endif //HMM_BAUM_WELCH_H_

