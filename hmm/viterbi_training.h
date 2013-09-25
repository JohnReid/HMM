/**
@file

Copyright John Reid 2007, 2008

*/

#ifndef HMM_VITERBI_TRAINING_H_
#define HMM_VITERBI_TRAINING_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/training.h"
#include "hmm/viterbi.h"


namespace hmm {



/** Update estimators of the parameters by applying viterbi algorithm and using emissions. */
template< typename OutputSeq >
double
viterbi_training_update(
	const model & m,
	const OutputSeq & O,
	unsigned_vec & q_star,
	double_array & a_estimator,
	double_array & b_estimator,
	double_vec & pi_estimator )
{
	//calculate most likely state sequence
	const double log_P_star = viterbi(m, O, q_star);

	const unsigned T = boost::size( O );

	//is there anything to do?
	if( 0 == T ) return log_P_star; //no

	pi_estimator[ q_star[0] ] += 1.0;

	//add emissions and transitions
	boost::optional< unsigned > last_q;
	for( unsigned t = 0; T != t; ++t )
	{
		const unsigned q = q_star[t];

		//do we have a last q? i.e. has there been a transition
		if( last_q ) a_estimator[*last_q][q] += 1.0;

		//is the output known?
		if( O[t] != m.converter.order_n_size ) b_estimator[q][O[t]] += 1.0;

		//remember this q for next time
		last_q = q;
	}

	return log_P_star;
}

template< typename SeqRange >
double
viterbi_training_multi_sequence_iteration(
	model & m,
	const model::prior & prior,
	const SeqRange & sequences )
{
	using namespace boost;

	BOOST_STATIC_ASSERT( (
		is_same<
			typename range_value< SeqRange >::type,
			output_seq
		>::value ) );

	prior.check_matches_model(m);

	double_array a_estimator( prior.a );
	double_array b_estimator( prior.b );
	double_vec pi_estimator( prior.pi );

	double LL = 0.0;
	BOOST_FOREACH( const typename range_value< SeqRange >::type & seq, sequences )
	{
		unsigned_vec q_star;
		const double log_P_star = viterbi_training_update(
			m,
			seq,
			q_star,
			a_estimator,
			b_estimator,
			pi_estimator );
		HMM_CHECK_IS_A_NUMBER( log_P_star );
		LL += log_P_star;
	}

	impl::update_model(
		m,
		pi_estimator,
		a_estimator,
		b_estimator );

	return LL;
}

/** 
Takes multiple sequences and performs Viterbi training.
Returns the log likelihood (last but one LL - re-run forward() to get the true last LL) and the number of iterations
*/
template< 
	typename SeqRange,
	typename Callback
>
training_result_tuple
viterbi_training_multi_sequence(
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
		viterbi_training_multi_sequence_iteration< SeqRange >,
		max_iterations,
		tolerance,
		callback );
}



} //namespace hmm

#endif //HMM_VITERBI_TRAINING_H_

