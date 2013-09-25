/**
@file

Copyright John Reid 2007

*/

#ifndef HMM_SAMPLE_H_
#define HMM_SAMPLE_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "hmm/defs.h"

namespace hmm {


inline
void
sample(
	const model & m,
	unsigned T,
	unsigned_vec & states,
	output_seq & output )
{
	double_array A;
	calculate_model_transition_matrix( m, A );

	double_array B;
	calculate_model_emission_matrix( m, B );

	double_vec pi;
	calculate_model_initial_distribution( m, pi );

	output.resize( T );
	states.resize( T );

	//choose a uniformly random order n-1 to use as initial previous obs
	unsigned previous_order_n_minus_1 = gsl_rng_uniform_int( myrrh::gsl::get_rng(), m.converter.order_n_minus_1_size );

	//for each base to sample
	for( unsigned t = 0; T != t; ++t )
	{
		states[ t ] = 
			0 == t 
				? myrrh::gsl::one_sample_from_multinomial( pi )
				: myrrh::gsl::one_sample_from_multinomial( A[ states[ t - 1 ] ] );
		typedef double_array::value_type::index_range range;
		unsigned order_0_output = myrrh::gsl::one_sample_from_multinomial( 
			B[states[t]][ 
				boost::indices[
					range(
						previous_order_n_minus_1 * m.converter.order_0_size, 
						(previous_order_n_minus_1 + 1) * m.converter.order_0_size) ] ] );
		HMM_VERIFY( order_0_output < m.converter.order_0_size );
		output[t] = order_0_output;
		//output[t] = previous_order_n_minus_1 * m.converter.order_0_size + order_0_output;
		//HMM_VERIFY( output[t] < m.converter.order_n_size );

		//the last order n-1 obs of this observation will be the previous order n-1 obs for the next one
		previous_order_n_minus_1 = m.converter.get_last_order_n_minus_1( output[t] );
	}
}


} //namespace hmm


#endif //HMM_SAMPLE_H_

