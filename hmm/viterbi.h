/**
@file

Copyright John Reid 2007, 2008

*/

#ifndef HMM_VITERBI_H_
#define HMM_VITERBI_H_

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







void
calculate_log_a(
	const model & m,
	double_array & log_a )
{
	const unsigned N = m.N;
	log_a.resize( boost::extents[ N ][ N ] );
	for( unsigned i = 0; m.N != i; ++i )
	{
		for( unsigned j = 0; N != j; ++j )
		{
			log_a[ i ][ j ] = myrrh::safe_log( m.a( i, j ) );
		}
	}
}

void
calculate_log_b(
	const model & m,
	double_array & log_b )
{
	const unsigned N = m.N;
	const unsigned M = m.M;
	log_b.resize( boost::extents[ N ][ M + 1 ] ); //M+1 for unknown observation
	for( unsigned i = 0; m.N != i; ++i )
	{
		for( unsigned k = 0; M + 1 != k; ++k )
		{
			log_b[ i ][ k ] = myrrh::safe_log( m.b( i, k ) );
		}
	}
}

/** The viterbi algorithm - returns log( P_star ). */
template< typename OutputSeq >
double
viterbi(
	const model & m,
	const OutputSeq & O,
	unsigned_vec & q_star )
{
	const unsigned N = m.N;
	const unsigned T = boost::size( O );

	q_star.resize( T );

	if( T == 0 ) return 1.0;

	double_array log_a( boost::extents[ 1 ][ 1 ] );
	calculate_log_a( m, log_a );

	double_array log_b( boost::extents[ 1 ][ 1 ] );
	calculate_log_b( m, log_b );

	double_array phi( boost::extents[ T ][ N ] );
	unsigned_array psi( boost::extents[ T ][ N ] );

	//initialisation
	for( unsigned i = 0; N != i; ++i )
	{
		phi[ 0 ][ i ] = myrrh::safe_log( m.pi(i) ) + log_b[i][ O[0] ];
		psi[ 0 ][ i ] = 0;
	}

	//more efficient if we just examine those states that can preceed any given state so calculate them here
	unsigned_vec_vec predecessors;
	calculate_predecessors( m, predecessors );

	//recursion
	//for each time point
	for( unsigned t = 1; T != t; ++t )
	{
		const unsigned output = O[t];
		double_array::reference phi_t = phi[t];
		double_array::reference phi_t_1 = phi[t-1];
		unsigned_array::reference psi_t = psi[t];
		unsigned_array::reference psi_t_1 = psi[t-1];

		//for each state
		for( unsigned j = 0; N != j; ++j )
		{
			double max = - std::numeric_limits< double >::max();
			unsigned max_state = N;
			//for each predecessor
			BOOST_FOREACH( unsigned i, predecessors[j] )
			//for( unsigned i = 0; N != i; ++i )
			{
				//if the predecessor is not a possible state then don't bother
				if( N == psi_t_1[i] ) continue;

				//otherwise check if this is the best state path to this state
				const double tmp = phi_t_1[i] + log_a[i][j];
				if( tmp >= max )
				{
					max = tmp;
					max_state = i;
				}
			}
			phi_t[j] = max + log_b[j][output];
			psi_t[j] = max_state;
		}
	}

	//termination
	double_array::const_reference::const_iterator max_element = std::max_element( phi[T-1].begin(), phi[T-1].end() );
	const double P_star = *max_element;
	q_star[T-1] = max_element - phi[T-1].begin();

	//path (state sequence) backtracking
	for( unsigned t = T - 1; 0 != t; --t )
	{
		const unsigned best_state = psi[t][ q_star[t] ];
		if( N == best_state )
			throw std::logic_error( "No Viterbi state path possible for this model and this sequence" );
		q_star[t-1] = best_state;
	}

	return P_star;
}


} //namespace hmm


#ifdef _MSC_VER
#pragma float_control(pop)
#endif //_MSC_VER



#endif //HMM_VITERBI_H_

