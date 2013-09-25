/**
@file

Copyright John Reid 2006

*/

#ifndef HMM_PSSM_MODEL_H_
#define HMM_PSSM_MODEL_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/defs.h"
#include "hmm/model.h"

#include "myrrh/math.h"

namespace hmm {


namespace impl {

struct state_map
{
	unsigned L;
	state_map( unsigned L ) : L( L ) { }
	inline unsigned state( unsigned l, bool positive ) const { return 1 + l + ( positive ? 0 : L ); }
	inline unsigned l( unsigned state ) const { return ( state - 1 ) % L; }
	inline bool positive( unsigned state ) const { return state > 1 + L; }
};

inline
void
fill_in_pssm_model_structure(
	model & m )
{
	if( m.N % 2 != 1 ) throw std::logic_error( "Expecting model to have an odd number of states" );

	const unsigned L = ( m.N - 1 ) / 2; //# bases in pssm
	const unsigned M = 4;
	const impl::state_map _state_map( L );


	//for each pssm base and the background distribution we have 4 emission parameters
	//each pssm base has one transition parameter, the background state has 3
	//and one parameter for each initial start position
	HMM_VERIFY( m.P == (L + 1) * 4 + L + 3 + m.N );

	unsigned p = 0;

	//the pi (initial distribution) parameters
	for( unsigned i = 0; m.N != i; ++i ) 
	{
		m.pi_parameterisation[ i ] = p;
		++p;
	}

	//the background state emission parameters
	for( unsigned k = 0; M != k; ++k )
	{
		m.b_parameterisation[ 0 ][ k ] = p;
		++p;
	}

	//the pssm state emission parameters
	for( unsigned l = 0; L != l; ++l )
	{
		const unsigned positive_state = _state_map.state( l, true );
		const unsigned negative_state = _state_map.state( l, false );

		for( unsigned k = 0; M != k; ++k )
		{
			m.b_parameterisation[ positive_state ][ k ] = p;
			m.b_parameterisation[ negative_state ][ M - k - 1 ] = p;
			++p;
		}
	}

	//the background state transition parameters
	const unsigned first_positive_state = _state_map.state( 0, true );
	const unsigned first_negative_state = _state_map.state( L - 1, false );
	m.a_parameterisation[ 0 ][ 0 ] = p;
	++p;
	m.a_parameterisation[ 0 ][ first_positive_state ] = p;
	++p;
	m.a_parameterisation[ 0 ][ first_negative_state ] = p;
	++p;

	//the pssm state transition parameters
	for( unsigned l = 0; L != l; ++l )
	{
		const unsigned positive_state = _state_map.state( l, true );
		const unsigned negative_state = _state_map.state( l, false );

		const unsigned next_positive_state = L - 1 == l ? 0 : _state_map.state( l + 1, true );
		const unsigned next_negative_state = 0 == l ? 0 : _state_map.state( l - 1, false );

		m.a_parameterisation[ positive_state ][ next_positive_state ] = p;
		m.a_parameterisation[ negative_state ][ next_negative_state ] = p;
		++p;
	}

}

} //namespace impl



/**
Builds a hmm that models a sequence with some binding sites in it according to the given pssm position distributions.

The HMM models binding sites on the positive and negative strands.
*/
model::ptr
build_pssm_model(
	double p_binding_site,
	const double_array & position_distributions )
{
	using namespace impl;

	//the number of bases in the pssm
	const unsigned L = position_distributions.size();
	const state_map _state_map( L );

	//we have one background state and 2 states for each base in the pssm (+ve and -ve strand)
	const unsigned N = 1 + L * 2;
	//state 0 is the background state
	//states 1..L are the +ve strand states
	//states L+1..L*2 are the -ve strand states

	//we can output any of 4 nucleotides
	const unsigned M = 4;

	//for each pssm base and the background distribution we have 4 emission parameters
	//each pssm base has one transition parameter, the background state has 3
	//we have one pi parameter for each state
	const unsigned P = (L + 1) * 4 + L + 3 + N;

	model::ptr result( new model( N, M, P ) );

	//fill in the structure - tying parameters together
	fill_in_pssm_model_structure( *result );

	//set the initial distribution to something sensible
	result->theta[ *( result->pi_parameterisation[ 0 ] ) ] = 1.0 - p_binding_site;
	for( unsigned i = 1; N != i; ++ i ) result->theta[ *( result->pi_parameterisation[ i ] ) ] = p_binding_site / double( N - 1 );

	//the background state emission parameters
	for( unsigned k = 0; M != k; ++k ) result->theta[ *( result->b_parameterisation[ 0 ][ k ] ) ] = 0.25;

	//the pssm state emission parameters
	for( unsigned l = 0; L != l; ++l )
	{
		const unsigned positive_state = _state_map.state( l, true );
		for( unsigned k = 0; M != k; ++k ) result->theta[ *( result->b_parameterisation[ positive_state ][ k ] ) ] = position_distributions[ l ][ k ];
	}

	//the background state transition parameters
	const unsigned first_positive_state = _state_map.state( 0, true );
	const unsigned first_negative_state = _state_map.state( L - 1, false );
	result->theta[ *( result->a_parameterisation[ 0 ][ 0 ] ) ] = 1.0 - p_binding_site;
	result->theta[ *( result->a_parameterisation[ 0 ][ first_positive_state ] ) ] = p_binding_site / 2;
	result->theta[ *( result->a_parameterisation[ 0 ][ first_negative_state ] ) ] = p_binding_site / 2;

	//the pssm state transition parameters
	for( unsigned l = 0; L != l; ++l )
	{
		const unsigned positive_state = _state_map.state( l, true );
		const unsigned next_positive_state = L - 1 == l ? 0 : _state_map.state( l + 1, true );
		result->theta[ *( result->a_parameterisation[ positive_state ][ next_positive_state ] ) ] = 1.0;
	}

	check_consistent( *result );

	return result;
}


/**
Builds a hmm that models a sequence with some binding sites in it according to a random pssm specified by a dirichlet pssm_prior.

The HMM models binding sites on the positive and negative strands.
*/
model::ptr
build_random_pssm_model(
	double p_binding_site,
	unsigned K,
	double pssm_prior = 1e-2 )
{
	double_array position_distributions( boost::extents[ K ][ 4 ] );
	const double_vec dirichlet_params( 4, pssm_prior );
	for( unsigned k = 0; K != k; ++k ) 
	{
		double_array::reference position_distributions_k = position_distributions[ k ];
		myrrh::gsl::draw_from_dirichlet( dirichlet_params, position_distributions_k );
	}

	return build_pssm_model( p_binding_site, position_distributions );
}

/**
Copies and shifts a HMM that models a pssm +/- shift_offset bases.
*/
model::ptr
shift_pssm_model(
	const model & original,
	int shift_offset,
	const double_vec & fill_in_distribution )
{
	const unsigned M = 4;
	HMM_VERIFY( M == fill_in_distribution.size() );

	model::ptr shifted( new model( original.N, original.M, original.P ) );

	impl::fill_in_pssm_model_structure( *shifted );

	//the background state emission parameters
	for( unsigned k = 0; M != k; ++k ) shifted->theta[ *( shifted->b_parameterisation[ 0 ][ k ] ) ] = original.b( 0, k );

	//pi for background state
	shifted->pi_parameterisation[ 0 ] = original.pi_parameterisation[ 0 ];

	//the background state transition parameters
	const unsigned L = ( original.N - 1 ) / 2; //# bases in pssm
	const impl::state_map _state_map( L );
	const unsigned first_positive_state = _state_map.state( 0, true );
	const unsigned first_negative_state = _state_map.state( L - 1, false );
	shifted->theta[ *( shifted->a_parameterisation[ 0 ][ 0 ] ) ] = original.a( 0, 0 );
	shifted->theta[ *( shifted->a_parameterisation[ 0 ][ first_positive_state ] ) ] = original.a( 0, first_positive_state );
	shifted->theta[ *( shifted->a_parameterisation[ 0 ][ first_negative_state ] ) ] = original.a( 0, first_negative_state );

	//for each base in the new pssm
	for( unsigned l = 0; L != l; ++l )
	{
		const unsigned positive_state = _state_map.state( l, true );
		const unsigned negative_state = _state_map.state( l, false );

		//which base in the original pssm does l correspond to
		const int original_l = l + shift_offset;
		const unsigned original_positive_state = _state_map.state( (L + original_l) % L, true );
		const unsigned original_negative_state = _state_map.state( (L + original_l) % L, false );
		HMM_VERIFY( 0 < original_positive_state );
		HMM_VERIFY( 0 < original_negative_state );

		//does this base, l, in the shifted pssm correspond to a base in the original one?
		const bool has_original_state = 0 <= original_l && original_l < int( L );

		//emission parameters
		if( has_original_state ) for( unsigned k = 0; M != k; ++k ) shifted->theta[ *( shifted->b_parameterisation[ positive_state ][ k ] ) ] = original.b( original_positive_state, k );
		else for( unsigned k = 0; M != k; ++k ) shifted->theta[ *( shifted->b_parameterisation[ positive_state ][ k ] ) ] = fill_in_distribution[ k ];

		//pi - initial distribution
		shifted->pi_parameterisation[ positive_state ] = original.pi_parameterisation[ original_positive_state ];
		shifted->pi_parameterisation[ negative_state ] = original.pi_parameterisation[ original_negative_state ];
		if( ! has_original_state )
		{
			//reset old initial distribution - (1 - pi(background))/N
			shifted->theta[ *( shifted->pi_parameterisation[ positive_state ] ) ] 
				= shifted->theta[ *( shifted->pi_parameterisation[ negative_state ] ) ] 
				= ( 1.0 - shifted->theta[ *( shifted->pi_parameterisation[ 0 ] ) ] ) / ( 2 * original.N );
		}

		//the state transition parameter
		const unsigned next_positive_state = L - 1 == l ? 0 : _state_map.state( l + 1, true );
		shifted->theta[ *( shifted->a_parameterisation[ positive_state ][ next_positive_state ] ) ] = 1.0;
	}

	//scale the initial distribution into a proper distribution
	normalise_parameters( shifted->pi_parameterisation, shifted->theta );

	check_consistent( *shifted );

	return shifted;
}


} //namespace hmm

#endif //HMM_PSSM_MODEL_H_

