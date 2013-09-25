/**
@file

Copyright John Reid 2007

*/

#ifndef HMM_MODEL_CONVERSION_H_
#define HMM_MODEL_CONVERSION_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/model.h"
#include "hmm/model_state.h"

#include <myrrh/math.h>

namespace hmm {

template< typename T >
struct ptr_value_less_than
{
	bool operator()( T p1, T p2 ) const { return p1.get() < p2.get(); }
};

/**
Maps states to indices.
*/
typedef std::map< model_state::ptr, unsigned, ptr_value_less_than< model_state::ptr > > model_state_2_index_map;

/**
Calculates indices for all states in the model and populates the map.
*/
inline 
void calculate_indices_for_states( 
	const model_defined_by_states & model,
	model_state_2_index_map & index_map )
{
	index_map.clear();
	for( unsigned i = 0; model.N() != i; ++i ) index_map[ model.states[ i ] ] = i;
}

/**
Gets the index of a state or throws an error.
*/
inline 
unsigned get_index_of_state( const model_state_2_index_map & index_map, model_state::ptr state )
{
	model_state_2_index_map::const_iterator i = index_map.find( state );
	if( index_map.end() == i ) throw std::logic_error( "Could not find state in map" );
	return i->second;
}


/**
Creates a standard model from a model defined by its states.
*/
inline
model::ptr
model_states_2_model(
	const model_defined_by_states & states_model )
{
	const unsigned N = states_model.N();
	const unsigned M = states_model.M();
	const unsigned P = states_model.parameters.size();

	//allocate result
	model::ptr _result( new model( N, M, P, states_model.converter.order ) );
	model & result( *_result );

	//calculate an index for each state
	model_state_2_index_map state_2_index;
	calculate_indices_for_states( states_model, state_2_index );

	//copy the parameters
	for( unsigned p = 0; P != p; ++p ) result.theta[ p ] = states_model.parameters[ p ];

	//define the transitions and emissions for each state
	for( unsigned i = 0; N != i; ++i )
	{
		model_state::ptr state = states_model.states[ i ];
		HMM_VERIFY( state_2_index[ state ] == i );

		//for initial distribution
		result.pi_parameterisation[ i ] = state->pi;

		//for each transition
		BOOST_FOREACH( model_state::successor & successor, state->successors )
		{
			if( states_model.states.end() == std::find( states_model.states.begin(), states_model.states.end(), successor.state ) ) throw std::logic_error( "Cannot find state at all" );
			const unsigned j = get_index_of_state( state_2_index, successor.state );

			if( successor.a ) result.a_parameterisation[ i ][ j ] = successor.a;
		}

		//for each emission
		for( unsigned k = 0; M != k; ++k )
		{
			if( state->b[ k ] ) result.b_parameterisation[ i ][ k ] = state->b[ k ];
		}
	}

	normalise_model_parameters( result );

	return _result;
}

/**
Creates a model defined by its states from a standard model.
*/
inline
model_defined_by_states::ptr
model_2_model_states(
	const model & model )
{
	//allocate the result
	model_defined_by_states::ptr _result( new model_defined_by_states( model.M, model.converter.order ) );
	model_defined_by_states & result( *_result );

	//copy the parameters
	std::copy( model.theta.begin(), model.theta.end(), std::back_inserter( result.parameters ) );
	HMM_VERIFY( model.theta.size() == result.parameters.size() );

	//add each state
	for( unsigned i = 0; model.N != i; ++i ) 
		result.states.push_back( 
			model_state::ptr(
				new model_state( 
					result.parameters, 
					model.M,
					model.pi_parameterisation[ i ] ) ) );

	//add the transitions
	for( unsigned i = 0; model.N != i; ++i ) 
	{
		for( unsigned j = 0; model.N != j; ++j )
		{
			if( model.a_parameterisation[ i ][ j ] )
			{
				result.states[ i ]->successors.push_back(
					model_state::successor(
						result.states[ j ],
						model.a_parameterisation[ i ][ j ] ) ) ;
			}
		}
	}

	//add the emissions
	for( unsigned i = 0; model.N != i; ++i ) 
	{
		for( unsigned k = 0; model.M != k; ++k )
		{
			if( model.b_parameterisation[ i ][ k ] ) result.states[ i ]->b[ k ] = model.b_parameterisation[ i ][ k ];
		}
	}

	return _result;
}



} //namespace hmm


#endif //HMM_MODEL_CONVERSION_H_

