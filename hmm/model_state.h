/**
@file

Copyright John Reid 2007

Code to define a HMM by its states and their connections. I.e. a graph oriented data structure c.f. the model oriented data
structure in model.h

*/

#ifndef HMM_MODEL_STATE_H_
#define HMM_MODEL_STATE_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/defs.h"
#include "hmm/parameters.h"
#include "hmm/markov_order.h"

#include <myrrh/defs.h>
#include <myrrh/math.h>

namespace hmm {

using myrrh::double_vec;
using myrrh::double_array;

using myrrh::unsigned_vec;
using myrrh::unsigned_vec_vec;
using myrrh::unsigned_array;


struct model_state
{
	typedef boost::shared_ptr< model_state > ptr;
	typedef std::vector< ptr > ptr_vec;

	struct successor
	{
		typedef std::vector< successor > vec;

		inline 
		successor( ptr state, optional_param_idx a ) 
			: state( state )
			, a( a )
		{
			if( ! state ) throw std::logic_error( "Null successor state" );
			if( ! a ) throw std::logic_error( "No parameter for this transition" );
		}

		ptr state;									/**< Succeeding state. */
		optional_param_idx a;						/**< Transition probability. */

		bool operator==( const successor & x ) const { return state == x.state; }
		bool operator!=( const successor & x ) const { return state != x.state; }
		bool operator<( const successor & x ) const { return state < x.state; }
	};

	double_vec & parameters;						/**< Parameters. */
	optional_param_idx pi;							/**< Chance of being in this state initially. */
	successor::vec successors;						/**< States we can translate to. */
	param_idx_vec b;								/**< Emission probabilities. */

	inline model_state( 
		double_vec & parameters, 
		unsigned M /* output_alphabet_size */, 
		optional_param_idx pi = optional_param_idx() ) 
		: parameters( parameters )
		, pi( pi )
		, b( M )
	{ }

	inline unsigned M() const { return b.size(); }

	inline void add_successor(
		ptr state, 
		optional_param_idx a ) 
	{
		successors.push_back( successor( state, a ) );
	}

	inline double a( const successor & s ) { return s.a ? parameters[ *( s.a ) ] : 0.0; }

	inline bool is_consistent() const
	{
		return true;
	}
};

struct model_defined_by_states
	: boost::noncopyable
{
	typedef boost::shared_ptr< model_defined_by_states > ptr;
	typedef markov_order_n_converter< unsigned > markov_converter;

	markov_converter converter;
	model_state::ptr_vec states;
	double_vec parameters;

	model_defined_by_states( unsigned M, unsigned markov_order = 0 ) 
		: converter( markov_converter::order_0_size_from_order_n_size( markov_order, M ), markov_order ) 
	{ }

	inline unsigned N() const { return states.size(); }
	inline unsigned M() const { return converter.order_n_size; }

	inline model_state::ptr add_state( optional_param_idx pi = optional_param_idx() ) 
	{
		states.push_back( model_state::ptr( new model_state( parameters, M(), pi ) ) );
		return states.back();
	}

	inline optional_param_idx add_parameter( double value = 0.0 ) 
	{
		parameters.push_back( value );
		return optional_param_idx( parameters.size() - 1 );
	}
};


/**
Finds all the predecessors of a given state in the model and populates the vector.
*/
inline 
void fill_with_predecessors(
	const model_defined_by_states & model,
	model_state::ptr state, 
	model_state::ptr_vec & predecessors ) 
{
	predecessors.clear();
	for( unsigned i = 0; model.N() != i; ++i )
	{
		model_state::ptr potential_predecessor = model.states[ i ];
		BOOST_FOREACH( const model_state::successor & successor, potential_predecessor->successors )
		{
			if( state == successor.state )
			{
				predecessors.push_back( potential_predecessor );
			}
		}
	}
}

inline
void 
check_consistent(
	const model_defined_by_states & model )
{
	BOOST_FOREACH( model_state::ptr state, model.states )
	{
		if( ! state ) throw std::logic_error( "Null state" );

		BOOST_FOREACH( const model_state::successor & successor, state->successors )
			if( model.states.end() == std::find( model.states.begin(), model.states.end(), successor.state ) )
				throw std::logic_error( "Could not find successor" );
	}
}


} //namespace hmm


#endif //HMM_MODEL_STATE_H_

