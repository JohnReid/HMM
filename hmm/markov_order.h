/**
@file

Copyright John Reid 2007

Code to convert between sequences of different markov orders.

I.e. the 0-order sequence A,C,G,T could be converted to the 1-order sequence N,AC,CG,GT and back again to N,C,G,T

*/

#ifndef HMM_MARKOV_ORDER_H_
#define HMM_MARKOV_ORDER_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER


#include "hmm/defs.h"

#include <myrrh/types.h>
#include <myrrh/math.h>

#include <deque>

namespace hmm {
namespace impl {


/** Calculates if the type is large enough to hold the data of markov order order. */
template< typename char_t >
bool
can_type_hold_markov_value_of_order( char_t alphabet_size, unsigned order )
{
	const double log_largest_char = ( order + 1 ) * std::log( double( alphabet_size ) );
	const double log_largest_possible = std::log( double( std::numeric_limits< char_t >::max() ) );
	return log_largest_char <= log_largest_possible;
}


} //namespace impl




template<
	typename CharT = unsigned
>
struct markov_order_n_converter
{
	typedef CharT char_t;

	unsigned order;								/**< The order of the markov observations. */
	char_t order_0_size;						/**< The size of the order 0 alphabet. */
	char_t order_n_size;						/**< The size of the order order alphabet. */
	char_t order_n_minus_1_size;				/**< The size of the order order-1 alphabet. */

	/** Serialise this model. */
	template< typename Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & order;
        ar & order_0_size;
        ar & order_n_size;
        ar & order_n_minus_1_size;
	}

	static
	char_t
	order_0_size_from_order_n_size( unsigned order, char_t order_n_size )
	{
		const double log_order_n_size = std::log( double( order_n_size ) );
		const double log_order_0_size = log_order_n_size / (order + 1);
		const double order_0_size = std::exp( log_order_0_size );
		const char_t result = char_t( myrrh::nearest_integer( order_0_size ) );
		if( order_n_size != order_n_size_from_order_0_size( order, result ) )
			throw std::logic_error( "order_n_size incorrectly specified" );
		return result;
	}

	static
	char_t
	order_n_size_from_order_0_size( unsigned order, char_t order_0_size )
	{
		return char_t( myrrh::int_power( order_0_size, order + 1 ) );
	}

	/** Useless default constructor to help with serialisation. */
	markov_order_n_converter()
		: order( 0 )
		, order_0_size( 0 )
		, order_n_size( 0 )
		, order_n_minus_1_size( 0 )
	{
	}

	/** Constructor. */
	markov_order_n_converter( char_t alphabet_size, unsigned order )
		: order( order )
		, order_0_size( alphabet_size )
		, order_n_size( order_n_size_from_order_0_size( order, alphabet_size ) )
		, order_n_minus_1_size( 0 == order ? 1 : order_n_size_from_order_0_size( order - 1, alphabet_size ) )
	{
		if(
			! impl::can_type_hold_markov_value_of_order< char_t >(
				alphabet_size,
				order ) )
		{
			throw std::logic_error( HMM_MAKE_STRING( "Cannot store characters of markov order " << order << " in this type" ) );
		}
	}

	bool operator==( const markov_order_n_converter< char_t > & rhs ) const {
		return
			order == rhs.order
			&& order_0_size == rhs.order_0_size
			&& order_n_size == rhs.order_n_size
			&& order_n_minus_1_size == rhs.order_n_minus_1_size
			;
	}

	/**
	Makes order + 1 markov order 0 observations into a one observation of order order.
	*/
	template< typename It >
	char_t
	convert_to_order_n( It begin, It end ) const
	{
		if( end - begin != int( order + 1 ) ) throw std::logic_error( HMM_MAKE_STRING( "Input sequence must have " << order + 1 << " characters" ) );
		char_t power = char_t( 1 );
		char_t result = 0;
		while( begin != end )
		{
			//is the input unknown?
			if( typename It::value_type( order_0_size ) == *begin ) return order_n_size;
			result += power * *begin;
			power *= order_0_size;
			++begin;
		}
		return result;
	}


	/**
	Gets the last order 0 observation from an order n input
	*/
	char_t
	get_last_order_0( char_t input ) const
	{
		return order_n_size == input ? order_0_size : (input % order_0_size);
	}

	/**
	Gets the previous order n-1 observation from an order n input
	*/
	char_t
	get_previous_order_n_minus_1( char_t input ) const
	{
		return order_n_size == input ? order_n_minus_1_size : (input / order_0_size);
	}

	/**
	Gets the last order n-1 observation from an order n input
	*/
	char_t
	get_last_order_n_minus_1( char_t input ) const
	{
		return order_n_size == input ? order_n_minus_1_size : (input % order_n_minus_1_size);
	}

	/**
	Makes one markov order n observation into its n + 1 constituent observations of order 0. (outputs the observations in reverse order)
	*/
	template< typename OutputIt >
	void
	convert_to_order_0_reversed(
		char_t input,
		OutputIt output_it ) const
	{
		//is the input unknown?
		if( order_n_size == input )
		{
			//yes
			for( unsigned i = 0; order + 1 != i; ++i )
			{
				*output_it = order_0_size; //the output is unknown
				++output_it;
			}
		}
		else
		{
			//no
			for( unsigned i = 0; order + 1 != i; ++i )
			{
				*output_it = input % order_0_size;
				++output_it;
				input /= order_0_size;
			}
		}
	}

	/**
	Converts a 0-order sequence to a sequence of order n.
	*/
	template<
		typename InputIt,
		typename OutputIt
	>
	void
	make_markov_order_n(
		InputIt begin,
		InputIt end,
		OutputIt output_it ) const
	{
		//keeps history of previous characters
		std::deque< char_t > history( order, order_0_size );

		while( begin != end )
		{
			history.push_front( *begin );
			*output_it = convert_to_order_n( history.begin(), history.end() );
			++output_it;
			++begin;
			history.pop_back();
		}
	}

	/**
	Converts a order-n sequence to a sequence of order 0.
	*/
	template<
		typename InputIt,
		typename OutputIt
	>
	void
	make_markov_order_0(
		InputIt begin,
		InputIt end,
		OutputIt output_it ) const
	{
		while( begin != end )
		{
			*output_it = get_last_order_0( *begin );
			++output_it;
			++begin;
		}
	}

	/**
	Is the input sequence a valid order n sequence?
	*/
	template< typename InputIt >
	bool
	is_valid_order_n_sequence(
		InputIt begin,
		InputIt end ) const
	{
		char_t last_order_n = order_n_size;
		while( begin != end )
		{
			if( char_t( *begin ) > order_n_size ) return false;

			//if both this and last obs are known then check they match
			if( order_n_size != char_t( *begin )
				&& order_n_size != last_order_n
				&& get_last_order_n_minus_1(last_order_n) != get_previous_order_n_minus_1(*begin) )
			{
				return false;
			}
			last_order_n = *begin;
			++begin;
		}
		return true;
	}

	/**
	Is the input sequence a valid order 0 sequence?
	*/
	template< typename InputIt >
	bool
	is_valid_order_0_sequence(
		InputIt begin,
		InputIt end ) const
	{
		for( ;begin != end; ++begin) if( *begin > order_0_size ) return false;
		return true;
	}

	/**
	Number of known bases of order n sequence.
	*/
	template<
		typename InputIt
	>
	unsigned
	num_known_bases_order_n(
		InputIt begin,
		InputIt end ) const
	{
		unsigned count = 0;
		for( ; begin != end; ++begin ) if( order_n_size != *begin ) ++count;
		return count;
	}

	/**
	Number of known bases of order 0 sequence.
	*/
	template<
		typename InputIt
	>
	unsigned
	num_known_bases_order_0(
		InputIt begin,
		InputIt end ) const
	{
		unsigned count = 0;
		for( ; begin != end; ++begin ) if( order_0_size != *begin ) ++count;
		return count;
	}
};

} //namespace hmm





#endif //HMM_MARKOV_ORDER_H_

