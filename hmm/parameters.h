/**
@file

Copyright John Reid 2007

*/

#ifndef HMM_PARAMETERS_H_
#define HMM_PARAMETERS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "hmm/defs.h"

#include <myrrh/math.h>




namespace hmm {

typedef boost::optional< unsigned > optional_param_idx;
typedef std::vector< optional_param_idx > param_idx_vec;
typedef boost::multi_array< optional_param_idx, 2 > param_idx_array;


/**
Takes a vector of optional indexes into a parameter vector.

Calculates the sum of those parameters.
*/
template< 
	typename Parameters,
	typename Parameterisation
>
double sum_parameters( 
	const Parameterisation & parameterisation,
	Parameters & parameters )
{
	using namespace boost;

	//what is the sum of the indexed parameters?
	double sum = 0.0;
	BOOST_FOREACH( const typename range_value< Parameterisation >::type & p, parameterisation ) {
	    if( p ) {
	        sum += parameters[ *p ];
	    }
	}

	HMM_VERIFY( MYRRH_ISFINITE( sum ) );

	return sum;
}


/**
Takes a vector of optional indexes into a parameter vector.

Scales the indexed parameters such the values indexed in the original vector sum to 1.0

If they sum to 0.0 currently then make a uniform distribution.
*/
template< 
	typename Parameters,
	typename Parameterisation
>
void normalise_parameters( 
	const Parameterisation & parameterisation,
	Parameters & parameters )
{
	using namespace boost;

	//what is the sum of the indexed parameters?
	const double sum = sum_parameters( parameterisation, parameters );

	//the parameters sum to 0
	if( 0.0 == sum ) 
	{
		//how many references to parameters are there?
		unsigned references = 0;
		BOOST_FOREACH( const typename range_value< Parameterisation >::type & p, parameterisation ) if( p ) ++references;

		//is there something to normalise?
		if( 0 != references )
		{
			const double uniform = 1.0 / double( references );
			BOOST_FOREACH( const typename range_value< Parameterisation >::type & p, parameterisation ) if( p ) parameters[ *p ] = uniform;
		}
	}
	else
	{
		BOOST_FOREACH( const typename range_value< Parameterisation >::type & p, parameterisation ) if( p ) parameters[ *p ] /= sum;
	}
}



} //namespace hmm


#endif //HMM_PARAMETERS_H_

