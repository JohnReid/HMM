/**
@file

Copyright John Reid 2007

*/

#ifndef HMM_MODEL_H_
#define HMM_MODEL_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/defs.h"
#include "hmm/markov_order.h"
#include "hmm/parameters.h"

#include <myrrh/defs.h>
#include <myrrh/math.h>
#include <myrrh/serialize/multi_array.h>
#include <myrrh/serialize/optional.h>

namespace hmm {




using myrrh::double_vec;
using myrrh::double_array;

using myrrh::unsigned_vec;
using myrrh::unsigned_vec_vec;
using myrrh::unsigned_array;


template< typename T, std::size_t D >
void
fill_multi_array(
	typename boost::multi_array< T, D > & a,
	const T & value )
{
	std::fill(
		a.origin(),
		a.origin() + a.num_elements(),
		value );
}



struct model
{
	typedef boost::shared_ptr< model > ptr;
	typedef markov_order_n_converter< unsigned > markov_converter;

	/** Prior for model. */
	struct prior
	{
		double_array a;							/**< Prior for transition probabilities. */
		double_array b;							/**< Prior for emission probabilities. */
		double_vec pi;							/**< Prior for initial state distribution. */

		prior( unsigned N, unsigned M )
			: a( boost::extents[ N ][ N ] )
			, b( boost::extents[ N ][ M ] )
			, pi( N, 0.0 )
		{
			fill_multi_array( a, 0.0 );
			fill_multi_array( b, 0.0 );
		}

		unsigned N() const { return a.shape()[0]; }
		unsigned M() const { return b.shape()[1]; }

		bool matches_model( const model & m ) const { return N() == m.N && M() == m.M; }
		void check_matches_model( const model & m ) const { if( ! matches_model(m) ) throw std::logic_error( "Prior does not match model" ); }
	};

	markov_converter converter;					/**< Converts between emissions. */
	unsigned N;									/**< # of states in the model. */
	unsigned M;									/**< # of characters in the output alphabet of the model. */
	unsigned P;									/**< # of parameters. */
	double_vec theta;							/**< Parameters for emission and transition probabilities. */
	param_idx_vec pi_parameterisation;			/**< Parameters for initial state distribution. */
	double pi_normaliser;						/**< 1.0/sum(pi) to normalise pi distribution. */
	param_idx_array a_parameterisation;			/**< How transitions are coupled to parameters. */
	double_vec a_normaliser;					/**< 1.0/sum(a[i]) to normalise transition distribution. */
	param_idx_array b_parameterisation;			/**< How emissions are coupled to parameters. */
	double_array b_normaliser;					/**< 1.0/sum(b[i]) to normalise emission distribution. */

	/** Useless default constructor to help with serialisation. */
	model()
		: N( 0 )
		, M( 0 )
		, P( 0 )
		, theta( 0, 0.0 )
		, pi_parameterisation( 0 )
		, pi_normaliser( 1. )
		, a_parameterisation( boost::extents[ 0 ][ 0 ] )
		, a_normaliser( 0, 1.0 )
		, b_parameterisation( boost::extents[ 0 ][ 0 ] )
		, b_normaliser( boost::extents[ 0 ][ 0 ] )
	{
	}

	/** Constructs a model of the given sizes. */
	model( unsigned N, unsigned M, unsigned P, unsigned markov_order = 0 )
		: converter( markov_converter::order_0_size_from_order_n_size( markov_order, M ), markov_order )
		, N( N )
		, M( M )
		, P( P )
		, theta( P, 0.0 )
		, pi_parameterisation( N )
		, pi_normaliser( 1.0 )
		, a_parameterisation( boost::extents[ N ][ N ] )
		, a_normaliser( N, 1.0 )
		, b_parameterisation( boost::extents[ N ][ M ] )
		, b_normaliser( boost::extents[ N ][ num_emission_distributions() ] )
	{
		fill_multi_array( b_normaliser, 1.0 );
	}

	/** Serialise this model. */
	template< typename Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & converter;
        ar & N;
        ar & M;
        ar & P;
        ar & theta;
        ar & pi_parameterisation;
        ar & pi_normaliser;
        ar & a_parameterisation;
        ar & a_normaliser;
        ar & b_parameterisation;
        ar & b_normaliser;
    }

	//copy constructor and operator= should work fine as compiler implemented

	/** The number of emission distributions (due to the emission markov order). */
	inline unsigned num_emission_distributions() const
	{
		return converter.order_n_minus_1_size;
	}

	/** Get the parameter for the optional idx. */
	inline double get_parameter( optional_param_idx idx ) const
	{
		return idx ? theta[ *idx ] : 0.0;
	}

	/** Throws logic_error if i is not a valid state index. */
	inline void check_state( unsigned i ) const
	{
		if( i >= N ) throw std::logic_error( HMM_MAKE_STRING( "State " << i << " outside bounds [0," << N << ")" ) );
	}

	/** Throws logic_error if k is not a valid output index. */
	inline void check_output( unsigned k ) const
	{
		if( k >= M ) throw std::logic_error( HMM_MAKE_STRING( "Output character " << k << " outside bounds [0," << M << ")" ) );
	}

	/** p( transition from state i to state j ) */
	inline double a( unsigned i, unsigned j ) const
	{
		const double result = get_parameter( a_parameterisation[ i ][ j ] ) * a_normaliser[ i ];
		HMM_VERIFY( myrrh::is_probability( result ) );
		return result;
	}

	/** p( output k from state i ) */
	inline double b( unsigned i, unsigned k ) const
	{
		if( M == k ) return 1.0; //k == M is unknown observation, so p = 1
		HMM_VERIFY( k < M );
		const double result = get_parameter( b_parameterisation[ i ][ k ] ) * b_normaliser[ i ][ converter.get_previous_order_n_minus_1(k) ];
		HMM_VERIFY( myrrh::is_probability( result ) );
		return result;
	}

	/** p( in state i at time 0 ) */
	inline double pi( unsigned i ) const
	{
		const double result = get_parameter( pi_parameterisation[ i ] ) * pi_normaliser;
		HMM_VERIFY( myrrh::is_probability( result ) );
		return result;
	}

	/** Same as a but arguments are checked for validity. */
	inline double a_checked( unsigned i, unsigned j ) const
	{
		check_state( i );
		check_state( j );
		return a( i, j );
	}

	/** Same as b but arguments are checked for validity. */
	inline double b_checked( unsigned i, unsigned k ) const
	{
		check_state( i );
		check_output( k );
		return b( i, k );
	}

	void
	set_transition(
		unsigned i,
		unsigned j,
		double p )
	{
		check_state( i );
		check_state( j );
		if( ! a_parameterisation[i][j] ) throw std::logic_error( HMM_MAKE_STRING("No parameter for transition from state "<<i<<" to state "<<j) );
		theta[ *(a_parameterisation[i][j]) ] = p;
	}

	void
	set_emission(
		unsigned i,
		unsigned k,
		double p )
	{
		check_state( i );
		check_output( k );
		if( ! b_parameterisation[i][k] ) throw std::logic_error( HMM_MAKE_STRING("No parameter for emission "<<k<<" from state "<<i) );
		theta[ *(b_parameterisation[i][k]) ] = p;
	}

	void
	set_initial(
		unsigned i,
		double p )
	{
		check_state( i );
		if( ! pi_parameterisation[i] ) throw std::logic_error( HMM_MAKE_STRING("No parameter for state "<<i<<" in initial distribution") );
		theta[ *(pi_parameterisation[i]) ] = p;
	}

};


inline void check_no_nans( const double_array & A )
{
#ifdef _DEBUG
	for( unsigned i = 0; A.shape()[0] != i; ++i )
		for( unsigned j = 0; A.shape()[1] != j; ++j )
			HMM_VERIFY( MYRRH_ISFINITE( A[i][j] ) );
#endif
}

inline void check_no_nans( const double_vec & v )
{
#ifdef _DEBUG
	BOOST_FOREACH( double x, v ) HMM_VERIFY( MYRRH_ISFINITE( x ) );
#endif
}

inline
void
calculate_model_transition_matrix(
	const model & m,
	double_array & A )
{
	const unsigned N = m.N;
	A.resize( boost::extents[ N ][ N ] );
	for( unsigned i = 0; N != i; ++i )
		for( unsigned j = 0; N != j; ++j )
			A[ i ][ j ] = m.a( i, j );
	check_no_nans( A );
}


inline
void
calculate_model_emission_matrix(
	const model & m,
	double_array & B )
{
	const unsigned N = m.N;
	const unsigned M = m.M;
	B.resize( boost::extents[ N ][ M ] );
	for( unsigned i = 0; N != i; ++i ) {
		for( unsigned k = 0; M != k; ++k ) {
			B[ i ][ k ] = m.b( i, k );
		}
	}
	check_no_nans( B );
}

/**
Calculates the model emission matrix including an entry for unknown output set to 1.0
*/
inline
void
calculate_model_emission_matrix_with_unknown(
	const model & m,
	double_array & B )
{
	const unsigned N = m.N;
	const unsigned M = m.M;
	B.resize( boost::extents[ N ][ M+1 ] );
	for( unsigned i = 0; N != i; ++i ) {
		for( unsigned k = 0; M != k; ++k ) {
			B[ i ][ k ] = m.b( i, k );
		}
		B[ i ][ M ] = 1.0;
	}
	check_no_nans( B );
}

inline
void
calculate_model_initial_distribution(
	const model & m,
	double_vec & pi )
{
	const unsigned N = m.N;
	pi.resize( N );
	for( unsigned i = 0; N != i; ++i ) pi[ i ] = m.pi( i );
	check_no_nans( pi );
}

inline
void
check_consistent(
	const model & m )
{
	double_array A( boost::extents[ m.N ][ m.N ] );
	calculate_model_transition_matrix( m, A );
	unsigned i = 0;
	BOOST_FOREACH( double_array::reference row, A )
	{
		const double sum = std::accumulate( row.begin(), row.end(), 0.0 );
		if( ! myrrh::is_close( 1.0, sum, 0.01 ) ) throw std::logic_error( MYRRH_MAKE_STRING( "HMM not consistent: A row "<<i<<" sum = "<<sum ) );
		unsigned j = 0;
		BOOST_FOREACH( double a, row )
		{
			if( ! myrrh::is_probability( a ) ) throw std::logic_error( MYRRH_MAKE_STRING( "HMM not consistent: A("<<i<<","<<j<<") = "<<a ) );
			++j;
		}
		++i;
	}

	double_array B( boost::extents[ m.N ][ m.M ] );
	calculate_model_emission_matrix( m, B );
	i = 0;
	BOOST_FOREACH( double_array::reference row, B )
	{
		typedef double_array::value_type::const_iterator it;
		it begin = row.begin();
		it end = begin + m.converter.order_0_size;
		while( begin != row.end() )
		{
			const double sum = std::accumulate( begin, end, 0.0 );
			if( ! myrrh::is_close( 1.0, sum, 0.01 ) ) throw std::logic_error( MYRRH_MAKE_STRING( "HMM not consistent: B row "<<i<<" sum = "<<sum ) );
			begin = end;
			end += m.converter.order_0_size;
		}
		unsigned j = 0;
		BOOST_FOREACH( double b, row )
		{
			if( ! myrrh::is_probability( b ) ) throw std::logic_error( MYRRH_MAKE_STRING( "HMM not consistent: B("<<i<<","<<j<<") = "<<b ) );
			++j;
		}
		++i;
	}

	double sum = 0.0;
	for( i = 0; m.N != i; ++i ) sum += m.pi( i );
	if( ! myrrh::is_close( 1.0, sum, 0.01 ) ) throw std::logic_error( MYRRH_MAKE_STRING( "HMM not consistent: pi sum = " << sum ) );
}


namespace detail {

inline double normaliser_from_sum( double sum )
{
	if( 0.0 == sum )
		throw std::logic_error( "normaliser_from_sum: sum == 0.0" );
	const double one_over_sum = 1.0 / sum;
	if( ! MYRRH_ISFINITE( one_over_sum ) )
		throw std::logic_error( HMM_MAKE_STRING( "normaliser_from_sum: one_over_sum is not finite: sum="<<sum ) );
	return one_over_sum;
}

}

inline
void
normalise_model_parameters( model & m )
{
	//normalise parameters for a
	for( unsigned i = 0; m.N != i; ++i ) {
		try {
			m.a_normaliser[i] = detail::normaliser_from_sum( sum_parameters( m.a_parameterisation[i], m.theta ) );
		} catch( const std::exception & e ) {
			throw std::logic_error( MYRRH_MAKE_STRING( "Could not normalise transition parameters for state "<<i<<": "<<e.what() ) );
		}
	}
	//for( unsigned i = 0; m.N != i; ++i ) normalise_parameters( m.a_parameterisation[ i ], m.theta );

	//normalise parameters for b
	for( unsigned i = 0; m.N != i; ++i )
	{
		try {
			typedef param_idx_array::value_type::const_iterator it;
			it begin = m.b_parameterisation[i].begin();
			it end = begin + m.converter.order_0_size;
			for( unsigned q = 0; begin != m.b_parameterisation[i].end(); ++q )
			{
				m.b_normaliser[i][q] = detail::normaliser_from_sum(
					sum_parameters(
						boost::make_iterator_range( begin, end ),
						m.theta ) );
				begin = end;
				end += m.converter.order_0_size;
			}
		} catch( const std::exception & e ) {
			throw std::logic_error( MYRRH_MAKE_STRING( "Could not normalise emission parameters for state "<<i<<": "<<e.what() ) );
		}
	}
	//for( unsigned i = 0; m.N != i; ++i ) normalise_parameters( m.b_parameterisation[ i ], m.theta );

	//normalise pi
	try {
		if( m.N ) m.pi_normaliser = detail::normaliser_from_sum( sum_parameters( m.pi_parameterisation, m.theta ) );
	} catch( const std::exception & e ) {
		throw std::logic_error( MYRRH_MAKE_STRING( "Could not normalise initial state parameters: "<<e.what() ) );
	}
	//normalise_parameters( m.pi_parameterisation, m.theta );
}



inline
model::ptr
create_random_fully_connected_model(
	const unsigned N,
	const unsigned M,
	unsigned markov_order = 0 )
{
	using namespace boost;
	using namespace myrrh::gsl;

	const unsigned P = N * N + N * M + N;

	model::ptr result( new model( N, M, P, markov_order ) );

	unsigned p = 0; //index into the parameters

	//sample pi
	const double_vec pi_dirichlet_params( N, 1.0 );
	double_vec pi( N );
	draw_from_dirichlet( pi_dirichlet_params, pi );
	for( unsigned i = 0; N != i; ++i )
	{
		result->pi_parameterisation[ i ] = p;
		result->theta[ p ] = pi[ i ];
		++p;
	}

	//sample a and b
	const double_vec a_dirichlet_params( N, 1.0 );
	const double_vec b_dirichlet_params( M, 1.0 );
	for( unsigned i = 0; N != i; ++i )
	{
		boost::iterator_range< double_vec::iterator > draw_range_a( result->theta.begin() + p, result->theta.begin() + p + N );
		draw_from_dirichlet( a_dirichlet_params, draw_range_a );
		for( unsigned j = 0; N != j; ++j ) result->a_parameterisation[ i ][ j ] = p++;

		boost::iterator_range< double_vec::iterator > draw_range_b( result->theta.begin() + p, result->theta.begin() + p + M );
		draw_from_dirichlet( b_dirichlet_params, draw_range_b );
		for( unsigned k = 0; M != k; ++k ) result->b_parameterisation[ i ][ k ] = p++;
	}

	HMM_VERIFY( P == p );

	normalise_model_parameters( *result );
	check_consistent( *result );

	return result;
}

/**
Examines a model and calculates which states can follow each state. Stores result in unsigned_vec_vec.

I.e. if j = successors[ i ][ x ] for some x, then j is a successor of state i.
*/
inline
void
calculate_successors(
	const model & m,
	unsigned_vec_vec & successors )
{
	const unsigned N = m.N;
	successors.resize( N );
	for( unsigned i = 0; N != i; ++i )
	{
		successors[ i ].clear();
		for( unsigned j = 0; N != j; ++j ) if( m.a_parameterisation[ i ][ j ] ) successors[ i ].push_back( j );
	}
}


/**
Examines a model and calculates which states can preceed each state. Stores result in unsigned_vec_vec.

I.e. if j = predecessors[ i ][ x ] for some x, then state j preceeds state i.
*/
inline
void
calculate_predecessors(
	const model & m,
	unsigned_vec_vec & predecessors )
{
	const unsigned N = m.N;
	predecessors.resize( N );
	for( unsigned j = 0; N != j; ++j )
	{
		predecessors[ j ].clear();
		for( unsigned i = 0; N != i; ++i ) if( m.a_parameterisation[ i ][ j ] ) predecessors[ j ].push_back( i );
	}
}


} //namespace hmm


#endif //HMM_MODEL_H_

