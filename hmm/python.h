/**
@file

Copyright John Reid 2006, 2007, 2008

*/

#ifndef HMM_PYTHON_H_
#define HMM_PYTHON_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "boost/python/detail/wrap_python.hpp"

#include <myrrh/python/convert.h>
#include <myrrh/python/tuple.h>

#include "hmm/defs.h"
#include "hmm/markov_order.h"

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/refcount.hpp>
#include <boost/python/suite/indexing/vector.hpp>


/*
#ifdef HMM_BUILD_PYTHON_EXPORT
# define HMM_EXPORT_FN_SPEC __declspec( dllexport )
#else
# define HMM_EXPORT_FN_SPEC __declspec( dllimport )
#endif
*/
#define HMM_EXPORT_FN_SPEC


namespace hmm {
namespace impl {

typedef std::vector< output_seq > sequence_vec;
typedef boost::shared_ptr< output_seq > output_seq_shared_ptr;
typedef boost::shared_ptr< sequence_vec > sequence_vec_shared_ptr;
typedef std::vector< double > double_vec;
typedef boost::shared_ptr< double_vec > double_vec_ptr;

sequence_vec_shared_ptr
extract_or_convert_sequences( boost::python::object seqs_python );


output_seq_shared_ptr
extract_or_convert_sequence( boost::python::object seq_python );

} //namespace impl


template< typename T >
std::string 
as_string( const T & t )
{
	return HMM_MAKE_STRING( t );
}

HMM_EXPORT_FN_SPEC void export_model();
HMM_EXPORT_FN_SPEC void export_model_by_states();
HMM_EXPORT_FN_SPEC void export_markov_order();
HMM_EXPORT_FN_SPEC void export_simple_pssm();


void std_exception_translator( std::logic_error const & x );
void c_string_translator( const char * x );
void string_translator( const std::string & x );


inline 
boost::python::tuple
boost_to_python_tuple( const boost::tuples::null_type & ) 
{
	return boost::python::make_tuple();
}

template< typename H, typename T >
inline 
boost::python::object 
boost_to_python_tuple( const boost::tuples::cons< H, T > & x ) 
{
	return boost::python::make_tuple( x.get_head() ) + boost_to_python_tuple( x.get_tail() );
}


template< typename T >
struct tupleconverter
{
	static PyObject * convert( T const & x )
	{
		return boost::python::incref( boost_to_python_tuple( x ).ptr() );
	}
};


namespace impl {
using namespace boost::python;
using namespace myrrh::python;

template<
	typename CharT = unsigned
>
struct python_markov_order_n_converter : markov_order_n_converter< CharT >
{
	typedef std::vector< CharT > vec;
	typedef boost::shared_ptr< vec > vec_ptr;
	typedef std::deque< CharT > deque;
	typedef markov_order_n_converter< CharT > base_t;
	
	python_markov_order_n_converter( unsigned alphabet_size, unsigned n ) : base_t( alphabet_size, n ) { }

	object convert_from_order_n_observation( unsigned order_n_obs ) const
	{
		deque output_seq_C;
		base_t::convert_to_order_0_reversed( order_n_obs, std::front_inserter( output_seq_C ) );
		return myrrh::python::convert_to_python( output_seq_C );
	}

	unsigned convert_to_order_n_observation( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );
		if( ! base_t::is_valid_order_0_sequence( input_seq_C.begin(), input_seq_C.end() ) ) throw std::logic_error( "Is not a order 0 sequence" );
		if( base_t::order+1 != input_seq_C.size() ) throw std::logic_error( 
			HMM_MAKE_STRING( "Input sequence is wrong length, expecting exactly " << base_t::order+1 << " observations" ) );

		vec output_seq_C;
		make_markov_order_n(
			input_seq_C.begin(),
			input_seq_C.end(),
			back_inserter( output_seq_C ) );

		return output_seq_C.back();
	}

	vec_ptr to_order_n( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );
		if( ! base_t::is_valid_order_0_sequence( input_seq_C.begin(), input_seq_C.end() ) ) throw std::logic_error( "Is not a order 0 sequence" );

		vec_ptr output_seq_C( new vec );
		make_markov_order_n(
			input_seq_C.begin(),
			input_seq_C.end(),
			back_inserter( *output_seq_C ) );

		return output_seq_C;
	}

	vec_ptr to_order_0( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );
		if( ! base_t::is_valid_order_n_sequence( input_seq_C.begin(), input_seq_C.end() ) ) throw std::logic_error( "Is not a order n sequence" );

		vec_ptr output_seq_C( new vec );
		make_markov_order_0(
			input_seq_C.begin(),
			input_seq_C.end(),
			back_inserter( *output_seq_C ) );

		return output_seq_C;
	}

	unsigned num_known_bases_order_n_python( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );

		return num_known_bases_order_n( input_seq_C.begin(), input_seq_C.end() );
	}

	unsigned is_valid_order_n_sequence( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );

		return base_t::is_valid_order_n_sequence( input_seq_C.begin(), input_seq_C.end() );
	}

	unsigned is_valid_order_0_sequence( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );

		return base_t::is_valid_order_0_sequence( input_seq_C.begin(), input_seq_C.end() );
	}

	unsigned num_known_bases_order_0_python( object input_seq_python ) const
	{
		vec input_seq_C;
		convert_from_python( input_seq_python, input_seq_C );

		return num_known_bases_order_0( input_seq_C.begin(), input_seq_C.end() );
	}
};



} //namespace impl



} //namespace hmm

#endif //HMM_PYTHON_H_


