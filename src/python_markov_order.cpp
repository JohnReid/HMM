/**
@file

Copyright John Reid 2007, 2012

*/

#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"
#include "hmm/python.h"

#include <myrrh/python/numpy.h>

namespace hmm {


HMM_EXPORT_FN_SPEC
void 
export_markov_order()
{
	using namespace boost::python;
	using namespace myrrh::python;

	using boost::python::arg;
	typedef impl::python_markov_order_n_converter< unsigned > python_converter;
	class_<
		python_converter
	>( 
		"MarkovOrderConverter",
		"Converts between sequences of order n and order 0",
		init< unsigned, unsigned >( ( arg( "alphabet_size" ), "order" ) )
	)
	.def_readonly(
		"n",
		&python_converter::order,
		"markov order for the converter" )
	.def_readonly(
		"order_0_size",
		&python_converter::order_0_size,
		"# characters in order 0 alphabet (also the value to represent missing data in an order 0 sequence)" )
	.def_readonly(
		"order_n_size",
		&python_converter::order_n_size,
		"# characters in order n alphabet (also the value to represent missing data in an order n sequence)" )
	.def_readonly(
		"order_n_minus_1_size",
		&python_converter::order_n_minus_1_size,
		"# characters in order n-1 alphabet" )
	.def(
		"num_known_bases_order_0",
		&python_converter::num_known_bases_order_0_python,
		"# known bases in a sequence of order 0",
		( arg( "sequence" ) ) )
	.def(
		"num_known_bases_order_n",
		&python_converter::num_known_bases_order_n_python,
		"# known bases in a sequence of order n",
		( arg( "sequence" ) ) )
	.def(
		"convert_from_order_n_observation",
		&python_converter::convert_from_order_n_observation,
		"Converts one order n observation to n+1 order 0 observations",
		( arg( "order_n_obs" ) ) )
	.def(
		"convert_to_order_n_observation",
		&python_converter::convert_to_order_n_observation,
		"Converts n+1 order 0 observations to one order n observation",
		( arg( "order_0_seq" ) ) )
	.def(
		"to_order_n",
		&python_converter::to_order_n,
		"Converts a sequence of markov order 0 to markov order n",
		( arg( "sequence" ) ) )
	.def(
		"to_order_0",
		&python_converter::to_order_0,
		"Converts a sequence of markov order n to markov order 0",
		( arg( "sequence" ) ) )
	.def(
		"make_order_0",
		&python_converter::get_last_order_0,
		"Returns the last order 0 observation that makes up the order n observation",
		( arg( "order_n_obs" ) ) )
	.def(
		"is_valid_order_0_sequence",
		&python_converter::is_valid_order_0_sequence,
		"Is the input sequence a valid order 0 sequence?",
		( arg( "order_0_obs" ) ) )
	.def(
		"is_valid_order_n_sequence",
		&python_converter::is_valid_order_n_sequence,
		"Is the input sequence a valid order n sequence?",
		( arg( "order_n_obs" ) ) )
	;
}

} //namespace hmm
