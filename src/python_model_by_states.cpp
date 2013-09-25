/**
@file

Copyright John Reid 2007, 2012

*/



#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"

#include "hmm/python.h"
#include "hmm/model_state.h"
#include "hmm/model_conversion.h"

namespace boost { namespace python { namespace indexing {
template<>
struct value_traits< hmm::model_state::ptr > : public simple_value_traits< hmm::model_state::ptr >
{
	static bool const equality_comparable = false;
	static bool const lessthan_comparable = false;
};
} } }

namespace hmm {

namespace impl {

impl::python_markov_order_n_converter< unsigned >
wrapped_converter( const model_defined_by_states & m )
{
	return impl::python_markov_order_n_converter< unsigned >( m.converter.order_0_size, m.converter.order );
}

} //namespace impl


template< typename PtrT >
bool python_pointer_eq( boost::python::object obj_1, boost::python::object obj_2 )
{
	boost::python::extract< PtrT > x1( obj_1 );
	boost::python::extract< PtrT > x2( obj_2 );
	return x1.check() && x2.check() && static_cast< PtrT >( x1 ).get() == static_cast< PtrT >( x2 ).get();
}

template< typename PtrT >
bool python_pointer_cmp( boost::python::object obj_1, boost::python::object obj_2 )
{
	boost::python::extract< PtrT > x1( obj_1 );
	boost::python::extract< PtrT > x2( obj_2 );
	return x1.check() && x2.check() && static_cast< PtrT >( x1 ).get() < static_cast< PtrT >( x2 ).get();
}

template< typename PtrT >
int python_pointer_hash( boost::python::object obj_1 )
{
	boost::python::extract< PtrT > x1( obj_1 );
	return boost::hash_value( static_cast< PtrT >( x1 ).get() );
}

HMM_EXPORT_FN_SPEC
void 
export_model_by_states()
{
	using namespace boost::python;
	using namespace indexing;
	using namespace myrrh::python;
	
	using boost::python::arg;
	class_< model_state::successor >( 
		"ModelStateSuccessor",
		"Defines a transition",
		init< model_state::ptr, optional_param_idx >( ( arg("state"), arg("parameterisation") ) )
	)
		.add_property( 
			"state",
			&model_state::successor::state,
			"The succeeding state" )
		.add_property( 
			"a",
			&model_state::successor::a,
			"The transition parameter index" )
		;
	
	class_< model_state::successor::vec >( "ModelStateSuccessorVec" )
		.def( container_suite< model_state::successor::vec >() )
	;

	class_< model_state >( 
		"ModelState",
		"A state in a hidden markov model",
		init< double_vec &, unsigned, optional_param_idx >( ( arg( "parameters" ), "M", "pi" ) )
	)
		.def( 
			"__eq__",
			python_pointer_eq< model_state::ptr >,
			"Equality through C++ identity" )
		.def( 
			"__hash__",
			python_pointer_hash< model_state::ptr >,
			"Hash" )
		.def_readwrite( 
			"b",
			&model_state::b,
			"Initial state probability" )
		.def_readwrite( 
			"pi",
			&model_state::pi,
			"Initial state probability" )
		.def_readwrite( 
			"successors",
			&model_state::successors,
			"The potential successors of this state" )
		.add_property( 
			"M",
			&model_state::M,
			"The size of the output alphabet for the state" )
		.def( 
			"add_successor",
			&model_state::add_successor,
			( arg( "state" ), "parameter" ),
			"Add a successor to this state" )
		.def( 
			"a",
			&model_state::a,
			( arg( "successor" ) ),
			"p(transition to given successor)" )
		;
	
	class_< model_state::ptr_vec >( "ModelStateVec" )
		.def( container_suite< model_state::ptr_vec >() )
	;

	register_ptr_to_python< model_state::ptr >();

	typedef void ( * check_consistent_fn )( const model_defined_by_states & );
	class_< model_defined_by_states, boost::noncopyable >( 
		"ModelByStates",
		"A hidden markov model defined by its states",
		init< unsigned, unsigned >( 
			( arg( "M" ), arg( "markov_order" ) = 0 ),
			"M is the output alphabet size and markov_order is the markov order of the model." 
		)
	)
		.def_readwrite( 
			"states",
			&model_defined_by_states::states,
			"The states in the model" )
		.def_readwrite( 
			"parameters",
			&model_defined_by_states::parameters,
			"The parameters of the model" )
		.add_property( 
			"N",
			&model_defined_by_states::N,
			"The number of states in the model" )
		.add_property( 
			"M",
			&model_defined_by_states::M,
			"The size of the output alphabet for the model"  )
		.add_property(
			"converter",
			impl::wrapped_converter,
			"Converts between order 0 emissions and order n emissions" )
		.def( 
			"add_state",
			&model_defined_by_states::add_state,
			( arg( "pi" ) = optional_param_idx() ),
			"Add one state to the model" )
		.def( 
			"add_parameter",
			&model_defined_by_states::add_parameter,
			( arg( "value" ) = 0.0 ),
			"Add one parameter to the model" )
		.def( 
			"check_consistent",
			check_consistent_fn( check_consistent ),
			"Checks the model is consistent" )
		;
	register_ptr_to_python< model_defined_by_states::ptr >();

	def(
		"model_states_2_model",
		model_states_2_model,
		"Converts a model defined by its states to a standard model" );

	def(
		"model_2_model_states",
		model_2_model_states,
		"Converts a standard model to a model defined by its states" );
}

} //namespace hmm
