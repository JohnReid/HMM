/**
@file

Copyright John Reid 2007, 2012

*/

#ifdef _MSC_VER
#pragma warning( disable : 4172 )
#endif //_MSVC




#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"

#include "hmm/python.h"
#include "hmm/model.h"
#include "hmm/forward_backward.h"
#include "hmm/viterbi.h"
#include "hmm/viterbi_training.h"
#include "hmm/baum_welch.h"
#include "hmm/pssm_model.h"
#include "hmm/sample.h"
#include "hmm/count_mers.h"

#include <myrrh/math.h>
#include <myrrh/serialisable.h>
#include <myrrh/python/numpy.h>
#include <myrrh/python/multi_array_to_numpy.h>
#include <myrrh/python/boost_range.h>

namespace myrrh {
namespace python {

numpy_converter converter;
std::string exposed_typechars;

} //namespace myrrh
} //namespace myrrh


namespace hmm {

using myrrh::python::make_boost_range;

namespace impl {


template< typename T >
boost::python::object as_writeable_numpy( T const& t )
{
	boost::python::object result = myrrh::python::convert_to_python( t );
	result.attr( "setflags" )( true );
	return result;
}


sequence_vec_shared_ptr convert_sequences( boost::python::object seqs_python )
{
	sequence_vec_shared_ptr result;
	const unsigned num_seqs = boost::python::len( seqs_python );
	result.reset( new sequence_vec( num_seqs ) );
	for( unsigned i = 0; num_seqs != i; ++i )
		myrrh::python::convert_from_python( seqs_python[ i ], (*result)[ i ] );
	return result;
}

output_seq_shared_ptr convert_sequence( boost::python::object seq_python )
{
	output_seq_shared_ptr result( new output_seq );
	myrrh::python::convert_from_python( seq_python, *result );
	return result;
}


sequence_vec_shared_ptr
extract_or_convert_sequences( boost::python::object seqs_python )
{
	boost::python::extract< sequence_vec_shared_ptr > e( seqs_python );
	return e.check() ? e() : convert_sequences( seqs_python );
}


output_seq_shared_ptr
extract_or_convert_sequence( boost::python::object seq_python )
{
	boost::python::extract< output_seq_shared_ptr > e( seq_python );
	//std::cout << "Sequence extraction "<<(e.check()?"succeeds":"fails")<<"\n";
	return e.check() ? e() : convert_sequence( seq_python );
}


} //namespace impl





/** Pickles vectors. */
template< typename Container >
struct vector_pickle_suite : boost::python::pickle_suite
{
	typedef typename Container::value_type value_type;

    static
    boost::python::tuple
    getstate(const Container & c)
    {
		return boost::python::tuple( c );
    }

    static
    void
    setstate(Container & c, boost::python::tuple state)
    {
		const unsigned l = boost::python::len(state);
		c.resize(l);
		for( unsigned i = 0; l != i; ++i )
		{
			value_type el = boost::python::extract< value_type >( state[i] );
			c[i] = el;
		}
    }
};





/** Pickles models. */
struct model_pickle_suite : boost::python::pickle_suite
{
    static
    boost::python::tuple
    getinitargs(const model & m)
    {
		return boost::python::make_tuple( m.N, m.M, m.P, m.converter.order );
    }

    static
    boost::python::tuple
    getstate(const model & m)
    {
		return boost::python::make_tuple(
			myrrh::python::convert_to_python( m.theta ),
			m.pi_parameterisation,
			m.pi_normaliser,
			m.a_parameterisation,
			myrrh::python::convert_to_python( m.a_normaliser ),
			m.b_parameterisation,
			myrrh::python::convert_to_python( m.b_normaliser )
		);
    }

    static
    void
    setstate(model & m, boost::python::tuple state)
    {
		if( 7 != boost::python::len( state ) ) throw std::logic_error( "Pickling error: Expecting 7 tuple items for model" );

		using boost::python::extract;
		myrrh::python::convert_from_python( state[0], m.theta );
		m.pi_parameterisation = extract< param_idx_vec >( state[1] );
		m.pi_normaliser = extract< double >( state[2] );
		m.a_parameterisation = ( const param_idx_array & )( extract< param_idx_array >( state[3] ) );
		myrrh::python::convert_from_python( state[4], m.a_normaliser );
		m.b_parameterisation = ( const param_idx_array & )( extract< param_idx_array >( state[5] ) );
		myrrh::python::convert_from_python( state[6], m.b_normaliser );
    }
};



/** Pickles optionals. */
template< typename T >
struct optional_pickle_suite : boost::python::pickle_suite
{
    static
    boost::python::tuple
	getstate(const boost::optional<T> & i)
    {
		return i ? boost::python::make_tuple( true, *i ) : boost::python::make_tuple( false );
    }

    static
    void
    setstate(boost::optional<T> & i, boost::python::tuple state)
    {
		if( state[0] ) i = T( boost::python::extract< T >( state[1] ) );
		else i = boost::optional<T>();
    }
};


optional_param_idx
get_transition_parameterisation(
	const model & m,
	unsigned i,
	unsigned j )
{
	m.check_state( i );
	m.check_state( j );
	return m.a_parameterisation[ i ][ j ];
}


void
set_transition_parameterisation(
	model & m,
	unsigned i,
	unsigned j,
	optional_param_idx p )
{
	m.check_state( i );
	m.check_state( j );
	m.a_parameterisation[ i ][ j ] = p;
}



optional_param_idx
get_emission_parameterisation(
	const model & m,
	unsigned i,
	unsigned k )
{
	m.check_state( i );
	m.check_output( k );
	return m.b_parameterisation[ i ][ k ];
}

void
set_emission_parameterisation(
	model & m,
	unsigned i,
	unsigned k,
	optional_param_idx p )
{
	m.check_state( i );
	m.check_output( k );
	m.b_parameterisation[ i ][ k ] = p;
}



model::ptr
wrapped_build_pssm_model(
	double p_binding_site,
	const boost::python::object & position_distributions )
{
	double_array tmp;
	myrrh::python::convert_from_python( position_distributions, tmp );

	return build_pssm_model(
		p_binding_site,
		tmp );
}

model::ptr
wrapped_shift_pssm_model(
	const model & original_model,
	int shift_offset,
	const boost::python::object & fill_in_distribution )
{
	double_vec tmp;
	myrrh::python::convert_from_python( fill_in_distribution, tmp );

	return shift_pssm_model(
		original_model,
		shift_offset,
		tmp );
}

void
convert_seq_from_python( const model & m, boost::python::object O_python, output_seq & O_C )
{
	myrrh::python::convert_from_python( O_python, O_C );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C.begin(), O_C.end() ) );
}



boost::python::tuple
wrapped_viterbi(
	const model & m,
	boost::python::object O_python )
{
	check_consistent( m );

	impl::output_seq_shared_ptr O_C = impl::extract_or_convert_sequence( O_python );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C->begin(), O_C->end() ) );

	unsigned_vec q_star;
	const double log_P_star = viterbi( m, *O_C, q_star );

	return boost::python::make_tuple( log_P_star, myrrh::python::convert_to_python( q_star ) );
}






// a seq of integral types without unsigned counterparts
#define HMM_PYTHON_basic_ints            (char)(short)(int)(long)(long long)

// generate a seq containing "signed t" and "unsigned t"
#define HMM_PYTHON_int_pair(r,data,t)      (signed t)(unsigned t)

// a seq of all the integral types
#define HMM_PYTHON_ints                                           \
    BOOST_PP_SEQ_FOR_EACH(HMM_PYTHON_int_pair, ~, HMM_PYTHON_basic_ints)

/** Does not copy sequences - uses new myrrh conversion to avoid copying. */
boost::python::tuple
wrapped_viterbi_2(
	const model & m,
	boost::python::object O_py )
{
	check_consistent( m );

	unsigned_vec q_star;
	double log_P_star;
	bool work_done = false;
	// generate the function to do the work for type T
	#define HMM_PYTHON_do_work(r,data,T) { \
		boost::python::extract< boost::multi_array_ref< T, 1 > > O_extractor( O_py ); \
		if( ! work_done && O_extractor.check() ) { \
			const boost::multi_array_ref< T, 1 > & O = O_extractor(); \
			HMM_VERIFY( m.converter.is_valid_order_n_sequence( O.begin(), O.end() ) ); \
			log_P_star = viterbi( m, O, q_star ); \
			work_done = true; \
	    } \
	}
	BOOST_PP_SEQ_FOR_EACH(HMM_PYTHON_do_work, ~, HMM_PYTHON_ints)
	#undef HMM_PYTHON_do_work
	if( ! work_done ) {
		throw std::logic_error( "Cannot convert python object to sequence for HMM library." );
	}

	return boost::python::make_tuple( log_P_star, myrrh::python::convert_to_python( q_star ) );
}

#undef HMM_PYTHON_ints
#undef HMM_PYTHON_int_pair
#undef HMM_PYTHON_basic_ints









boost::python::tuple
wrapped_forward(
	const model & m,
	boost::python::object O_python,
	double scaling_threshold )
{
	impl::output_seq_shared_ptr O_C = impl::extract_or_convert_sequence( O_python );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C->begin(), O_C->end() ) );

	double_array alpha;
	double_vec c;
	const double log_p_O = forward( m, *O_C, alpha, c, scaling_threshold );

	return boost::python::make_tuple( log_p_O, myrrh::python::convert_to_python( alpha ), myrrh::python::convert_to_python( c ) );
}

double
wrapped_LL(
	const model & m,
	boost::python::object O_python,
	double scaling_threshold )
{
	impl::output_seq_shared_ptr O_C = impl::extract_or_convert_sequence( O_python );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C->begin(), O_C->end() ) );

	double_array alpha;
	double_vec c;
	const double log_p_O = forward( m, *O_C, alpha, c, scaling_threshold );

	return log_p_O;
}

double
wrapped_LL_multi(
	const model & m,
	boost::python::object py_seqs,
	double scaling_threshold )
{
	double result = 0.0;
	double_array alpha;
	double_vec c;
	BOOST_FOREACH( impl::output_seq_shared_ptr c_O, make_boost_range< impl::output_seq_shared_ptr >( py_seqs ) ) {
	//BOOST_FOREACH( const output_seq & O_C, *seqs_C )
		HMM_VERIFY( m.converter.is_valid_order_n_sequence( c_O->begin(), c_O->end() ) );
		result += forward( m, *c_O, alpha, c, scaling_threshold );
	}

	return result;
}

boost::python::object
wrapped_backward(
	const model & m,
	boost::python::object O_python,
	boost::python::object c_python )
{
	impl::output_seq_shared_ptr O_C = impl::extract_or_convert_sequence( O_python );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C->begin(), O_C->end() ) );

	double_vec c_C;
	myrrh::python::convert_from_python< double_vec >( c_python, c_C );

	double_array beta;
	backward( m, *O_C, c_C, beta );

	return myrrh::python::convert_to_python( beta );
}

boost::python::object
wrapped_forward_backward(
	const model & m,
	boost::python::object O_python,
	double scaling_threshold  )
{
	impl::output_seq_shared_ptr O_C = impl::extract_or_convert_sequence( O_python );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C->begin(), O_C->end() ) );

	double_array alpha;
	double_vec c;
	const double log_p_O = forward( m, *O_C, alpha, c, scaling_threshold );
	double_array beta;
	backward( m, *O_C, c, beta );

	return boost::python::make_tuple(
		log_p_O,
		myrrh::python::convert_to_python( alpha ),
		myrrh::python::convert_to_python( beta ),
		myrrh::python::convert_to_python( c ) );
}

void
wrapped_baum_welch_updater(
	model & m,
	boost::python::object O_python,
	boost::python::object alpha_python,
	boost::python::object beta_python,
	boost::python::object c_python )
{
	impl::output_seq_shared_ptr O_C = impl::extract_or_convert_sequence( O_python );
	HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C->begin(), O_C->end() ) );

	double_array alpha_C;
	myrrh::python::convert_from_python( alpha_python, alpha_C );

	double_array beta_C;
	myrrh::python::convert_from_python( beta_python, beta_C );

	double_vec c_C;
	myrrh::python::convert_from_python( c_python, c_C );

	baum_welch_updater( m, *O_C, alpha_C, beta_C, c_C );
}

double
wrapped_baum_welch_iteration(
	model & m,
	boost::python::object py_seqs,
	boost::python::object prior_python )
{
#ifdef HMM_VERIFY_ON
	BOOST_FOREACH( const output_seq & O_C, make_boost_range< const output_seq & >( py_seqs ) )
		HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C.begin(), O_C.end() ) );
#endif //HMM_VERIFY_ON

	model::prior prior_C = prior_python ? boost::python::extract< model::prior >(prior_python) : model::prior(m.N, m.M);

	return baum_welch_multi_sequence_iteration( m, prior_C, make_boost_range< const output_seq & >( py_seqs ) );
}

template< typename T >
struct python_function_result_extractor_1
{
	boost::python::object _fn;
	python_function_result_extractor_1( boost::python::object fn ) : _fn( fn ) { }
	template< typename Arg1 >
	T operator()( Arg1 arg_1 ) const
	{
		return boost::python::extract< T >( _fn( arg_1 ) );
	}
};



training_result_tuple
wrapped_baum_welch(
	model & m,
	boost::python::object py_seqs,
	boost::python::object prior_python,
	unsigned max_iterations,
	double tolerance,
	boost::python::object callback )
{
#ifdef HMM_VERIFY_ON
	BOOST_FOREACH( const output_seq & O_C, make_boost_range< const output_seq & >( py_seqs ) )
		HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C.begin(), O_C.end() ) );
#endif //HMM_VERIFY_ON

	model::prior prior_C = prior_python ? boost::python::extract< model::prior >(prior_python) : model::prior(m.N, m.M);

	return
		callback
			? baum_welch_multi_sequence(
				m,
				prior_C,
				make_boost_range< const output_seq & >( py_seqs ),
				max_iterations,
				tolerance,
				python_function_result_extractor_1< bool >( callback ) )
			: baum_welch_multi_sequence(
				m,
				prior_C,
				make_boost_range< const output_seq & >( py_seqs ),
				max_iterations,
				tolerance,
				impl::null_callback() );
}

training_result_tuple
wrapped_viterbi_training(
	model & m,
	boost::python::object py_seqs,
	boost::python::object prior_python,
	unsigned max_iterations,
	double tolerance,
	boost::python::object callback )
{
#ifdef HMM_VERIFY_ON
	BOOST_FOREACH( const output_seq & O_C, make_boost_range< const output_seq & >( py_seqs ) )
		HMM_VERIFY( m.converter.is_valid_order_n_sequence( O_C.begin(), O_C.end() ) );
#endif //HMM_VERIFY_ON

	model::prior prior_C = prior_python ? boost::python::extract< model::prior >(prior_python) : model::prior(m.N, m.M);

	return
		callback
			? viterbi_training_multi_sequence(
				m,
				prior_C,
				make_boost_range< const output_seq & >( py_seqs ),
				max_iterations,
				tolerance,
				python_function_result_extractor_1< bool >( callback ) )
			: viterbi_training_multi_sequence(
				m,
				prior_C,
				make_boost_range< const output_seq & >( py_seqs ),
				max_iterations,
				tolerance,
				impl::null_callback() );
}

boost::python::object
transition_matrix( const model & m )
{
	double_array A;
	calculate_model_transition_matrix( m, A );
	return myrrh::python::convert_to_python( A );
}

boost::python::object
emission_matrix( const model & m )
{
	double_array B;
	calculate_model_emission_matrix( m, B );
	return myrrh::python::convert_to_python( B );
}

boost::python::object
initial_distribution( const model & m )
{
	double_vec pi;
	calculate_model_initial_distribution( m, pi );
	return myrrh::python::convert_to_python( pi );
}

boost::python::object
wrapped_sample( const model & m, unsigned T )
{
	unsigned_vec state_seq;
	output_seq output;
	sample( m, T, state_seq, output );
	return boost::python::make_tuple( myrrh::python::convert_to_python( state_seq ), myrrh::python::convert_to_python( output ) );
}

typedef boost::shared_ptr< std::vector< double > > double_vec_ptr;
double_vec_ptr
wrapped_draw_from_dirichlet( boost::python::object params_python )
{
	double_vec params_C;
	myrrh::python::convert_from_python( params_python, params_C );
	double_vec_ptr result( new double_vec( params_C.size() ) );
	myrrh::gsl::draw_from_dirichlet( params_C, *result );
	return result;
}

impl::python_markov_order_n_converter< unsigned >
wrapped_converter( const model & m )
{
	return impl::python_markov_order_n_converter< unsigned >( m.converter.order_0_size, m.converter.order );
}

struct count_mers_output
{
	boost::python::object callback;
	count_mers_output( boost::python::object callback ) : callback( callback ) { }

	template< typename Prefix >
	void operator()( const Prefix & prefix, mer_count count ) {
		callback( myrrh::python::convert_to_python( prefix ), count.total );
	}

	template< typename Prefix >
	void operator()( const Prefix & prefix, const mer_count_with_seq_membership & count ) {
		callback( myrrh::python::convert_to_python( prefix ), count.total, count.seq_idx_membership.total_set() );
	}
};

void
wrapped_count_mers( boost::python::object py_seqs, unsigned n, boost::python::object callback )
{
	count_mers_output output( callback );
	count_mers< mer_count >( make_boost_range< const output_seq & >( py_seqs ), n, output );
}


void
wrapped_count_mers_with_sequence( boost::python::object py_seqs, unsigned n, boost::python::object callback )
{
	count_mers_output output( callback );
	count_mers< mer_count_with_seq_membership >( make_boost_range< const output_seq & >( py_seqs ), n, output );
}

struct top_mers_callback {
	typedef markov_order_n_converter<> converter_t;

	boost::python::list _k_mers;
	converter_t _converter;
	output_seq _k_mer;

	top_mers_callback( unsigned k )
		: _converter( 4, k-1 )
		, _k_mer( k )
	{ }

	void operator()( converter_t::char_t prefix, const mer_count_with_seq_membership & count ) {
		_converter.convert_to_order_0_reversed( prefix, _k_mer.begin() );
		_k_mers.append( boost::python::make_tuple( myrrh::python::convert_to_python( _k_mer ), count.total, count.seq_idx_membership.total_set() ) );
	}
};

boost::python::list
top_mers_by_sequence_membership( boost::python::object py_seqs, unsigned k, unsigned n )
{
	typedef mer_count_with_seq_membership count;
	rev_comp_collapser< count > collapser( k );
	count_mers< count >( make_boost_range< const output_seq & >( py_seqs ), k, collapser );
	top_mers_callback callback( k );
	collapser.largest_counts( n, callback );
	return callback._k_mers;
}


void
set_emission_distribution( model & m, unsigned i, boost::python::object emission_dist_python )
{
	if( boost::python::len( emission_dist_python ) != int( m.M ) ) throw std::logic_error( "emission distribution must have model.M entries" );
	double_vec emission_dist_C;
	myrrh::python::convert_from_python( emission_dist_python, emission_dist_C );
	for( unsigned k = 0; 4 != k; ++k ) m.theta[ *( m.b_parameterisation[ i ][ k ] ) ] = emission_dist_C[ k ];
}

void seed_rng( unsigned long int s )
{
	gsl_rng_set( myrrh::gsl::get_rng(), s );
}


bool optional_param_idx_has_value( optional_param_idx idx ) { return idx; }
unsigned optional_param_idx_get_value( optional_param_idx idx ) {
	if( ! idx ) throw std::logic_error( "This optional parameter index has no value" );
	return *idx;
}
void optional_param_idx_set_value( optional_param_idx & idx, boost::python::object o )
{
	if( ! o ) idx = optional_param_idx();
	else idx = unsigned( boost::python::extract< unsigned >( o ) );
}
std::string optional_param_idx_str( optional_param_idx idx ) { if( idx ) return HMM_MAKE_STRING( *idx ); else return "<*>"; }

model::ptr
deserialise_model( const std::string & filename, bool binary ) {
	model::ptr m( new model );
	const boost::filesystem::path archive_filename( filename );
	if( binary ) {
		myrrh::deserialise< true >( m, archive_filename );
	} else {
		myrrh::deserialise< false >( m, archive_filename );
	}
	return m;
}

void
serialise_model( model::ptr m, const std::string & filename, bool binary ) {
	const boost::filesystem::path archive_filename( filename );
	if( binary ) {
		myrrh::serialise< true >( m, archive_filename );
	} else {
		myrrh::serialise< false >( m, archive_filename );
	}
}

HMM_EXPORT_FN_SPEC
void
export_model()
{
	using namespace boost::python;
	using namespace indexing;
	using namespace myrrh::python;

	//std::cout << "Exposing converters\n";
	import_array(); //make sure we have initialised numpy API
	expose_converters< npy_byte >();
	expose_converters< npy_ubyte >();
	expose_converters< npy_short >();
	expose_converters< npy_ushort >();
	expose_converters< npy_int >();
	expose_converters< npy_uint >();
	expose_converters< npy_long >();
	expose_converters< npy_ulong >();
	expose_converters< npy_longlong >();
	expose_converters< npy_ulonglong >();

	myrrh::python::register_tuple_converter< training_result_tuple >();
	//myrrh::python::register_converter< double_vec_ptr >();
	//numpy_multi_array_converter< double_array >::register_to_and_from_python();

	def( "seed_rng", seed_rng, "Seed the random number generator to make results reproducible" );
	def( "is_close", myrrh::is_close< double > );

	class_< double_vec, double_vec_ptr, boost::noncopyable >( "DoubleVec" )
	.def( container_suite< double_vec >() )
	.def_pickle( vector_pickle_suite< double_vec >() )
	;

	class_<
		output_seq,
		impl::output_seq_shared_ptr,
		boost::noncopyable
	>(
		"PreprocessedSequence",
		"Type that holds a pre-processed sequence to cut down on copying between language boundary.",
		no_init
	)
	.def(
		"as_numpy",
		impl::as_writeable_numpy< hmm::output_seq >,
		"Copy the sequence data into a numpy array."
	);
	def(
		"preprocess_sequence",
		impl::convert_sequence,
		"Pre-processes a sequence to cut down on copying between language boundary"
	);

	class_<
		optional_param_idx
	>(
		"OptionalParameterIndex",
		"Optional index for the parameterisation of one of the probabilities in a HMM. E.g. a = theta[ idx ]",
		init< unsigned >( "Create the optional parameter index with the given index value" )
	)
	.def(
		init<>( "Create the optional parameter index in an uninitialised state (i.e. no index)" ) )
	.def_pickle( optional_pickle_suite< unsigned >() )
	.def(
		"__nonzero__",
		optional_param_idx_has_value,
		"Does the parameterisation have a value?"
	)
	.def(
		"__str__",
		optional_param_idx_str,
		"String representation of the parameterisation"
	)
	.add_property(
		"idx",
		optional_param_idx_get_value,
		optional_param_idx_set_value,
		"The parameterisation's index into the parameters"
	)
	;

	class_<
		param_idx_vec
	>(
		"OptionalParameterIndexVec" )
	.def( container_suite< param_idx_vec >() )
	.def_pickle( vector_pickle_suite< param_idx_vec >() )
	;
	numpy_multi_array_converter< param_idx_array >::register_to_and_from_python();

	using boost::python::arg;
	class_<
		model
	>(
		"Model",
		"HMM",
		init< unsigned, unsigned, unsigned, unsigned >( ( arg("N"), arg("M"), arg("P"), arg("markov_order") = 0 ) )
	)
	.def( init< const model & >( arg("model"), "Copy constructs the model from another" ) )
	.def_pickle( model_pickle_suite() )
	.add_property(
		"converter",
		wrapped_converter,
		"Converts between order 0 emissions and order n emissions" )
	.def_readonly(
		"N",
		&model::N,
		"# states in the model" )
	.def_readonly(
		"M",
		&model::M,
		"# characters in the output alphabet of the model" )
	.def_readonly(
		"P",
		&model::P,
		"# parameters in the model" )
	.add_property(
		"A",
		transition_matrix,
		"Matrix of all transition probabilities, A[i,j] = p(transition from state i to state j)" )
	.add_property(
		"B",
		emission_matrix,
		"Matrix of all emission probabilities, A[i,k] = p(output symbol k from state i)" )
	.def(
		"a",
		&model::a_checked,
		"State transition distributions, p(transition from state i to state j)",
		( arg( "i" ), arg( "j" ) ) )
	.def(
		"b",
		&model::b_checked,
		"Emission distributions, p(output symbol k from state i)",
		( arg( "i" ), arg( "k" ) ) )
	.def(
		"set_transition",
		&model::set_transition,
		"Set the transition probability, p, from state i to j",
		( arg("i"), arg("j"), arg("p") ) )
	.def(
		"set_emission",
		&model::set_emission,
		"Set the emission probability, p, of output k from state i",
		( arg("i"), arg("k"), arg("p") ) )
	.def(
		"set_initial",
		&model::set_initial,
		"Set the probability, p, of initially being in state i",
		( arg("i"), arg("p") ) )
	.def(
		"get_transition_parameterisation",
		get_transition_parameterisation,
		"The parameterisation for the transition from state i to state j",
		( arg( "i" ), arg( "j" ) ) )
	.def(
		"set_transition_parameterisation",
		set_transition_parameterisation,
		"The parameterisation for the transition from state i to state j",
		( arg( "i" ), arg( "j" ) ) )
	.def(
		"get_emission_parameterisation",
		get_emission_parameterisation,
		"The parameterisation for the emission of symbol k from state i",
		( arg( "i" ), arg( "k" ) ) )
	.def(
		"set_emission_parameterisation",
		set_emission_parameterisation,
		"The parameterisation for the emission of symbol k from state i",
		( arg( "i" ), arg( "k" ) ) )
	.def(
		"set_emission_distribution",
		set_emission_distribution,
		"Set the emission distribution for state i",
		( arg( "i" ), arg( "emission_dist" ) ) )
	.def_readwrite(
		"theta",
		&model::theta,
		"The model's parameters" )
	.add_property(
		"pi",
		initial_distribution,
		"The distribution over initial states" )
	.def_readwrite(
		"pi_parameterisation",
		&model::pi_parameterisation,
		"The parameterisation for the distribution over initial states" )
	.def(
		"_check_consistent",
		check_consistent,
		"Make sure model is consistent" )
	.def(
		"normalise",
		normalise_model_parameters,
		"Scale parameters of model to make proper distributions" )
	.def(
		"viterbi",
		wrapped_viterbi,
		"The viterbi algorithm - returns log p( s ), s",
		( arg( "model" ), arg( "observations" ) ) )
	.def(
		"viterbi_2",
		wrapped_viterbi_2,
		"The viterbi algorithm - returns log p( s ), s",
		( arg( "model" ), arg( "observations" ) ) )
	.def(
		"LL",
		wrapped_LL,
		"Returns the log likelihood of the observations (one sequence)",
		( arg( "model" ), arg( "observations" ), arg( "scaling_threshold" ) = 1e-8 ) )
	.def(
		"LL_multi",
		wrapped_LL_multi,
		"Returns the log likelihood of the sequences",
		( arg( "model" ), arg( "sequences" ), arg( "scaling_threshold" ) = 1e-8 ) )
	.def(
		"forward",
		wrapped_forward,
		"The forward procedure - returns log(obs), alpha, c",
		( arg( "model" ), arg( "observations" ), arg( "scaling_threshold" ) = 1e-8 ) )
	.def(
		"backward",
		wrapped_backward,
		"The backward procedure - returns beta",
		( arg( "model" ), arg( "observations" ), arg( "scaling_factors" ) ) )
	.def(
		"forward_backward",
		wrapped_forward_backward,
		"Runs the forward and backward procedure - returns (LL, alpha, beta, c)",
		( arg( "model" ), arg( "observations" ), arg( "scaling_threshold" ) = 1e-8 ) )
	.def(
		"baum_welch_iteration",
		wrapped_baum_welch_iteration,
		"One iteration of the Baum-Welch procedure",
		(
			arg( "model" ),
			arg( "sequences" ),
			arg( "prior" ) = boost::python::object()
		) )
	.def(
		"baum_welch",
		wrapped_baum_welch,
		"The Baum-Welch procedure - returns (LL, # iterations)",
		(
			arg( "model" ),
			arg( "sequences" ),
			arg( "prior" ) = boost::python::object(),
			arg( "max_iterations" ) = 0,
			arg( "tolerance" ) = 1e-5,
			arg( "callback" ) = boost::python::object()
		) )
	.def(
		"viterbi_training",
		wrapped_viterbi_training,
		"The Viterbi training procedure - returns (LL, # iterations)",
		(
			arg( "model" ),
			arg( "sequences" ),
			arg( "prior" ) = boost::python::object(),
			arg( "max_iterations" ) = 0,
			arg( "tolerance" ) = 1e-5,
			arg( "callback" ) = boost::python::object()
		) )
	.def(
		"sample",
		wrapped_sample,
		"Sample one sequence of length T from the model",
		( arg( "model" ), arg( "T" ) ) )
	.def(
		"serialise",
		serialise_model,
		"Serialise the model to the given filename in a binary or text format.",
		( arg( "model" ), arg( "filename" ), arg( "binary" ) = true ) )
	.def(
		"deserialise",
		deserialise_model,
		"Deserialise the model from the given filename in a binary or text format.",
		( arg( "filename" ), arg( "binary" ) = true ) )
	.staticmethod( "deserialise" )
	;
	register_ptr_to_python< model::ptr >();

	class_<
		model::prior
	>(
		"ModelPrior",
		"Prior for the parameters of a HMM",
		init< unsigned, unsigned >( ( arg("N"), arg("M") ) )
	)
	.add_property(
		"A",
		access_and_convert< model::prior, double_array, &model::prior::a >,
		convert_and_set< model::prior, double_array, &model::prior::a >,
		"Prior for the model's transitions" )
	.add_property(
		"B",
		access_and_convert< model::prior, double_array, &model::prior::b >,
		convert_and_set< model::prior, double_array, &model::prior::b >,
		"Prior for the model's emissions" )
	.add_property(
		"pi",
		access_and_convert< model::prior, double_vec, &model::prior::pi >,
		convert_and_set< model::prior, double_vec, &model::prior::pi >,
		"Prior for the model's initial state distribution" )
	.add_property(
		"N",
		&model::prior::N,
		"# states in the model" )
	.add_property(
		"M",
		&model::prior::M,
		"# characters in the output alphabet" )
	.def(
		"_check_matches_model",
		&model::prior::check_matches_model,
		"Check the shape of the prior matches the model",
		arg( "model" ) )
	.def(
		"_matches_model",
		&model::prior::matches_model,
		"Does the shape of the prior match the model?",
		arg( "model" ) )
	;

	def(
		"build_pssm_model",
		wrapped_build_pssm_model,
		"Build a HMM that models binding sites for a simple PSSM",
		(
			arg( "p_binding_site" ),
			arg( "position_distributions" )
		) );

	def(
		"shift_pssm_model",
		wrapped_shift_pssm_model,
		"Copies and shifts a HMM that models a pssm +/- shift_offset bases.",
		(
			arg( "original_model" ),
			arg( "position_distributions" ),
			arg( "fill_in_distribution" )
		) );

	def(
		"count_mers",
		wrapped_count_mers,
		"Counts all the mers of a given length, n.",
		(
			arg( "sequences" ),
			arg( "n" ),
			arg( "callback" )
		) );

	def(
		"count_mers_with_sequence_counts",
		wrapped_count_mers_with_sequence,
		"Counts all the mers of a given length, n. Include how many sequences each appears in.",
		(
			arg( "sequences" ),
			arg( "n" ),
			arg( "callback" )
		) );

	def(
		"top_mers_by_sequence_membership",
		top_mers_by_sequence_membership,
		"Counts all k-mers in all the sequences, sorts them by # sequences they occur in then by total occurences and returns the top n.",
		(
			arg( "sequences" ),
			arg( "k" ),
			arg( "n" )
		) );

	def(
		"build_fully_connected_hmm",
		create_random_fully_connected_model,
		"Build a HMM that is fully connected and has no tied parameters",
		( arg( "N" ), arg( "M" ), arg( "markov_order" ) ) );

	def(
		"dirichlet_draw",
		wrapped_draw_from_dirichlet,
		"Draws from a dirichlet distribution" );
}

} //namespace hmm
