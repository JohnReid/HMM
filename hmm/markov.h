#ifndef HMM_MARKOV_H_
#define HMM_MARKOV_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER


#include "hmm/defs.h"

#include <myrrh/types.h>
#include <myrrh/math.h>
#include <myrrh/serialize/multi_array.h>

namespace hmm {


/// Probabilities on usual scale [0,1]
template< typename Float >
struct probability_scale {
    typedef Float float_t;
    static Float identity() { return Float( 1. ); }
    static Float multiply( Float p1, Float p2 ) { return p1 * p2; }
};



/// Probabilities on log scale
template< typename Float >
struct log_probability_scale {
    typedef Float float_t;
    static Float identity() { return Float( 0. ); }
    static Float multiply( Float p1, Float p2 ) { return p1 + p2; }
    static Float from_probability( Float p ) { return std::log( p ); }
    static Float to_probability( Float p ) { return std::exp( p ); }
};




//forward decl
template<
    size_t Order,
    size_t AlphabetSize = 4,
    typename Float = double,
    typename Scale = log_probability_scale< Float >
> struct complete_markov_model;





namespace detail {

/// Null model. Used as specialisation for the lower order model of a 0-order model.
template< typename Scale >
struct null_model {
    template< typename It >
    typename Scale::float_t evaluate( It begin, It end ) const {
        return Scale::identity();
    }

    template< typename InputIt, typename OutputIt >
    void
    _consume_lower_orders(
        InputIt begin,
        typename std::iterator_traits< InputIt >::difference_type length,
        OutputIt output ) const
    {
    }


    //
    // Serialization
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) { }

    template< typename InputIt, typename OutputIt >
    InputIt
    evaluate_in_context(
        InputIt context,
        InputIt begin,
        InputIt end,
        OutputIt output
    ) const {
        return begin;
    }
};

/// Select lower order Markov model.
template<
    size_t Order,
    size_t AlphabetSize,
    typename Float,
    typename Scale
>
struct lower_order_selector {
    typedef complete_markov_model< Order - 1, AlphabetSize, Float, Scale > type;
};

/// Select lower order Markov model (specialisation).
template< size_t AlphabetSize, typename Float, typename Scale >
struct lower_order_selector< 0, AlphabetSize, Float, Scale > {
    typedef null_model< Scale > type;
};

/** Traits of the Markov model. */
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct markov_traits {
    static size_t order() { return Order; } ///< The order of the Markov model.
    static size_t alphabet_size() { return AlphabetSize; } ///< The alphabet size.
    typedef Float float_t; ///< The type used to store (log) probabilities.
    typedef Scale scale_t; ///< The type that defines the scale we store probabilities on.
};


/** Calculates extents for boost multi array. */
template< size_t NumDims, size_t AlphabetSize >
struct extents_builder {
    typedef typename boost::detail::multi_array::extent_gen< NumDims > return_t; /**< Functor return type. */

    /** Calculates extents for boost multi array. */
    return_t operator()() const {
        return extents_builder< NumDims-1, AlphabetSize >()()[AlphabetSize];
    }
};

//specialise when AlphabetSize = 0
/** Calculates extents for boost multi array. */
template< size_t AlphabetSize >
struct extents_builder< 0, AlphabetSize > {
    typedef boost::detail::multi_array::extent_gen< 0 > return_t; /**< Functor return type. */

    /** Calculates extents for boost multi array. */
    return_t operator()() const {
        return boost::extents;
    }
};



/// Access an element of a multi_array indexed by the values in the iterator.
template< typename Reference, typename MultiArray, typename IndexIt >
Reference
access_multi_array_element(
    boost::type< Reference >,
    MultiArray & a,
    IndexIt it
) {
    typename MultiArray::index offset = 0;
    const typename MultiArray::index * strides = a.strides();
    for( typename MultiArray::size_type n = 0; n != MultiArray::dimensionality; ++n ) {
        MYRRH_ASSERT( *it < a.shape()[n] );
        offset += ( unsigned int )( *it ) * strides[n];
        ++it;
    }

    return a.origin()[offset];
}


/// Initialise counts for a complete Markov model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct
add_pseudo_counts {
    void
    operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & model, Float pseudo_count ) const {
        std::for_each(
            model.mm.storage.data(),
            model.mm.storage.data() + model.mm.storage.num_elements(),
            boost::lambda::_1 += pseudo_count
        );
        add_pseudo_counts< Order-1, AlphabetSize, Float, Scale >()( model.lower_orders, pseudo_count );
    }
};


/// Specialisation for 0-order
template< size_t AlphabetSize, typename Float, typename Scale >
struct
add_pseudo_counts< 0, AlphabetSize, Float, Scale > {
    void
    operator()( complete_markov_model< 0, AlphabetSize, Float, Scale > & model, Float pseudo_count ) const {
        std::for_each(
            model.mm.storage.data(),
            model.mm.storage.data() + model.mm.storage.num_elements(),
            boost::lambda::_1 += pseudo_count
        );
    }
};



/// Add the Markov counts for the sequence to the counts array.
template< typename InputIt, typename CountsArray >
void
add_counts( CountsArray & counts, InputIt begin, InputIt end ) {
    typedef typename CountsArray::element element;
    if( CountsArray::dimensionality < 1 ) {
        throw std::logic_error( "Need to have dimensionality at least 1." );
    }
    size_t order = CountsArray::dimensionality - 1;
    if( end - begin <= typename std::iterator_traits< InputIt >::difference_type( order ) ) {
        return; //sequence not long enough, nothing to do
    }
    const InputIt last = end - order;
    for( InputIt i = begin; last != i; ++i ) {
        detail::access_multi_array_element( boost::type< element & >(), counts, i ) += element( 1 );
    }
}


/// Add counts to complete Markov model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct
add_counts_to_complete {
    template< typename It >
    void
    operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & model, It begin, It end ) const {
        add_counts( model.mm.storage, begin, end );
        add_counts_to_complete< Order-1, AlphabetSize, Float, Scale >()( model.lower_orders, begin, end );
    }
};


/// Specialisation for 0-order
template< size_t AlphabetSize, typename Float, typename Scale >
struct
add_counts_to_complete< 0, AlphabetSize, Float, Scale > {
    template< typename It >
    void
    operator()( complete_markov_model< 0, AlphabetSize, Float, Scale > & model, It begin, It end ) const {
        add_counts( model.mm.storage, begin, end );
    }
};



/// Normalise a range. That is, scale values so that sum == 1.
template< typename It >
void
normalise( It begin, It end ) {
    typedef typename It::value_type value_type;
    const value_type sum = std::accumulate( begin, end, value_type( 0 ) );
    for( ; end != begin; ++begin ) {
        *begin /= sum;
    }
}

/// Normalise a subarray view
template<
    typename ValueType,
    std::size_t NumDims
>
struct
normalise_subarray {
    typedef boost::detail::multi_array::sub_array< ValueType, NumDims > sub_array;

    void
    operator()( sub_array a ) const {
        for( typename sub_array::iterator i = a.begin(); a.end() != i; ++i ) {
//      typedef typename multi_array_t::template array_view< Dims-1 >::type subarray_view;
//      BOOST_FOREACH( subarray_view sa, a ) {
            normalise_subarray< ValueType, NumDims-1 >()( *i );
        }
    }
};

/// Specialisation for 1 dimension
template <
    typename ValueType
>
struct
normalise_subarray< ValueType, 1 > {
    typedef boost::detail::multi_array::sub_array< ValueType, 1 > sub_array;

    void
    operator()( sub_array a ) const {
        normalise( a.begin(), a.end() );
    }
};


/// Normalise a multi_array
template<
    typename ValueType,
    std::size_t NumDims,
    typename Allocator
>
struct
normalise_multi_array {
    typedef boost::multi_array< ValueType, NumDims, Allocator > multi_array_t;

    void
    operator()( multi_array_t & a ) const {
        for( typename multi_array_t::iterator i = a.begin(); a.end() != i; ++i ) {
            normalise_subarray< ValueType, NumDims-1 >()( *i );
        }
    }
};

/// Specialisation for 1 dimension
template <
    typename ValueType,
    typename Allocator
>
struct
normalise_multi_array< ValueType, 1, Allocator > {
    void
    operator()( boost::multi_array< ValueType, 1, Allocator > & a ) const {
        normalise( a.begin(), a.end() );
    }
};


/// Normalise counts in complete Markov model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct
normalise_counts {
    void
    operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & model ) const {
        normalise_multi_array< Float, Order+1, std::allocator< Float > >()( model.mm.storage );
        normalise_counts< Order-1, AlphabetSize, Float, Scale >()( model.lower_orders );
    }
};


/// Specialisation for 0-order
template< size_t AlphabetSize, typename Float, typename Scale >
struct
normalise_counts< 0, AlphabetSize, Float, Scale > {
    void
    operator()( complete_markov_model< 0, AlphabetSize, Float, Scale > & model ) const {
        normalise( model.mm.storage.begin(), model.mm.storage.end() );
    }
};



/// Scale the array
template< typename Scale, typename Array >
void
scale_array( Array & a ) {
    for( typename Array::element * i = a.data(); a.data() + a.num_elements() != i; ++i ) {
        *i = Scale::from_probability( *i );
    }
}


/// Convert complete Markov model to correct scale.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct
convert_to_scale {
    void
    operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & model ) const {
        scale_array< Scale >( model.mm.storage );
        convert_to_scale< Order-1, AlphabetSize, Float, Scale >()( model.lower_orders );
    }
};


/// Specialisation for 0-order
template< size_t AlphabetSize, typename Float, typename Scale >
struct
convert_to_scale< 0, AlphabetSize, Float, Scale > {
    void
    operator()( complete_markov_model< 0, AlphabetSize, Float, Scale > & model ) const {
        scale_array< Scale >( model.mm.storage );
    }
};


/// Convert probabilities to scale.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct
convert_to_scale_if_needed {
    void
    operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & model ) {
        detail::convert_to_scale< Order, AlphabetSize, Float, Scale >()( model );
    }
};

/// Specialisation for probability scale.
template< size_t Order, size_t AlphabetSize, typename Float >
struct
convert_to_scale_if_needed< Order, AlphabetSize, Float, probability_scale< Float > > {
    void
    operator()( complete_markov_model< Order, AlphabetSize, Float, probability_scale< Float > > & model ) {
        //do nothing
    }
};



} //namespace detail




/** A Markov model that models characters depending on Order previous characters. */
template<
    size_t Order,
    size_t AlphabetSize = 4,
    typename Float = double,
    typename Scale = log_probability_scale< Float >
>
struct markov_model {

    typedef detail::markov_traits< Order, AlphabetSize, Float, Scale > traits;
    typedef boost::multi_array< Float, Order+1 > storage_t; ///< Stores the (log) probabilities.

    /** The extents of the counts boost multi_array. */
    static typename detail::extents_builder< Order+1, AlphabetSize >::return_t get_extents() {
        return detail::extents_builder< Order+1, AlphabetSize >()();
    }

    //
    // Members
    storage_t storage;

    //
    // Serialization
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template< typename Archive >
    void serialize(Archive & ar, const unsigned int version) {
        ar & storage;
    }


    markov_model( Float init_value = Float( 0. ) ) : storage( detail::extents_builder< Order+1, AlphabetSize >()() ) {
        std::fill( storage.data(), storage.data() + storage.num_elements(), init_value );
    }

    template< typename It >
    float_t evaluate( It begin ) const {
        return detail::access_multi_array_element( boost::type< float_t >(), storage, begin );
    }
};



// forward decl.
template< size_t LowerOrder, size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct _get_lower_order_model;


/** A Markov model that models characters depending on as many preceding characters as are available. */
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct complete_markov_model {

    //
    // Typedefs
    typedef complete_markov_model< Order, AlphabetSize, Float, Scale > self_t;
    typedef boost::shared_ptr< self_t >                                ptr;
    typedef detail::markov_traits< Order, AlphabetSize, Float, Scale > traits;

    //
    // Members
    markov_model< Order, AlphabetSize, Float, Scale > mm;
    typename detail::lower_order_selector< Order, AlphabetSize, Float, Scale >::type lower_orders;

    //
    // Serialization
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template< typename Archive >
    void serialize(Archive & ar, const unsigned int version) {
        ar & mm;
        ar & lower_orders;
    }

    /// Evaluate the sequence given by [begin, end) and output cumulative probabilities into the output iterator
    template< typename InputIt, typename OutputIt >
    void
    evaluate( InputIt begin, InputIt end, OutputIt output ) const {
        const typename std::iterator_traits< InputIt >::difference_type length = end - begin;
        const size_t order = traits::order();
        _consume_lower_orders( begin, length, output );
        if( order < size_t( length ) ) {
            end -= order;
            for( ++begin; begin != end; ++begin ) {
                const double likelihood = mm.evaluate( begin );
                *output = likelihood;
                ++output;
            }
        }
    }

    /// Evaluate the sequence given by [begin, end) but preceded by data starting at context.
    /// Output cumulative probabilities into the output iterator
    template< typename InputIt, typename OutputIt >
    InputIt
    evaluate_in_context(
        InputIt context,
        InputIt begin,
        InputIt end,
        OutputIt output
    ) const {
        // if we cannot evaluate the beginning at this order, go to a lower order
        if( context + traits::order() > begin) {
            begin = lower_orders.evaluate_in_context( context, begin, std::min( end, context + traits::order() ), output );
        }

        // evaluate what we can at this order
        begin -= traits::order();
        end -= traits::order();
        while( begin < end ) {
            const double likelihood = mm.evaluate( begin );
            *output = likelihood;
            ++output;
            ++begin;
        }

        return begin + traits::order();
    }

    template< size_t LowerOrder >
    complete_markov_model< LowerOrder, AlphabetSize, Float, Scale > &
    get_lower_order_model() {
        return _get_lower_order_model< LowerOrder, Order, AlphabetSize, Float, Scale >()( *this );
    }

    template< typename InputIt, typename OutputIt >
    void
    _consume_lower_orders(
        InputIt begin,
        typename std::iterator_traits< InputIt >::difference_type length,
        OutputIt output ) const
    {
        lower_orders._consume_lower_orders( begin, length, output );
        if( typename std::iterator_traits< InputIt >::difference_type( traits::order() ) < length ) {
            const double likelihood = mm.evaluate( begin );
            *output = likelihood;
            ++output;
        }
    }
};

/// Get the lower order model.
template< size_t LowerOrder, size_t Order, size_t AlphabetSize, typename Float, typename Scale >
struct _get_lower_order_model {
    complete_markov_model< LowerOrder, AlphabetSize, Float, Scale > &
    operator()( complete_markov_model< Order, AlphabetSize, Float, Scale > & self ) {
        return self.lower_orders.template get_lower_order_model< LowerOrder >();
    }
};

/// Specialise for LowerOrder == Order
template< size_t LowerOrder, size_t AlphabetSize, typename Float, typename Scale >
struct _get_lower_order_model< LowerOrder, LowerOrder, AlphabetSize, Float, Scale > {
    complete_markov_model< LowerOrder, AlphabetSize, Float, Scale > &
    operator()( complete_markov_model< LowerOrder, AlphabetSize, Float, Scale > & self ) {
        return self;
    }
};


/// Initialise every count in the Markov model's storage.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
void
add_pseudo_counts( complete_markov_model< Order, AlphabetSize, Float, Scale > & model, Float init_value ) {
    detail::add_pseudo_counts< Order, AlphabetSize, Float, Scale >()( model, init_value );
}


/// Add counts from sequence to complete model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale, typename It >
void
add_counts_to_complete( complete_markov_model< Order, AlphabetSize, Float, Scale > & model, It begin, It end ) {
    detail::add_counts_to_complete< Order, AlphabetSize, Float, Scale >()( model, begin, end );
}


/// Normalise counts in complete model.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
void
normalise_counts( complete_markov_model< Order, AlphabetSize, Float, Scale > & model ) {
    detail::normalise_counts< Order, AlphabetSize, Float, Scale >()( model );
}


/// Convert probabilities to scale.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale >
void
convert_to_scale( complete_markov_model< Order, AlphabetSize, Float, Scale > & model ) {
    detail::convert_to_scale_if_needed< Order, AlphabetSize, Float, Scale >()( model );
}


/// Initialise the model from the sequence.
template< size_t Order, size_t AlphabetSize, typename Float, typename Scale, typename It >
void
initialise_from_sequence(
    complete_markov_model< Order, AlphabetSize, Float, Scale > & model,
    Float pseudo_count,
    It begin,
    It end
) {
    add_pseudo_counts( model, pseudo_count );
    add_counts_to_complete( model, begin, end );
    normalise_counts( model );
    convert_to_scale( model );
}


/**
 * Functor that takes values and adds them to the back of the seq.
 */
template< typename BackInsertionSeq >
struct push_back_cumulative {
    BackInsertionSeq & seq;
    push_back_cumulative( BackInsertionSeq & seq ) : seq( seq ) { }
    void operator()( double x ) {
        seq.push_back( seq.empty() ? x : x + seq.back() );
    }
};


/// Push back cumulative
template< typename BackInsertionSeq >
push_back_cumulative< BackInsertionSeq >
make_push_back_cumulative( BackInsertionSeq & seq ) {
    return push_back_cumulative< BackInsertionSeq >( seq );
}




} //namespace hmm

#endif //HMM_MARKOV_H_
