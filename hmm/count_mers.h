/**
@file

Copyright John Reid 2006

*/

#ifndef HMM_COUNT_MERS_H_
#define HMM_COUNT_MERS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER

#include "hmm/defs.h"
#include "hmm/model.h"

#include "myrrh/math.h"
#include "myrrh/bit_set.h"

#include <queue>

namespace hmm {
namespace impl {

template< typename OutputFn >
struct count_mers_output_fn
{
	OutputFn & output_fn;
	
	count_mers_output_fn( OutputFn & output_fn ) : output_fn( output_fn ) { }

	template< typename NodePtr, typename Prefix >
	void operator()( NodePtr node, const Prefix & prefix )
	{
		if( node->value.total ) output_fn( prefix, node->value );
	}
};

} //namespace impl





struct mer_count {
	unsigned total; /**< Total count of all occurences in the sequences. */

	mer_count( unsigned num_seqs ) : total( 0 ) { }

	void operator()( unsigned seq_idx ) {
		++total;
	}

	bool operator<( const mer_count & rhs ) const {
		return total < rhs.total;
	}

	mer_count & operator|=( const mer_count & rhs ) {
		total += rhs.total;
		return *this;
	}
};


struct mer_count_with_seq_membership : mer_count {
	myrrh::bit_set<> seq_idx_membership; /**< Which sequences this string belongs to. */

	mer_count_with_seq_membership( unsigned num_seqs ) : mer_count( num_seqs ), seq_idx_membership( num_seqs ) { }

	void operator()( unsigned seq_idx ) {
		mer_count::operator()( seq_idx );
		seq_idx_membership.set( seq_idx );
	}

	bool operator<( const mer_count_with_seq_membership & rhs ) const {
		unsigned num_seqs = seq_idx_membership.total_set();
		unsigned rhs_num_seqs = rhs.seq_idx_membership.total_set();
		if( num_seqs == rhs_num_seqs ) return mer_count::operator<(rhs);
		return num_seqs < rhs_num_seqs;
	}

	mer_count_with_seq_membership & operator|=( const mer_count_with_seq_membership & rhs ) {
		mer_count::operator|=( rhs );
		seq_idx_membership |= rhs.seq_idx_membership;
		return *this;
	}
};





template< 
	typename CountValue,
	typename SeqRange,
	typename OutputFn
>
void
count_mers(
	const SeqRange & sequences,
	unsigned n,
	OutputFn & output )
{
	using namespace boost;

	typedef typename range_value< SeqRange >::type sequence_t;
	typedef typename range_value< sequence_t >::type char_t;

	typedef myrrh::trie::trie< char_t, CountValue > trie;
	typedef typename trie::node_ptr node_ptr;
	trie t;

	unsigned num_seqs = boost::size( sequences );
	unsigned seq_idx = 0;
	BOOST_FOREACH( const sequence_t seq, sequences )
	{
		typedef typename range_const_iterator< sequence_t >::type const_iterator;
		const_iterator begin = const_begin( seq );
		const_iterator end = const_end( seq );
		const_iterator mer_end = begin;
		for( unsigned i = 0; n != i && end != mer_end; ++i ) ++mer_end;
		if( mer_end - begin == int( n ) )
		{
			while( true )
			{
				node_ptr p = t.insert( begin, mer_end, CountValue( num_seqs ) );
				p->value( seq_idx );
				if( end == mer_end ) break;
				++begin;
				++mer_end;
			}
		}
		++seq_idx;
	}

	impl::count_mers_output_fn< OutputFn > output_fn( output );
	t.traverse_prefixes( output_fn );
}







/**
Gets the complement of a nucleotide.
*/
struct nucleo_complement
	: std::unary_function< char, char >
{
	inline char operator()( char c ) const {
		if( 4 == c ) return 4;
		if( 4 < c || 0 > c ) throw std::logic_error( HMM_MAKE_STRING( "Character must be between 0 and 4: '" << c << "'" ) );
		return 3 - c;
	}
};



/**
Returns an iterator that iterates over the complement of the original iterator.
*/
template< typename Iterator >
boost::transform_iterator< nucleo_complement, Iterator >
make_complement_iterator( Iterator it )
{
	return boost::make_transform_iterator( it, nucleo_complement() );
}





/**
Returns an iterator that iterates over the reverse complement of the original iterator.
*/
template< typename BidirectionalIterator >
boost::reverse_iterator< boost::transform_iterator< nucleo_complement, BidirectionalIterator > >
make_reverse_complement_iterator( BidirectionalIterator x )
{
	return boost::make_reverse_iterator( make_complement_iterator( x ) );
}



/**
Defines the type of a reverse complement iterator.
*/
template< typename BidirectionalIterator >
struct reverse_complement_iterator {
	typedef boost::reverse_iterator< boost::transform_iterator< nucleo_complement, BidirectionalIterator > > type;
};




/**
Enables us to count k-mers collapsed with their reverse complements.
*/
template< typename Value >
struct rev_comp_collapser
{
	typedef Value value_t;
	typedef markov_order_n_converter<> converter_t;
	typedef std::map< converter_t::char_t, value_t > counts_t;
	typedef typename counts_t::iterator counts_it;
	typedef typename counts_t::const_iterator counts_const_it;

	converter_t _converter;
	counts_t _counts;

	rev_comp_collapser( unsigned n ) : _converter( 4, n-1 ) { }

	template< typename Prefix >
	void
	operator()( const Prefix & prefix, const Value & value ) 
	{
		HMM_VERIFY( boost::size( prefix ) == _converter.order+1 );

		//calculate a number that represents the reverse complement of the prefix
		const converter_t::char_t mer_idx = _converter.convert_to_order_n(
			boost::begin( prefix ),
			boost::end( prefix ) );

		//look for this prefix in the map
		counts_it count = _counts.find( mer_idx );

		//did we find it?
		if( _counts.end() == count ) {
			// we have not seen the reverse complement of our prefix so insert our counts as new entry
			const converter_t::char_t rev_comp_idx = _converter.convert_to_order_n(
				make_reverse_complement_iterator( boost::end( prefix ) ),
				make_reverse_complement_iterator( boost::begin( prefix ) ) );
			_counts.insert( typename counts_t::value_type( rev_comp_idx, value ) );
		} else {
			//update the counts we have for the reverse complement of this prefix with these counts.
			count->second |= value;
		}
	}

	struct compare_second {
		bool operator()( counts_const_it x, counts_const_it y ) const {
			return x->second < y->second;
		}
	};

	/**
	Output the prefixes with the largest counts.
	*/
	template< typename OutputFn >
	void
	largest_counts( unsigned num_to_output, OutputFn output_fn )
	{
		//copy all the iterators into a priority queue sorted by count
		std::priority_queue< counts_const_it, std::vector< counts_const_it >, compare_second > Q;
		counts_const_it end = _counts.end();
		for( counts_const_it i = _counts.begin(); end != i; ++i ) {
			Q.push( i );
		}

		//output the top so many...
		while( num_to_output && ! Q.empty() ) {
			output_fn( Q.top()->first, Q.top()->second );
			Q.pop();
			--num_to_output;
		}
	}
};





} //namespace hmm

#endif //HMM_COUNT_MERS_H_

