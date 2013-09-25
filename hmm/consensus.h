/**
@file

Copyright John Reid 2006

*/

#ifndef HMM_CONSENSUS_H_
#define HMM_CONSENSUS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include "hmm/defs.h"

#include <myrrh/math.h>

namespace hmm {
namespace consensus {

/** Maps from consensus code indices to which characters in original alphabet the codes represent. */
template< typename Char = unsigned, typename Code = unsigned >
struct codes
{
	typedef Char char_t;
	typedef Code code_t;

	char_t _alphabet_size;
	code_t _num_codes;

	codes( char_t alphabet_size )
		: _alphabet_size( alphabet_size )
		, _num_codes( myrrh::int_power(2, alphabet_size) )
	{
	}

	void check_char(char_t c) const
	{
		if(c >= _alphabet_size)
			throw std::logic_error("Character not in alphabet");
	}
	void check_code(code_t c) const {
		if(c >= _num_codes)
			throw std::logic_error("Invalid code");
	}
	char_t alphabet_size() const { return _alphabet_size; }
	code_t num_codes() const { return _num_codes; }
	bool matches(char_t c, code_t code) const
	{
		check_char(c);
		check_code(code);
		return bool(1 & (code) >> c);
	}
	code_t code_for_char(char_t c) const
	{
		check_char(c);
		const code_t result = (1 << c);
		check_code(result);
		return result;
	}

	template< typename CharRange >
	code_t
	code_for_chars( const CharRange & chars ) const
	{
		code_t result = 0;
		BOOST_FOREACH( char_t c, chars )
		{
			check_char(c);
			result += (1 << c);
		}
		check_code(result);
		return result;
	}

	template< typename CharIt >
	void
	chars_for_code( code_t code, CharIt output ) const
	{
		check_code(code);
		for( unsigned i = 0; _alphabet_size != i; ++i )
		{
			if( code & 1 )
			{
				*output = i;
				++output;
			}
			code >>= 1;
		}
	}

	/** output all the codes that have this char in them */
	template< typename CodeIt >
	void
	codes_for_char( char_t c, CodeIt output ) const
	{
		check_char(c);
		codes_for_char_recurse( c, 0, 0, output );
	}

	unsigned
	num_chars_in_code( code_t code ) const
	{
		check_code(code);
		unsigned result = 0;
		for( unsigned i = 0; _alphabet_size != i; ++i )
		{
			if( code & 1 )
				++result;
			code >>= 1;
		}
		return result;
	}

	template< std::size_t L, typename SeqIt >
	void
	count_consensus_L_words( boost::multi_array<unsigned, L> & counts, SeqIt seq_begin, SeqIt seq_end ) const
	{
		//resize and zero the array
		counts.resize( std::vector<unsigned>(L, _num_codes) );
		std::fill( counts.data(), counts.data()+counts.num_elements(), 0 );

		//define iterators for the beginning and end of words
		SeqIt word_begin = seq_begin;
		SeqIt word_end = seq_begin;
		for( unsigned i = 0; seq_end != word_end && L != i; ++i ) ++word_end;
		//did we manage to find one L-word?
		if( word_end - seq_begin < int( L ) ) return; //no

		do
		{
			update_word_counts( counts, word_begin, word_end );
			++word_begin;
			++word_end;
		}
		while( word_end == seq_end );
	};

	template< std::size_t L, typename SeqIt >
	void
	update_word_counts( boost::multi_array<unsigned, L> & counts, SeqIt word_begin, SeqIt word_end ) const
	{
		update_word_counts_recurse<L, boost::multi_array<unsigned, L>, SeqIt>()( *this, counts, word_begin );
	}

private:
	//output all the codes that have this char in them - private recursion func
	template< typename CodeIt >
	void
	codes_for_char_recurse( char_t c, char_t recurse_char, code_t code, CodeIt output ) const
	{
		if( _alphabet_size == recurse_char )
		{
			*output = code;
			++output;
		}
		else
		{
			codes_for_char_recurse( c, recurse_char + 1, code + code_for_char(c), output );
			if( recurse_char != c ) codes_for_char_recurse( c, recurse_char + 1, code, output );
		}
	}

	template< unsigned L, typename CountArray, typename SeqIt >
	struct update_word_counts_recurse
	{
		void operator()( const codes<char_t, code_t> & codes, CountArray & counts, SeqIt word_begin )
		{
			const char_t c = *word_begin;
			++word_begin;
			for( code_t code = 0; codes._num_codes != code; ++code )
				if( codes.matches(c, code) ) {
					typedef typename CountArray::reference count_array_ref;
					count_array_ref ref = counts[code];
					update_word_counts_recurse<L-1, count_array_ref, SeqIt>()( codes, ref, word_begin );
				}
		}
	};
};

template< typename Char, typename Code >
template< typename CountArray, typename SeqIt >
struct codes< Char, Code >::update_word_counts_recurse< 1, CountArray, SeqIt >
{
	void operator()( const codes<Char, Code> & codes, CountArray & counts, SeqIt word_begin )
	{
		typedef Char char_t;
		typedef Code code_t;
		const char_t c = *word_begin;
		for( code_t code = 0; codes._num_codes != code; ++code )
			if( codes.matches(c, code) )
				++counts[code];
	}
};

} //namespace consensus
} //namespace hmm

#endif //HMM_CONSENSUS_H_

