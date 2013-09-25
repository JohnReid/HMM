/**
@file

Copyright John Reid 2007, 2008, 2012

*/

#ifndef HMM_DEFS_H_
#define HMM_DEFS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER


#include <myrrh/assert.h>

#ifdef NDEBUG

# define HMM_VERIFY( x )

#else //NDEBUG

# define HMM_VERIFY_ON
# define BOOST_ENABLE_ASSERT_HANDLER
# define MYRRH_ASSERT_THROWS
# define HMM_VERIFY( x ) MYRRH_ASSERT( x )

namespace boost {
inline
void
assertion_failed( char const * expr, char const * function, char const * file, long line )
{
	myrrh::assertion_throw( expr, function, file, line );
}
}

#endif //NDEBUG




#include <exception>





#include <myrrh/types.h>
#include <myrrh/trie.h>


#include <boost/config.hpp>





//need to be included first
#ifndef DONT_USE_BOOST_SERIALIZATION
#ifdef _MSC_VER
# pragma warning( push )
# pragma warning( disable : 4099 )
# pragma warning( disable : 4996 )
#endif //_MSC_VER
#  include <boost/archive/text_oarchive.hpp>
#  include <boost/archive/text_iarchive.hpp>
#  include <boost/archive/binary_oarchive.hpp>
#  include <boost/archive/binary_iarchive.hpp>
#  include <boost/serialization/map.hpp>
#  include <boost/serialization/set.hpp>
#  include <boost/serialization/shared_ptr.hpp>
#  include <boost/serialization/vector.hpp>
#  include <boost/serialization/list.hpp>
#ifdef _MSC_VER
# pragma warning( pop )
#endif //_MSC_VER
#endif



#include <boost/test/utils/wrap_stringstream.hpp>
#define HMM_MAKE_STRING(x) ( boost::wrap_stringstream().ref() << x ).str()


#include <boost/iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/operators.hpp>
#include <boost/array.hpp>
#include <boost/regex.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multi_array.hpp>
#include <boost/optional.hpp>

#include <string>
#include <vector>
#include <numeric>
#include <exception>

#include <float.h>



namespace hmm {

typedef unsigned output_character;
typedef std::vector< output_character > output_seq;


} //namespace hmm




#ifdef _MSC_VER
# define HMM_ISNAN( x ) _isnan( x )
#else
# define HMM_ISNAN( x ) isnan( x )
#endif


#endif //HMM_DEFS_H_

