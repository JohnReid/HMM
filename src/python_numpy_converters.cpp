/**
@file

Copyright John Reid 2007, 2012, 2015

*/

#ifdef _MSC_VER
#pragma warning( disable : 4172 )
#endif //_MSVC




#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"

#include "hmm/python.h"

#include <myrrh/python/multi_array_to_numpy.h>

namespace myrrh {
namespace python {

numpy_converter converter;
std::string exposed_typechars;

} //namespace myrrh
} //namespace myrrh


namespace hmm {

HMM_EXPORT_FN_SPEC
void
export_numpy_converters()
{
    using namespace boost::python;
    using namespace indexing;
    using namespace myrrh::python;

    std::cout << "Exposing converters\n";
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
    // expose multiarray to numpy functions
    expose_man_fns();
}

} // namespace hmm
