/**
@file

Copyright John Reid 2007, 2008, 2009, 2012, 2015

*/

#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"

#include <boost/python/detail/prefix.hpp>
#include "hmm/python.h"

using namespace boost;
using namespace std;

using namespace hmm;

int major_version() { return 2; }
int minor_version() { return 0; }
std::string version() { return HMM_MAKE_STRING( major_version() << "." << minor_version() ); }


void test_vec( const std::vector< unsigned int > & py_obj ) {
}

void test_ma_int1( const boost::multi_array_ref< int, 1 > & py_obj ) {
}

void test_ma_double1( const boost::multi_array_ref< double, 1 > & py_obj ) {
}

void test_ma_uint1( const boost::multi_array_ref< unsigned int, 1 > & py_obj ) {
}

void test_ma_uint2( const boost::multi_array_ref< unsigned int, 2 > & py_obj ) {
}


BOOST_PYTHON_MODULE( _hmm )
{
    using namespace boost::python;
    using boost::python::arg;

    // for debugging
    def( "test_ma_double1", test_ma_double1 );
    def( "test_ma_int1", test_ma_int1 );
    def( "test_ma_uint1", test_ma_uint1 );
    def( "test_ma_uint2", test_ma_uint2 );
    def( "test_vec", test_vec );

    scope().attr("__doc__") =
        "Python interface to C++ hidden Markov model library.\r\n"
        "Implements HMMs of higher Markov orders.";

    scope().attr("__major_version__") = major_version();
    scope().attr("__minor_version__") = minor_version();
    scope().attr("__version__") = version();

#ifndef NDEBUG
    scope().attr("__debug__") = 1;
    std::cout << "WARNING: Debug version of _hmm module loaded. "
        "If you did not intend this then check your configuration!" << std::endl;
#else //_DEBUG
    scope().attr("__debug__") = 0;
#endif //_DEBUG

    export_numpy_converters();
    export_model();
    export_markov_order();
    export_model_by_states();
    export_simple_pssm();
}


