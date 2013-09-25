/**
@file

Copyright John Reid 2009

*/

#include <boost/multi_array.hpp>

int
main( int argc, char * argv[] ) {
    boost::multi_array< double, 2 > a;
    a.resize( boost::extents[2][3] );
    return 0;
}
