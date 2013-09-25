#include <time.h>
#include <iostream>
#include "boost/multi_array.hpp"

using namespace std;
using namespace boost;

int main()
{
    const size_t DimX = 200;
    const size_t DimY = 150;
    const size_t DimZ = 200;
    multi_array<int,3> V3(extents[DimZ][DimY][DimX]);
    const unsigned num_iters = 1000;

    for (size_t Z=0; Z!=DimZ; ++Z) {
        for (size_t Y=0; Y!=DimY; ++Y) {
            for (size_t X=0; X!=DimX; ++X) {
                V3[Z][Y][X] = Z + Y + X;
            }
        }
    }

    int iVal = 0;
    //unsigned long time1 = GetTickCount();
    //time_t time1 = time(0);
    clock_t time1 = clock();
    for( unsigned i = 0; num_iters != i; ++i ) {
        for (size_t Z=0; Z!=DimZ; ++Z) {
            for (size_t Y=0; Y!=DimY; ++Y) {
                for (size_t X=0; X!=DimX; ++X) {
                    iVal += V3[Z][Y][X];
                }
            }
        }
    }
    //unsigned long time2 = GetTickCount();
    //time_t time2 = time(0);
    clock_t time2 = clock();
    //cout << "Multi_array: " << (time2-time1) << endl;
    cout << "Multi_array: " << double(time2-time1)/CLOCKS_PER_SEC << endl;

    int * pInt = V3.data();
    //time1 = GetTickCount();
    //time1 = time(0);
    time1 = clock();
    for( unsigned i = 0; num_iters != i; ++i ) {
        for (size_t Z=0; Z!=DimZ; ++Z) {
            for (size_t Y=0; Y!=DimY; ++Y) {
                for (size_t X=0; X!=DimX; ++X) {
                    iVal += pInt[(Z*DimY+Y)*DimX+X];
                }
            }
        }
    }
    //time2 = GetTickCount();
    //time2 = time(0);
    time2 = clock();
    //cout << "Pointer....: " << (time2-time1) << endl;
    cout << "Pointer....: " << double(time2-time1)/CLOCKS_PER_SEC << endl;

    cout << "Sum of values: " << iVal << endl;

    return 0;
}

