// batch script


// add to makefile:
// batcher:
//      g++ -Wall -o batcher batcher.cpp


// usage: ./nest 10000     1.57        50  1
//               friction, drop angle, AR, boundary

#include <cstdlib>  // system
#include <iomanip>  // setprecision
#include <iostream>
#include <sstream>  // stringstream
#include <cmath>    // M_PI


const char* inputStr(float f, float d, float a, int b) {
    std::stringstream stream;

    stream  << "./nest ";

    stream  << std::fixed
            << std::setprecision(5) // five total characters. e.g. pi --> "3.14"
            << f << " "
            << d << " "
            << a << " ";

    stream  << std::fixed
            << std::setprecision(1)
            << b;

//    std::cout << stream.str().c_str() << std::endl;

    return stream.str().c_str();
}


int main(void) {

    float friction[1]   = { atan(0.) };//, atan(0.25), atan(0.5), atan(0.75), atan(1) };
    float dropangle[1]  = { 0. };//,
//                            1. * M_PI / 8.,
//                            2. * M_PI / 8.,
//                            3. * M_PI / 8.,
//                            4. * M_PI / 8.,
//                          };
    float AR[1]         = { 0.01 };//, 1., 3., 8., 25., 55. };
    int boundary[1]     = { 0 };//, 1 };


    for (int f = 0; f < 1; f++) {
        for (int d = 0; d < 1; d++) {
            for (int a = 0; a < 1; a++) {
                for (int b = 0; b < 1; b++) {

//                    printf("(friction, dropangle, AR, boundary) = (%f, %f, %f, %d)\n",
//                           friction[f], dropangle[d], AR[a], boundary[b]);

                    const char* call = inputStr(friction[f], dropangle[d],
                                                AR[a], boundary[b]);

//                    printf("%c\n", call[0]);

                    system(call);

                }
            }
        }
    }

    return 0;
}







