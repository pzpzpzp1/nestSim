// batch script


// add to makefile:
// batcher:
//      g++ -Wall -o batcher batcher.cpp


// usage: ./nest 10000     1.57        50  1
//               friction, drop angle, AR, boundary

#include <cstdlib>  // system
#include <iomanip>  // setprecision
#include <sstream>  // stringstream
#include <cmath>    // M_PI


const char* inputStr(float f, float d, float a, int b) {
    std::stringstream stream;

    stream  << "./nest";

    stream  << std::fixed
            << std::setprecision(5)
            << f << " "
            << d << " "
            << a << " ";

    stream  << std::fixed
            << std::setprecision(1)
            << b;

    return stream.str().c_str();
}


int main(void) {

    float friction[5]   = { atan(0.), atan(0.25), atan(0.5), atan(0.75), atan(1) };
    float dropangle[5]  = { 0.,
                            1. * M_PI / 8.,
                            2. * M_PI / 8.,
                            3. * M_PI / 8.,
                            4. * M_PI / 8.,
                          };
    float AR[6]         = { 0.01, 1., 3., 8., 25., 55. };
    int boundary[2]     = { 0, 1 };


    for (int f = 0; f < 5; f++) {
        for (int d = 0; d < 5; d++) {
            for (int a = 0; a < 6; a++) {
                for (int b = 0; b < 2; b++) {

                    printf("(friction, dropangle, AR, boundary) = (%f, %f, %f, %d)\n",
                           friction[f], dropangle[d], AR[a], boundary[b]);

                    const char* call = inputStr(friction[f], dropangle[d],
                                                AR[a], boundary[b]);
                    system(call);

                }
            }
        }
    }

    return 0;
}







