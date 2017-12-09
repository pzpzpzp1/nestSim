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
#include <stdio.h>

using namespace std;


int main(void) {

    char buffer[1000];

    float friction[1]   = { .64 };
    float dropangle[9]  = { 0.,
                            1. * M_PI / 16.,
                            2. * M_PI / 16.,
                            3. * M_PI / 16.,
                            4. * M_PI / 16.,
                            5. * M_PI / 16.,
                            6. * M_PI / 16.,
                            7. * M_PI / 16.,
                            8. * M_PI / 16.,
                          };
    float AR[3]         = { 8,25,55 };
    int boundary[1]     = { 0 };


    for (int f = 0; f < 1; f++) {
        for (int d = 0; d < 9; d++) {
            for (int a = 0; a < 6; a++) {
                for (int b = 0; b < 1; b++) {

                    printf("(friction, dropangle, AR, boundary) = (%f, %f, %f, %d)\n",
                           friction[f], dropangle[d], AR[a], boundary[b]);

                    sprintf(buffer,"./nest %f %f %f %d\n",friction[f], dropangle[d], AR[a], boundary[b]);

                    system(buffer);

                }
            }
        }
    }

    return 0;
}







