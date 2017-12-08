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

    float friction[9]   = { 0, .05, .1, .15,.2,.25,.3,.35,.4 };
    float dropangle[1]  = { //0.,
                            1. * M_PI / 8.//,
                            //2. * M_PI / 8.,
                            //3. * M_PI / 8.,
                            //4. * M_PI / 8.,
                          };
    float AR[1]         = { 55 };
    int boundary[1]     = { 1 };


    for (int f = 0; f < 9; f++) {
        for (int d = 0; d < 1; d++) {
            for (int a = 0; a < 1; a++) {
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







