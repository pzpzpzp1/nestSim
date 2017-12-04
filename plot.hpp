#ifndef PLOT_HPP
#define PLOT_HPP

// A little terminal plotter to stdout



// START

//    std::vector< std::vector<int> > data;
//    std::vector<int> col_0;
//        col_0.push_back(1);
//        col_0.push_back(1);
//        col_0.push_back(0);
//        col_0.push_back(1);
//    std::vector<int> col_1;
//        col_1.push_back(0);
//        col_1.push_back(1);
//        col_1.push_back(1);
//        col_1.push_back(1);
//    data.push_back(col_0);
//    data.push_back(col_1);
//
//    std::vector<int> x_axis;
//    x_axis.push_back(1);
//    x_axis.push_back(2);
//    x_axis.push_back(3);
//    x_axis.push_back(14);
//
//    std::vector<int> y_axis;
//    y_axis.push_back(1);
//    y_axis.push_back(200);
//
//    _plot(data, x_axis, y_axis);
//
//
//    return 0;




// END

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <assert.h>

template <class T>
void plot(std::vector<T>, T, bool);
template <class T>
void _plot(std::vector< std::vector<T> >, std::vector<int>, std::vector<int>);
int _width(int);

// xy offset
// xy range
// ...

// Plot a single vector of data assuming x = [1, 2, 3, ...];
template <class T>
void plot(std::vector<T> y, T y_step, bool fill) {
    int n = y.size();
    int y_max = ceil(*max_element(y.begin(), y.end()));
    int y_min = floor(*min_element(y.begin(), y.end()));



    return;
}

// Draw the actual plot from a vector of row data
template <class T>
void _plot(std::vector< std::vector<T> > data,
           std::vector<int> x_axis,
           std::vector<int> y_axis) {

    int y_axis_width = _width(y_axis[y_axis.size() - 1]);

    assert(x_axis.size() == data[0].size());
    assert(y_axis.size() == data.size());

    printf("\n");

    for (int row = data.size() - 1; row >= 0; row--) {
        // Print y-axis heading
        printf("%*d| ", y_axis_width, y_axis[row]);

        for (int col = 0; col < data[0].size(); col++) {
            if (data[row][col]) {
                printf("o ");
            }
            else {
                printf("  ");
            }
        }

        printf("\n");
    }

    // Print x-axis
    int n_spaces = y_axis_width + 1;
    int n_dashes = data[0].size() * 2;
    for (int i = 0; i < n_spaces; i++) { printf(" "); }
    for (int i = 0; i < n_dashes; i++) { printf("-"); }
//    printf("\n");

    // Print x-axis values
//    for (int i = 0; i < n_spaces; i++) { printf(" "); }
//    for (int col = 0; col < data[0].size(); col++) {
//        printf("%2d", x_axis[col]);
//    }

    printf("\n\n");

    return;
}

// Return the number of spaces needed to print an integer
int _width(int number) {
    return floor(log10(number)) + 1;
}

#endif // PLOT_HPP



























