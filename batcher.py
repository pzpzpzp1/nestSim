# a python 2.7.13 script to plot outputs from nestSim
#
#   usage: >> python plotter.py filename
#
#          where filename is a csv file, like 'example.csv' delimited by spaces
#          since apparently genfromtxt can handle hanging spaces but not commas
#
#               field1 field2 field3 ...
#               val1_1 val2_1 val3_1 ...
#               val1_2 val2_2 val3_2 ...
#               ...

import numpy as np
import matplotlib.pyplot as plt
import argparse


# return first cmdline argument as filename
def get_args():
    parser = argparse.ArgumentParser(description='Plot data from CSV file.')
    parser.add_argument('filename')
    args = parser.parse_args()
    return args.filename


if __name__ == "__main__":
    # import file as numpy record array
    filename = get_args()
    data = np.genfromtxt(filename, delimiter=' ', names=True)

    # ordered set of field names
    print( data.dtype.names)

    # vector of data from first field
    # print data['x']
    # print data[data.dtype.names[0]]     # same thing

    # single vector of data
    if len(data.dtype) == 1:
        plt.plot(data[data.dtype.names[0]])
        plt.xlabel('#')
        plt.ylabel(data.dtype.names[0])

    # plot first two fields
    else:
        plt.plot(data[data.dtype.names[0]], data[data.dtype.names[1]])
        plt.xlabel(data.dtype.names[0])
        plt.ylabel(data.dtype.names[1])

    plt.title('title')
    plt.grid(True)
    plt.savefig(filename[:-4] + '.png')
    plt.show()














