# a python 2.7.13 script to plot batch data
#
#   usage: >> python plot_data.py scalar_data mu dtheta=0.000 AR=1.000 bounded=0
#
# datatypes:
#   'mass_density'
#   'orientation_density'
#   'contact_density'
#   'scalar_data'
#
# fields:
#   'mu'
#   'dtheta'
#   'AR'
#   'bounded'


import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


# return first cmdline argument as filename
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('datatype') # which output data to look at (e.g. mass density)
    parser.add_argument('field')    # which field to look at (e.g. friction)
    parser.add_argument('ofield0')  # other 3 field values
    parser.add_argument('ofield1')
    parser.add_argument('ofield2')
    args = parser.parse_args()
    return args


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    # import file as numpy record array
    input = get_args()
    field = input.field # ['field']
    datatype = input.datatype # ['datatype']

    # substrings to check all other fields
    otherfields = [input.ofield0.replace('=', '_'),
                   input.ofield1.replace('=', '_'),
                   input.ofield2.replace('=', '_')]


    # where the data files live
    filepath = './outputdata/'
    filenames = os.listdir(filepath)
    vals = []
    files = []

    for filename in filenames:
        # pick out files matching the datatype and other fields
        if datatype in filename and otherfields[0] in filename and otherfields[1] in filename and otherfields[2] in filename:
            # print(filename)

            # get the value of the field
            flist = filename.split('_')
            for (ind, val) in enumerate(flist):
                if val == field:
                    if field == 'bounded':
                        fieldval = int(flist[ind + 1][0])
                    else:
                        fieldval = float(flist[ind + 1])

                    vals.append(fieldval)
                    files.append(filepath + filename)


    print(field + ' values: ')
    print(vals)
    print(' ')

    print('filenames: ')
    print(files)
    print(' ')

    data_array = [np.genfromtxt(file, delimiter=',')
                  for file in files]

    # plot the data
    for (ind, data) in enumerate(data_array):
        data_label = field + ' = ' + str(vals[ind])
        plt.plot(data, label=data_label)

    plt.xlabel('x_axis')
    plt.ylabel('y_axis')

    plt.title(datatype + ' as a function of ' + field)
    plt.grid(True)
    plt.legend()

    plt.show()








#    data = np.genfromtxt(filename, delimiter=' ', names=False)

    # ordered set of field names
    # print( data.dtype.names)

    # vector of data from first field
    # print data['x']
    # print data[data.dtype.names[0]]     # same thing

    # single vector of data
#    if len(data.dtype) == 1:
#        plt.plot(data[data.dtype.names[0]])
#        plt.xlabel('#')
#        plt.ylabel(data.dtype.names[0])
#
#    # plot first two fields
#    else:
#        plt.plot(data[data.dtype.names[0]], data[data.dtype.names[1]])
#        plt.xlabel(data.dtype.names[0])
#        plt.ylabel(data.dtype.names[1])
#
#    plt.title('title')
#    plt.grid(True)
#    plt.savefig(filename[:-4] + '.png')
#    plt.show()














