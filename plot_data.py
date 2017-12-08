# a python 2.7.13 script to plot batch data
#
#   usage: >> python plot_data.py scalar_data mu dtheta=0.000 AR=1.000 bounded=0
#
# datatypes:
#   'mass_density'
#   'orientation_density'
#   'contact_density'
#   dont use this: 'scalar_data' but use one of the two below:
#       'packing_fraction'
#       'stable_fraction'
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

    # deal with scalar data
    if datatype == 'packing_fraction':
        scalartype = 'packing_fraction'
        datatype = 'scalar_data'
    elif datatype == 'stable_fraction':
        scalartype = 'stable_fraction'
        datatype = 'scalar_data'

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

    # Sort the field values & filenames by field values
    vals, files = (list(t) for t in zip(*sorted(zip(vals, files))))

    print(field + ' values: ')
    print(vals)
    print(' ')

    print('filenames: ')
    print(files)
    print(' ')

    data_array = [np.genfromtxt(file, delimiter=',')
                  for file in files]


    # will need to handle data differently for scalar_data
    if datatype == 'scalar_data':
        # get packing fraction
        if scalartype == 'packing_fraction':
            data_array = [d[0] for d in data_array]

        # get stable fraction
        elif scalartype == 'stable_fraction':
            data_array = [d[1] for d in data_array]

        # plot the data
        plt.plot(vals, data_array)
        plt.xlabel(field)
        plt.ylabel(scalartype)

        plt.title(scalartype + ' as a function of ' + field)
        plt.grid(True)
        plt.show()


    # handle vector data
    else:
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










