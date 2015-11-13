import csv
import sys

import numpy as np
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt


def read_data_as_dict(file):
    '''
    Builds a dictionary out of the expression CSV.
    The dict has the following format:
    {'CELL_LINE_1':{
            'GENE1':1234, 'GENE2':2345, ...
        }
    'CELL_LINE_2':{
            'GENE1':1234, 'GENE2':2345, ...
        }
        .
        .
        .

    }
    :param file: Path to the csv file.
    :return: Dict (as above format).
    '''
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file)
        cells = next(csv_reader, [])
        cells.pop(0)
        dict = {}
        for cell in cells:
            dict[cell] = {}
        for c in csv_reader:
            for k, l in enumerate(c[1:]):
                dict[cells[k]][c[0]] = float(l)

    return dict


def convert_dict_to_matrix(cell_dict):
    '''
    Clustering algorithms require a numpy matrix.
    Take a dictionary (from the read_data_as_dict) and make it into a numpy matrix.
    :param cell_dict: Dict from read_data_as_dict.
    :return: a numpy matrix s.t. every row is a vector of exp. of genes
    '''
    matrix = []
    for c in cell_dict:
        arr = []
        for d in cell_dict[c]:
            arr.append(cell_dict[c][d])
        matrix.append(arr)
    return np.array(matrix)


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print("Please run as python clusterExpression.py <path_to_genx.csv>")
    # Read the csv as a 2d-dictionary/matrix
    dict_matrix = read_data_as_dict(sys.argv[1])
    matrix = convert_dict_to_matrix(dict_matrix)
    # Perform linkage / clustering using wards method
    # https://en.wikipedia.org/wiki/Ward%27s_method
    z = hac.linkage(matrix, method='ward')
    # Construct a dendrogram from that clustering
    d = hac.dendrogram(z, labels=list(dict_matrix.keys()))
    # Show a plot of the dendrogram.
    plt.show()
