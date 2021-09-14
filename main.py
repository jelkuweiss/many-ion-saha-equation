from functions import *
import numpy as np
import argparse
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Solve the Saha equation for a list of Elements')
group = parser.add_mutually_exclusive_group()

parser.add_argument('-e', '--elements', nargs='+', type=str, metavar='', help='List of element symbols')

group.add_argument('-d', '--densities', nargs='+', type=float, metavar='', help='List of element densities')
parser.add_argument('-t', '--temperature', type=float, metavar='', help='Temperature')

group.add_argument('-s', '--solarModel', type=str, metavar='', help='File name of the Solar Model')
parser.add_argument('-p', '--positions', nargs='+', type=int, metavar='', help='Position of Temperature and '
                                                                               'then elements in Solar Model ')
args = parser.parse_args()

if args.densities:
    M, max_ne = augmented_matrix_builder(args.elements, args.densities, args.temperature)
    print("Matrix Dimensions: %s" % len(M))
    s, actual_ne = matrix_solver_for_ne(M, max_ne)
    percentile_io = round((actual_ne*100)/max_ne, 1)
    print("The total ionization in the system is:", percentile_io, '%')
    print(s)
    print("--- %s seconds ---" % (time.time() - start_time))
elif args.solarModel:
    input_data = np.loadtxt(args.solarModel, dtype=np.float64)
    rows, columns = input_data.shape
    output_data = open('results.txt', 'a')
    for i in range(rows):
        T = input_data[i, args.positions[0]]
        dens = []
        for j in args.positions[1:]:
            dens += [input_data[i, j]]
        M, max_ne = augmented_matrix_builder(args.elements, dens, T)
        S, actual_ne = matrix_solver_for_ne(M, max_ne)
        input_data[i, args.positions].tofile(output_data, sep=' ', format='%s')
        output_data.write(" ")
        S.tofile(output_data, sep=' ', format='%s')
        output_data.write("\n")
    print("--- Executed in %s seconds ---" % (time.time() - start_time))
