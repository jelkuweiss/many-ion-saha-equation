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
parser.add_argument('-p', '--elementPosition', nargs='+', type=int, metavar='', help='Position of elements in Solar Model file')
args = parser.parse_args()

if args.densities:
    M, max_ne = augmented_matrix_builder(args.elements, args.densities, args.temperature)
    s = matrix_solver_for_ne(M, max_ne)
    print(s)
    print("--- %s seconds ---" % (time.time() - start_time))
elif args.solarModel:
    print('helo')
    input_data = np.loadtxt(args.solarModel, dtype=np.float64)
    rows, columns = input_data.shape
    output_data = open('results.txt', 'a')
    for i in range(rows):
        dens = []
    # ---- UNDER CONSTRUCTION ---
    # for i in range(rows):
    # for i in range(rows):
    #    nH = (input_data[i, 6])*4.31013646*(10**18)/(0.93878405573*10**9)
    #    nHe = (input_data[i, 7])*4.31013646*(10**18)/(0.372784341*10**10)
    # T = (input_data[i, 2])/11604.51812
    # M = augmented_matrix_builder(["H", "He"], [nH, nHe], T)
    # s = matrix_solver_for_ne(M,10**12)
    # input_data[i, 0:8].tofile(output_data, sep=' ', format='%s')
    # output_data.write(" ")
    # s.tofile(output_data, sep=' ', format='%s')
    # output_data.write("\n")
    # print(i/rows * 100)
    # ---- UNDER CONSTRUCTION ---
    print("--- Executed in %s seconds ---" % (time.time() - start_time))

# M = saha.augmented_matrix_builder(["H", "He"], [1662424176, 719087680], 1.549e+07)
# s = saha.matrix_solver_for_ne(M, 3100599537)
