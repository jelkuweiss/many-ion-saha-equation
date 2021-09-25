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
parser.add_argument('-f', '--fileName', type=str, metavar='', help='State the name of the file to solve the results '
                                                                   'into')
args = parser.parse_args()

if args.densities:
    solution = saha_solver(args.elements, args.densities, args.temperature)
    print("Temperature:", solution[0], 'eV')
    print("Total Number Densities:", solution[1:len(args.elements)+1], 'eV3')
    print("Number Densities of Ionised States:", solution[len(args.elements)+1:len(solution)-2], 'eV3')
    print("Maximum Electron Number Density:", solution[len(solution)-2], 'eV3')
    print("Actual Electron Number Density:", solution[len(solution)-1], 'eV3')
    print("The total ionization in the system is:", round((solution[len(solution)-1] * 100) / solution[len(solution)-2], 1), '%')
    print("--- %s seconds ---" % (time.time() - start_time))
elif args.solarModel:
    input_data = np.loadtxt(args.solarModel, dtype=np.float64)
    rows, columns = input_data.shape
    output_data = open(args.fileName, 'w')
    for i in range(rows):
        T = input_data[i, args.positions[0]]
        dens = []
        for j in args.positions[1:]:
            dens += [input_data[i, j]]
        solution = saha_solver(args.elements, dens, T)
        print(*solution, sep=", ", file=output_data)
    print("--- Executed in %s seconds ---" % (time.time() - start_time))
