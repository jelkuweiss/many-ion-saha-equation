from functions import *
import numpy as np


filename = "solar_density_data.txt"
input_data = np.loadtxt(filename, dtype=np.float64)
rows, columns = input_data.shape
output_data = open('results_new.txt', 'a')

for i in range(rows):
    nH = (input_data[i, 6])*4.31013646*(10**18)/(0.93878405573*10**9)
    nHe = (input_data[i, 7])*4.31013646*(10**18)/(0.372784341*10**10)
    T = (input_data[i, 2])/11604.51812
    M = augmented_matrix_builder(["H", "He"], [nH, nHe], T)
    s = matrix_solver_for_ne(M,10**12)
    input_data[i, 0:8].tofile(output_data, sep=' ', format='%s')
    output_data.write(" ")
    s.tofile(output_data, sep=' ', format='%s')
    output_data.write("\n")
    print(i/rows * 100)

#M = saha.augmented_matrix_builder(["H", "He"], [1662424176, 719087680], 1.549e+07)
#s = saha.matrix_solver_for_ne(M, 3100599537)