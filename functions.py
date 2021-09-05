#Importing
import numpy as np
from sympy import *
import math
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# List of dictionaries for element properties
# First dictionary contains the statistical weight of each ionised state of the atom
# Second one is the energy difference xi between the states
# Third one is the number of ionised states of each element

el_prop = [
    {
        "H": [2, 1],
        "He": [1, 2, 1],
        "C": [1, 1, 1, 1, 1, 1, 1],#Fix@@
        "N": [1, 1, 1, 1, 1, 1, 1, 1],
        "O": [1, 1, 1, 1, 1, 1, 1, 1, 1],
        "Fe": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    },
    {
        "H": [13.59],
        "He": [24.58, 54.4],
        "C": [11.26, 24.38, 47.89, 64.49, 392.09, 489.99],
        "N": [14.53, 29.60, 47.44, 77.47, 97.89, 552.07, 667.05],
        "O": [13.62, 35.12, 54.94, 77.41, 113.90, 138.12, 739.33, 871.41],
        "Fe": [7.90, 16.20, 30.65, 54.91, 75.00, 98.99, 124.98, 151.06, 233.6, 262.10, 290.9, 330.8, 361.0, 392.2, 456.2, 489.31, 1262.7, 1357.8, 1460, 1575.6, 1687.0, 1798.4, 1950.4, 2045.759, 8828.1879, 9277.6818],
    },
    {
        "H": 1,
        "He": 2,
        "C": 6,
        "N": 7,
        "O": 8,
        "Fe": 26,
    }
]


# Global Constants
h = 2*math.pi
kb = 1
me = 511e3 #[eV]


# This function represents the right hand side of the saha equation
# It will in the matrix builder to fill in for the values
# T is the temperature, xi is the energy difference of the states
# z0 and z1 are the statistical weights of the states (lowest and highest, respectively)
# (lam means lamda, the usual symbol associated with this side of the equation)
def lam(T, xi, z0, z1):
    return (z1/z0) * ((2*math.pi*me*kb*T)**1.5 / (h**3)) * math.exp((-xi) / (kb*T))


# This function builds the matrix that represents the system of many equations
# It takes in the list of elements, as well as their respective densities and the temperature in the region
# el_names is a list of names ['H','He']
# el_dens is a list of densities [5000,6000]
# T is the temperature
def augmented_matrix_builder(el_names, el_dens, t):
    # We need the total ionisation levels from all our elements to know the dimension of the matrix
    total_io_levels = 0
    for el in el_names:
        total_io_levels += el_prop[2][el]
    dim_h = 1 + len(el_names) + total_io_levels  # Horizontal dimension of M
    global nE
    nE = symbols('n_e')

    max_ne = 0
    M = []
    row_buffer = []
    skip_counter_1 = 0
    skip_counter_2 = 0

    for i, el in enumerate(el_names):
        for j in range(el_prop[2][el]):
            # Insert the saha equation line for each energy transition
            row_buffer = (skip_counter_1) * [0] + [-lam(t, el_prop[1][el][j], el_prop[0][el][j], el_prop[0][el][j + 1]),
                                                   nE] + (dim_h - skip_counter_1 - 2) * [0]
            M.append(row_buffer)
            # Skip a column next time
            skip_counter_1 += 1
        # Skip extra column after each element
        skip_counter_1 += 1

        # Total density row
        row_buffer = (skip_counter_2) * [0] + (el_prop[2][el] + 1) * [1] + (
                    dim_h - 1 - skip_counter_2 - el_prop[2][el] - 1) * [0] + [el_dens[i]]
        M.append(row_buffer)
        skip_counter_2 += el_prop[2][el] + 1

        # We also calculate the maximum energy density since we are looping
        max_ne += el_dens[i]*el_prop[2][el]

    # Electron Density row
    row_buffer = []
    for el in el_names:
        for k in range(el_prop[2][el] + 1):
            row_buffer += [k]
    row_buffer += [nE]
    M.append(row_buffer)

    # Should we return a list of lists (as we do here) or a sympy matrix
    return M, max_ne

# This function solves an augmented matrix containing an electron density symbol(sympy)
def matrix_solver_for_ne(M,max_ne):
    # Transform the list of lists into a sympy matrix
    M = Matrix(M)

    # Plotting
    #ne_values = range(min_ne, max_ne+1, round(max_ne/100))
    #det_values = []
    #for val in ne_values:
    #    det_values.append(determ_augm.subs(nE, val))
    #plt.scatter(ne_values, det_values)

    # The root scalar will go through the function trying every electron density to solve the augmented matrix
    def determinant_polynomial(ne):
        m = M.subs(nE, ne)
        return m.det()
    sol = root_scalar(determinant_polynomial, method='bisect', bracket=(0, max_ne))

    # Substitute ne for the root we found
    M.subs(nE, sol.root)

    # Convert it now to a numpy matrix and take the system part(not augmented) out
    MM = np.array(np.array(M.subs(nE, sol.root)), np.float64)
    M1 = MM[0:(len(MM[0])-1), 0:(len(MM[0])-1)]
    M2 = MM[0:(len(MM[0])-1),   (len(MM[0])-1)]

    # The solution in form of (nh0,nh1,nhe0,nhe1,nhe2)
    xx = np.linalg.solve(M1, M2)

    return xx
