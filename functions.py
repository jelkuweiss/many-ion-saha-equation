# Importing
import numpy as np
from sympy import *
import math
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# List of dictionaries for element properties
# First dictionary contains the statistical weight of each ionised state of the atom
# Second one is the energy difference xi between the states
# Third one is the number of ionised states of each element
# Fourth one is the mass in natural units of the elements

el_prop = [
    {
        "H": [2, 1],
        "He": [1, 2, 1],
        "He3": [1, 2, 1],
        "C": [1, 1, 1, 1, 1, 1, 1],  # Fix@@
        "N": [1, 1, 1, 1, 1, 1, 1, 1],
        "O": [1, 1, 1, 1, 1, 1, 1, 1, 1],
        "Ne": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "Si": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        "Fe": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    },
    {
        "H": [13.59],
        "He": [24.58, 54.4],
        "He3": [24.58, 54.4],
        "C": [11.26, 24.38, 47.89, 64.49, 392.09, 489.99],
        "N": [14.53, 29.60, 47.44, 77.47, 97.89, 552.07, 667.05],
        "O": [13.62, 35.12, 54.94, 77.41, 113.90, 138.12, 739.33, 871.41],
        "Ne": [21.56, 40.96, 63.42, 97.19, 126.25, 157.93, 207.27, 239.10, 1195.81, 1362.20],
        "Si": [8.15, 16.35, 33.49, 45.14, 166.77, 205.28, 246.57, 303.59, 351.28, 401.38, 476.27, 523.42, 2437.66,
               2673.18],
        "Fe": [7.90, 16.20, 30.65, 54.91, 75.00, 98.99, 124.98, 151.06, 233.6, 262.10, 290.9, 330.8, 361.0, 392.2,
               456.2, 489.31, 1262.7, 1357.8, 1460, 1575.6, 1687.0, 1798.4, 1950.4, 2045.759, 8828.1879, 9277.6818],
    },
    {
        "H": 1,
        "He": 2,
        "He3": 2,
        "C": 6,
        "N": 7,
        "O": 8,
        "Ne": 10,
        "Si": 14,
        "Fe": 26,
    },
    {
        "H": 940190710.547166,
        "He": 3733991792.654492,
        "He3": 3.01602932265 * 9.328908102 * (10 ** 8),
        "C": 12 * 9.328908102 * (10 ** 8),
        "N": 14.00307400446 * 9.328908102 * (10 ** 8),
        "O": 15.99491461960 * 9.328908102 * (10 ** 8),
        "Ne": 19.9924401762 * 9.328908102 * (10 ** 8),
        "Si": 27.9769265350 * 9.328908102 * (10 ** 8),
        "Fe": 55.9349363 * 9.328908102 * (10 ** 8),
    },
]

# Global Constants
h = 2 * math.pi
kb = 1
me = 511e3  # [eV]


# This function represents the right hand side of the saha equation
# It will in the matrix builder to fill in for the values
# T is the temperature, xi is the energy difference of the states
# z0 and z1 are the statistical weights of the states (lowest and highest, respectively)
# (lam means lamda, the usual symbol associated with this side of the equation)
def lam(T, xi, z0, z1):
    return (z1 / z0) * ((2 * math.pi * me * kb * T) ** 1.5 / (h ** 3)) * math.exp((-xi) / (kb * T))


# This is the main function. It takes a temperature [K] and mass density [g/cm3] and transforms them to eV,
# then solves the Saha Equation.
# Output is of the form [T(eV), nH(eV3), nHe(eV3), ..., nH0(eV3), nH1(eV3), nHe0(eV3), ..., nE(eV3)]
def saha_solver(el_names, el_dens, t):
    # We need the total ionisation levels from all our elements to know the dimension of the matrix
    total_io_levels = 0
    for el in el_names:
        total_io_levels += el_prop[2][el]
    dim_h = 1 + len(el_names) + total_io_levels  # Horizontal dimension of M
    nE = symbols('n_e')

    # Before starting we need to transform the masses given in g/cm3 into number densities in eV^4 then eV^3 by
    # diving with the mass of the element
    for ii, ell in enumerate(el_dens):
        # Take the mass, transform them into eV^4, divide by the atomic mass to get a number density in eV^3
        el_dens[ii] = ell * 4.29553 * (10 ** 18) / el_prop[3][el_names[ii]]

    # Also convert the Temperature to eV
    t = t / 11604.51812

    # maximum electron number density which we get in the below loop by assuming all elements lost their electrons (
    # total ionisation). So just by multiplying the number density of each element with its maximum number of
    # electrons and summing for all elements
    max_ne = 0
    # Initialising the Matrix [M] list. We build it first as a list of lists and later on transform it to Sympy Matrix
    M = []
    row_buffer = []  # used for the loop only
    skip_counter_1 = 0  # used for the loop only
    skip_counter_2 = 0  # used for the loop only

    # The loop here is building the M lists
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
        max_ne += el_dens[i] * el_prop[2][el]

    # Electron Density row
    row_buffer = []
    for el in el_names:
        for k in range(el_prop[2][el] + 1):
            row_buffer += [k]
    row_buffer += [nE]
    M.append(row_buffer)

    # Transform M into a Sympy matrix
    M = Matrix(M)

    # The root scalar will go through the function trying every electron density to solve the augmented matrix
    def determinant_polynomial(ne):
        m = M.subs(nE, ne)
        return m.det()

    sol = root_scalar(determinant_polynomial, method='bisect', bracket=(0, max_ne))

    # Substitute ne for the root we found
    M.subs(nE, sol.root)

    # Convert it now to a numpy matrix and take the system part(not augmented) out
    MM = np.array(np.array(M.subs(nE, sol.root)), np.float64)
    M1 = MM[0:(len(MM[0]) - 1), 0:(len(MM[0]) - 1)]
    M2 = MM[0:(len(MM[0]) - 1), (len(MM[0]) - 1)]

    # The solution in form of (nh0,nh1,nhe0,nhe1,nhe2,...)
    ionised_number_densities = np.linalg.solve(M1, M2)

    # We start building the list we want to return, by first putting in the temperature [eV]
    return_list = [t]
    # Now we append the total number densities [eV3]
    for x in el_dens:
        return_list.append(x)
    # Now append the ionisation number densities [eV3]
    for x in ionised_number_densities:
        return_list.append(x)
    # Now append the electron number density [eV3] and the maximal one too
    return_list.append(max_ne)
    return_list.append(sol.root)

    return return_list
