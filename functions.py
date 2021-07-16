#Importing
import numpy as np
from sympy import *
import math

#List of dictionaries for element properties
#first dictionary contains the statistical weight of each ionised state of the atom
#Second one is the energy difference xi between the states
#Third one is the number of ionised states of each element

element_properties = [
    {
        "H": [2, 1],
        "He":[1, 2, 1],
        "C": [1, 1, 1, 1, 1, 1, 1],#Fix
        "N": [1, 1, 1, 1, 1, 1, 1, 1],
        "O": [1, 1, 1, 1, 1, 1, 1, 1, 1],
    },
    {
        "H": [13.59],
        "He":[24.58, 54.4],
        "C": [11.26, 24.38, 47.89, 64.49, 392.09, 489.99],
        "N": [14.53, 29.60, 47.44, 77.47, 97.89, 552.07, 667.05],
        "O": [13.62, 35.12, 54.94, 77.41, 113.90, 138.12, 739.33, 871.41],
    },
    {
        "H": 1,
        "He":2,
        "C": 6,
        "N": 7,
        "O": 8,
    }
]


# Global Constants
h =  2*math.pi
kb = 1
me = 511e3 #[eV]


#This function represents the right hand side of the saha equation
#It will in the matrix builder to fill in for the values
#T is the temperature, xi is the energy difference of the states
#z0 and z1 are the statistical weights of the states (lowest and highest, respectively)
# (lam means lamda, the usual symbol associated with this side of the equation)
def lam(T,xi,z0,z1):
    return (z1/z0) * ((2*math.pi*me*kb*T)**1.5 / (h**3)) * math.exp((-xi) / (kb*T))


#This function builds the matrix that represents the system of many equations
#It takes in the list of elements, as well as their respective densities and the temperature in the region
#el_names is a list of names ['H','He']
#el_dens is a list of densities [5000,6000]
#T is the temperature
def augmented_matrix_builder(el_names,el_dens,T):
    #First we need to get the dimension of the vector in the matrix:
    #Get the total ionisation levels
    total_io_levels = 0
    for element in el_names:
        total_io_levels = total_io_levels + element_properties[2][element]
    #Dimesion of rows of augmented matrix
    dim_v = 1 + len(el_names) + total_io_levels

    #Define the symbol of the energy density
    nE = symbols('n_e')

    #Now we build the Matrix:
    #If you take the example of Hydrogen and Helium (1 + 2 saha equations; 2 total density equations; 1 electron density equation)
    #you will find that you can order them in block matrices (all but the electron density one)
    #We will use this fact to construct the matrix with loops
    M = [] #The matrix
    row = [] #Row buffer to build the matrix
    column_skip_counter = 0 #skip counter to help us build a column matrix (placing zeroz:)

    for e in range(len(el_names)):

        element = el_names[e]
        dim_block = element_properties[2][element] + 1 #size of the block matrix of each element
        block_column_skip_counter = 0 #skip counter for each element block matrix (they are mostly upper triangular)

        for i in range(dim_block - 1): #the -1 is because the last line of the block is just a list of 1's
            column_skip_counter = column_skip_counter + 1

            #first we do all the skipping we have to do to reach the place of our element block
            for k in range(0,column_skip_counter-1):
                row.append(0)

            #Then we build the block lines: first by skipping enough to make the block upper triangular
            for k in range(column_skip_counter,column_skip_counter+block_column_skip_counter-1):
                row.append(0)

            #Then we place the two elements
            row.append(-lam(T, element_properties[1][element][i], element_properties[0][element][i], element_properties[0][element][i+1]))
            row.append(nE)

            #fill the rest with zeros
            for k in range(column_skip_counter+block_column_skip_counter + 2 -1, dim_v-1): #-1 here since the last column is built seperately
                row.append(0)

            #Append to the matrix, clear the row, and do it again
            M.append(row)
            row = []

        #Now the last line in the block is the total density line which is the same for all elements
        if e > 0:
            previous_elements_block_size = 0
            for n in range(e):
                previous_elements_block_size = previous_elements_block_size + element_properties[2][el_names[n]] +1
            for k in range(0, previous_elements_block_size):
                row.append(0)
        for k in range(dim_block):
            row.append(1)
        for k in range(column_skip_counter+dim_block-1, dim_v-1):
            row.append(0)
        M.append(row)
        row = []
        #column_skip_counter = column_skip_counter + 1


    #Finally append the electron density row
    row = []
    for element in el_names:
        dim_block = element_properties[2][element] + 1 #calculate it again
        for i in range(dim_block):
            row.append(i)
    M.append(row)
    row = []

    #remains to create the right hand side vector containing the total densities and one electron density
    col = []
    basic_counter = 0
    for element in el_names:
        dim_block = element_properties[2][element] + 1  # calculate it again!
        for i in range(dim_block-1):
            col.append(0)
        col.append(el_dens[basic_counter])
        basic_counter = basic_counter + 1
    col.append(nE)

    #Now we add this column to the matrix @@ complete please
    for l in range(len(M)):
        M[l].append(col[l])

    #Print the matrix (only in testing phase)@@
    for l in M:
        print(l)

    return M

#This function solves an augmented matrix containing an electron density symbol(sympy)
def matrix_solver_for_ne(M):

    #First we define this function and expression to get the electron density from optimisation
    expr = M.det()
    def eq(ne):
        return expr.subs(nE, ne)

    # Maximum electron density to limit the optimiser
    max_ne = @@
    min_ne = 0

    # Solving for electron density
    sol = root_scalar(eq, method='brentq', bracket=(min_ne, max_ne))

    # Subsitute into the matrix
    M.subs(nE, sol.root)

    # Convert it now to a numpy matrix and take the system part(not augmented) out
    MM = np.array(np.array(M.subs(nE, sol.root)), np.float64)
    # print(MM)
    M1 = MM[0:5, 0:5]
    M2 = MM[0:5, 5]

    # The solution in form of (nh0,nh1,nhe0,nhe1,nhe2)
    xx = np.linalg.solve(M1, M2)

    # Now write to the export file (with some of the imported data for ease of use)
    imported_data[i, 0:8].tofile(data_export, sep=' ', format='%s')
    data_export.write(" ")
    xx.tofile(data_export, sep=' ', format='%s')
    data_export.write("\n")

    # Printing the percentage of calculations done
    print((i * 100) / row)