#Importing
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from sympy import *
import math

#List of dictionaries for element properties
#first dictionary contains the statistical weight of each ionised state of the atom
#Second one is the energy difference xi between the states
#Third one is the number of ionised states of each element
element_properties = [
    {
        "H" : [2,1],
        "He": [1,2,1],
    },
    {
        "H" : [13.6],
        "He": [24.6,54.4],
    },
    {
        "H" : 1,
        "He": 2,
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

    for element in el_names:

        dim_block = element_properties[2][element] + 1 #size of the block matrix of each element
        block_block_skip_counter = 0 #skip counter for each element block matrix (they are mostly upper triangular)

        for i in range(dim_block - 1): #the -1 is because the last line of the block is just a list of 1's
            column_skip_counter = column_skip_counter + 1

            #first we do all the skipping we have to do to reach the place of our element block
            for k in range(0,column_skip_counter):
                row.append(0)

            #Then we build the block lines: first by skipping enough to make the block upper triangular
            for k in range(column_skip_counter,column_skip_counter+block_block_skip_counter):
                row.append(0)

            #Then we place the two elements
            row.append(-lam(T, element_properties[1][element][i], element_properties[0][element][i], element_properties[0][element][i+1]))
            row.append(nE)

            #fill the rest with zeros
            for k in range(column_skip_counter+block_block_skip_counter + 2, dim_v):
                row.append(0)

            #Append to the matrix, clear the row, and do it again
            M.append(row)
            row.clear()

        #Now the last line in the block is the total density line which is the same for all elements
        for j in range(dim_v):
            row.append(1)
            M.append(row)
            row.clear()

                
    #Finally append the electron density row
    row = []
    for element in el_names:
        dim_block = element_properties[2][element] + 1 #calculate it again
        for i in range(dim_block):
            row.append(i)
    M.append(row)
    row.clear()

    #Print the matrix (only in testing phase)@@
    print(M)

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