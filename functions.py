#Importing
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from sympy import *
import math


#el_names is a list of names ['H','He']
#el_dens is a list of densities [5000,6000]
#T is the temperature
#ALL IN NATURAL UNITS
def solver(el_names,el_dens,T):
    #Known general data
    el_io_levels = {'H':1, 'He':2, 'C':3}
     
    #First we need to get the dimension of the vector in the matrix:
    #Get the total ionisation levels
    total_io_levels = 0
    for element in el_names:
        total_io_levels = total_io_levels + el_io_levels[element]
    #Dimesion of rows of augmented matrix
    dim_v = 1 + len(el_names) + total_io_levels
    
    #Now we build the Matrix
    M = []
    row = []
    column_skip_counter = 0
    for element in el_names:
        dim_block = el_io_levels[element]+1 #size of the block matrix of each element
        for i in range(dim_block):
            column_skip_counter = column_skip_counter + 1
            for j in range(dim_v):
                #first we do all the skipping we have to do
                for k in range(0,column_skip_counter):
                    row.append(0)
                for k in range(column_skip_counter,dim_block):
                    
                row.append@@
                
    #Finally append the electron density row            
    M.append@@

    #Print the matrix (only in testing phase)
    print(M)