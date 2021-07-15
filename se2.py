#Importing
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from sympy import *
import math


# Global Constants
h =  2*math.pi
kb = 1
me = 511e3 #[eV]
#Ionisation energies (eV)
xi_H1  = 13.6
xi_He1 = 24.6
xi_He2 = 54.4
#Partition Functions (taken from independent solutions for
# each element, so might not be as obvious as they are here...)
zh0  = 2
zh1  = 1
zhe0 = 1
zhe1 = 2
zhe2 = 1

#Define the right side of the Saha Equation, which depends on 
# the Temperature [T], the difference in energies of the states [xi]
# and the fraction of the partition functions G(r+1)/G(r) [z0,z1]
def lam(T,xi,z0,z1):
    return (z1/z0) * ((2*math.pi*me*kb*T)**1.5 / (h**3)) * math.exp((-xi) / (kb*T))

#Defining a temperature vector
N = 50
Temps = [];
for i in range(N):
    Temps.append((i+1000)*100000)

data_export = open('results2.txt', 'a')

#Loop over all the rows
for i in range(50000):
    #Get the value for that radius of the Sun
    nH = (0.36209)*4.31013646*(10**18)/(0.93878405573*10**9)
    nHe = (0.62194)*4.31013646*(10**18)/(0.372784341*10**(10))
    T = (Temps[i])/11604.51812

    #Making the electron density into a symbolic variable
    nE = symbols('n_e')

    #The matrix representing the system of equations and their value
    M = Matrix([[-lam(T,xi_H1,zh0,zh1), nE,                        0,                        0,   0,   0],
                [                    0,  0, -lam(T,xi_He1,zhe0,zhe1),                       nE,   0,   0],
                [                    0,  0,                        0, -lam(T,xi_He2,zhe1,zhe2),  nE,   0],
                [                    1,  1,                        0,                        0,   0,  nH],
                [                    0,  0,                        1,                        1,   1, nHe],
                [                    0,  1,                        0,                        1,   2,  nE]])


    #Building this equation in order to take it to zero later to find the value of n_e
    expr = M.det()
    def eq(ne):
        return expr.subs(nE,ne)
    
    #Maximum electron density to limit the optimiser
    max_ne = nH + 2*nHe;
    #print("The maximum electron density:")
    #print(max_ne)

    #minimum is obvious
    min_ne = 0

    #Solving for electron density
    sol = root_scalar(eq, method='brentq', bracket=(min_ne,max_ne))
    #print(sol)
    #print("\nThe ionisation percentage:")
    #print((sol.root * 100)/max_ne)

    #Subsitute into the matrix
    M.subs(nE,sol.root)

    #Convert it now to a numpy matrix and take the system part(not augmented) out
    MM = np.array(np.array(M.subs(nE,sol.root)), np.float64)
    #print(MM)
    M1 = MM[0:5,0:5]
    M2 = MM[0:5,5]
    #The solution in form of (nh0,nh1,nhe0,nhe1,nhe2)
    xx = np.linalg.solve(M1, M2)

    #Now write to the export file
    xx.tofile(data_export, sep=' ', format='%s')
    data_export.write("\n")

    print((i*100)/N)

print("Done")