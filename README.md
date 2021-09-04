# Saha Equation in Many-Ion Plasmas

The Saha equation is a thermodynamic realtionship between the ionisation levels of an element, their energy differences and temperature.
This equation can be used to compare the ratios of the multiple ionisation levels of given elements.

The goal of this repository is to provide a python solver for the Saha equation that can take a state (total densities of the elements as well as the temperature of the system) and return the exact densities of each ionised state of these elements.
# Usage
### Default
#### Input:
The default way to run the code is to give it the elements you are interested in, their densities, and the temperature of the system. This can be done directly from the terminal as follows

`python main.py -e H He -d 1662424176 719087680 -t 1.549e+09`

The example above is solving a system of Hydrogen and Helium, with respective densities of 1662424176 and 719087680 eV^4, at a temperature of 1.549e+09 eV.

#### Output:

The output here is simply a line in the terminal with the densities in the order of input, and least to most ionised (exp: H0, H1, He0, He1, He2)

`[7.29086779e+00 1.66242417e+09 3.45778562e-09 3.15369993e+00 7.19087677e+08]`

### Solar Model
#### Input:
Another way of running this program is by using a Solar Model file (or any other model). Such a file should contain rows of different temperatures and densities. A sample model can be found here https://wwwmpa.mpa-garching.mpg.de/~aldos/SSM/AGSS09/model_agss09.dat

To run the code you have to still give it the element symbols you want, with the filename of the model and the position of the temperature and respective element columns. Here is an example

`python main.py -e H He -s solar_data.txt -p 3 7 8`

Which means that we want to solve for Hydrogen and Helium, and that the third column is the temperature data, and the 7th and 8th are the respective densities of the elements.

#### Output:

By default the code will output the columns that you specified from the model file, followed by the ionisation densities in the order explained above, into a .txt file names results. For now this can only be changed by altering the code manually. 