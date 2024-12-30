# Saha Equation in Many-Ion Plasmas

The Saha equation is a thermodynamic realtionship between the ionisation levels of an element, their energy differences and temperature.
This equation can be used to compare the ratios of the multiple ionisation levels of given elements.

The goal of this repository is to provide a python solver for the Saha equation that can take a state (total densities
of the elements as well as the temperature of the system) and return the exact densities of each ionised state of these
elements.

# Unit Convention
Given that the code will be mainly used for solar studies, the input data will be assumed to be in CGS units (cm-g-s) and the temperature is assumed in Kelvin.
The code will then convert it into natural units (hbar = c = kb = 1). It is important to stick to these units because
some of the computations include mass which is provided in the code dictionary in natural units (1Da = 9.3289e8 eV). 

Here are some shortcuts to convert:

1eV = 1.78e-33 g ; 1eV^-1 = 1.97e-5 cm

1g = 0.5618e33 eV ; 1cm = 0.5076e5 eV^-1

1 g/cm^3 = 4.2955e18 eV^4

# Installation

Here the task is simple. You just have to download the code and have the appropriate packages installed on your machine. These packages are listed in the `virtual_env.yml` file along with their versions. 

If you use Anaconda/Miniconda, you can automate this process by using the env create command in this folder:

`conda env create -f virtual_env.yml`

Once this is done the code should run without any issues.

# Usage
### Default Usage
#### Input:
The default way to run the code is to give it the elements you are interested in, their densities, and the temperature of the system. This can be done directly from the terminal as follows

`python main.py -e H He -d 1662424176 719087680 -t 1.549e+09`

The example above is solving a system of Hydrogen and Helium, with respective densities of 1662424176 and 719087680 eV^4, at a temperature of 1.549e+09 eV.

#### Output:

The output here is descriptive (given that its just one computation). 

`Temperature: 133482.49224845882 eV`\
`Total Number Densities: [7.595260026104078e+18, 8.272280373370961e+17] eV3`\
`Number Densities of Ionised States: [7.540700734348527e+18, 5.4559291755552696e+16, 8.037987768016884e+17, 2.3261011402521604e+16, 168249132885872.66] eV3`\
`Maximum Electron Number Density: 9.24971610077827e+18 eV3`\
`Actual Electron Number Density: 7.815680142384608e+16 eV3`\
`The total ionization in the system is: 0.8 %`\
`--- 0.15461421012878418 seconds ---`

### Solar Model
#### Input:
Another way of running this program is by using a Solar Model file (or any other model). Such a file should contain rows of different temperatures and densities. A sample model can be found here https://wwwmpa.mpa-garching.mpg.de/~aldos/SSM/AGSS09/model_agss09.dat (now migrated to https://aliga.ice.csic.es/personal/aldos/Solar_Data_files/struct_b16_agss09.dat)

To run the code you have to still give it the element symbols you want, with the filename of the model and the position of the temperature and respective element columns. Here is an example

`python main.py -e H He -s solar_data.txt -p 2 6 7 -f results.txt`

Which means that we want to solve for Hydrogen and Helium, and that the third column is the temperature data, and the 7th and 8th are the respective densities of the elements. We also include the name of the file to save the results to.

#### Output:

The command output will be just the time elapsed. The actual solver results will go into the file with the specified name above. They will be comma seperated and in the following order

`[T(eV), nH(eV3), nHe(eV3), ..., nH0(eV3), nH1(eV3), nHe0(eV3), ..., nE_Max(eV3), nE(eV3)]`

# Available elements
Below is the list of items included in the code dictionary, and their corresponding symbol.
* Hydrogen (H)
* Helium (He)
* Helium 3 Isotope (He3)
* Carbon (C)
* Nitrogen (N)
* Oxygen (O)
* Neon (Ne)
* Silicon (Si)
* Iron (Fe)

# Issues and errors

If you face any issues or errors while working with this code, do not hesitate to post it on "https://github.com/jelkuweiss/many-ion-saha-equation/issues" and I will get back to it ASAP.

# Reference

To cite the material in this repository, you can simply use the widget in the about section. It will generate a citation automatically.
