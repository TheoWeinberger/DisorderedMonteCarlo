#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to plot thermodynamic data for disordered ising data                                                 #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing thermodynamic data for disordered lattice                                #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces plots of C, m, Q, Susceptibility, Q Susceptibility ,e, sublattice magnetisation                  #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from matplotlib import rc
# activate latex text rendering
rc('text', usetex=True)

#code to plot data from disordered lattice simulation

matplotlib.rcParams.update({'font.size': 60})
plt.rcParams["font.family"] = "Times New Roman"
matplotlib.rcParams.update({'errorbar.capsize': 6}) 

#read in data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_Disorder_L40.csv')

#plot specific heat on two different scales
Plot_SpecificHeat = pd.DataFrame(data, columns = ['Temperature', 'Specific Heat', 'Specific Heat SD'])
Plot_SpecificHeat['log(T)'] = np.log10(Plot_SpecificHeat['Temperature'])

Plot_SpecificHeat.plot.scatter(x = 'Temperature', y = 'Specific Heat', yerr = 'Specific Heat SD', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel(r"Specific Heat ($k_B$)") 
plt.show()

Plot_SpecificHeat.plot.scatter(x = 'log(T)',y = 'Specific Heat', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel("log(Temperature)") 
plt.ylabel("Specific Heat") 
plt.show()

#plot average spin
Plot_AverageSpin = pd.DataFrame(data, columns = ['Temperature', 'Average Spin', 'Spin SD'])
Plot_AverageSpin.plot.scatter(x = 'Temperature', y = 'Average Spin', yerr = 'Spin SD', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel("Average Spin") 
plt.show()

#plot Edwards-Anderson Q parameter 
Plot_Q = pd.DataFrame(data, columns = ['Temperature', 'Average Q', 'Q SD'])
Plot_Q.plot.scatter(x = 'Temperature', y = 'Average Q', yerr = 'Q SD', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel("Q") 
plt.show()

#plot Suceptibility on two different scales
Plot_Susceptibility = pd.DataFrame(data, columns = ['Temperature', 'Susceptibility', 'Susceptibility SD'])
Plot_Susceptibility.plot.scatter(x = 'Temperature', y = 'Susceptibility', yerr = 'Susceptibility SD', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel("Susceptibility") 
plt.show()

Plot_Susceptibility['log(T)'] = np.log10(Plot_SpecificHeat['Temperature'])
Plot_Susceptibility.plot.scatter(x = 'log(T)', y = 'Susceptibility', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel("log(Temperature)") 
plt.ylabel("Susceptibility") 
plt.show()

#plot Edwards Anderson susceptibility
Plot_Susceptibility_Q = pd.DataFrame(data, columns = ['Temperature', 'Susceptibility Q', 'Susceptibility Q SD'])
Plot_Susceptibility_Q.plot.scatter(x = 'Temperature', y = 'Susceptibility Q', yerr = 'Susceptibility Q SD', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel("Susceptibility Q") 
plt.show()

#plot Average Energy
Plot_AverageEnergy = pd.DataFrame(data, columns = ['Temperature', 'Average Energy', 'Energy SD'])
Plot_AverageEnergy.plot.scatter(x = 'Temperature', y = 'Average Energy', yerr = 'Energy SD', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel(r"Average Energy ($|J|$)") 
plt.show()

#plot sublattice magnetisation of one sublattice
Plot_AverageM1 = pd.DataFrame(data, columns = ['Temperature', 'Average M1', 'M1 SD'])
Plot_AverageM1.plot.scatter(x = 'Temperature', y = 'Average M1', marker = '+', s = 200, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel("M1") 
plt.show()