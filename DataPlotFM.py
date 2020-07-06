#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to plot thermodynamic data for FM ising data                                                         #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing thermodynamic data for FM lattice                                        #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces plots of C, m, Susceptibility,e                                                                  #                            
#############################################################################################################
#############################################################################################################


import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from matplotlib import rc
# activate latex text rendering
rc('text', usetex=True)


#formatting parameters
matplotlib.rcParams.update({'font.size': 40})
matplotlib.rcParams.update({'errorbar.capsize': 2}) 
plt.rcParams["font.family"] = "Times New Roman"

#read in data to be plotted
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_FM_L20.csv')

#plot specific heat on two different scales
Plot_SpecificHeat = pd.DataFrame(data, columns = ['Temperature', 'Specific Heat', 'Specific Heat SD'])
Plot_SpecificHeat['log(T)'] = np.log10(Plot_SpecificHeat['Temperature'])

Plot_SpecificHeat.plot.scatter(x = 'Temperature', y = 'Specific Heat', yerr = 'Specific Heat SD', marker = '+', s = 100, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature", fontsize = 45) 
plt.ylabel("Specific Heat", fontsize = 45) 
plt.show()

Plot_SpecificHeat.plot.scatter(x = 'log(T)',y = 'Specific Heat', yerr = 'Specific Heat SD', marker = '+', s = 100, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("log(Temperature)", fontsize = 45) 
plt.ylabel("Specific Heat", fontsize = 45) 
plt.show()

#plot average spin
Plot_AverageSpin = pd.DataFrame(data, columns = ['Temperature', 'Average Spin', 'Spin SD'])
Plot_AverageSpin.plot.scatter(x = 'Temperature', y = 'Average Spin', yerr = 'Spin SD', marker = '+', s = 100, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature", fontsize = 45) 
plt.ylabel("Average Spin", fontsize = 45) 
plt.show()

#plot susceptibility on two different scales
Plot_Susceptibility = pd.DataFrame(data, columns = ['Temperature', 'Susceptibility', 'Susceptibility SD'])
Plot_Susceptibility.plot.scatter(x = 'Temperature', y = 'Susceptibility', yerr = 'Susceptibility SD', marker = '+', s = 100, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature", fontsize = 45) 
plt.ylabel("Susceptibility", fontsize = 45) 
plt.show()

Plot_Susceptibility['log(T)'] = np.log10(Plot_SpecificHeat['Temperature'])
Plot_Susceptibility.plot.scatter(x = 'log(T)', y = 'Susceptibility', yerr = 'Susceptibility SD', marker = '+', s = 100, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("log(Temperature)", fontsize = 45) 
plt.ylabel("Susceptibility", fontsize = 45) 
plt.show()

#plote average energy
Plot_AverageEnergy = pd.DataFrame(data, columns = ['Temperature', 'Average Energy', 'Energy SD'])
Plot_AverageEnergy.plot.scatter(x = 'Temperature', y = 'Average Energy', yerr = 'Energy SD', marker = '+', s = 100, color='black')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature", fontsize = 45) 
plt.ylabel("Average Energy", fontsize = 45) 
plt.show()