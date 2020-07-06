#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to plot my raw AFM/FM ising model data against the data from the Beath et al model to test how well  #
# my monte carlo simulation works                                                                           #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# AFM/FM final data and AFM/FM Beath data from CSV file                                                     #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# produces plots of raw AFM/FM data for the average energy and susceptiblity as a comparison/benchmark for  #
# my monte carlo simulation                                                                                 #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)

#set plot parameters
matplotlib.rcParams.update({'font.size': 60})
matplotlib.rcParams.update({'errorbar.capsize': 6}) 
plt.rcParams["font.family"] = "Times New Roman"

#read in data to be plotted
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_FM_L20.csv')
BeathData = pd.read_csv(r'C:\Users\Theo\PartIII\BeathFM.csv')

#Define separate data frames containing raw data and Beath data for susceptibility
Plot_Susceptibility = pd.DataFrame(data, columns = ['Temperature', 'Susceptibility', 'Susceptibility SD'])
Plot_SusceptibilityBeath = pd.DataFrame(BeathData, columns = ['T Susceptibility', 'Susceptibility'])

#plot raw data as an errorbar plot and Beath data as a line plot
plt.errorbar(Plot_Susceptibility['Temperature'], Plot_Susceptibility['Susceptibility'], yerr = Plot_Susceptibility['Susceptibility SD'], fmt = 'o', marker = '+', color='black', markersize = 15)
plt.plot(np.array(Plot_SusceptibilityBeath['T Susceptibility']), np.array(Plot_SusceptibilityBeath['Susceptibility']), color = 'black', markersize = 15 )

#plot formatting
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel("Susceptibility") 
plt.show()



#Define separate data frames containing raw data and Beath data for average energy
Plot_AverageEnergy = pd.DataFrame(data, columns = ['Temperature', 'Average Energy', 'Energy SD'])
Plot_AverageEnergyBeath = pd.DataFrame(BeathData, columns = ['T Energy', 'Energy'])

#plot raw data as an errorbar plot and Beath data as a line plot
plt.errorbar(Plot_AverageEnergy['Temperature'], Plot_AverageEnergy['Average Energy'], yerr = Plot_AverageEnergy['Energy SD'], fmt = 'o', marker = '+', color='black', markersize = 15)
plt.plot(np.array(Plot_AverageEnergyBeath['T Energy']), np.array(Plot_AverageEnergyBeath['Energy']), color = 'black', markersize = 15)

#plot formatting
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel(r"Average Energy ($|J|$)") 
plt.show()
