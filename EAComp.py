#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to plot and compare sublattice magnetisation data against the EA order parameter                     #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing thermodynamic data for disordered lattice                                #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces plots of Q and M1 against temperature and the overal susceptibiltiy                              #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np

#formatting parameters
matplotlib.rcParams.update({'font.size': 40})
plt.rcParams["font.family"] = "Times New Roman"

#read in data
data = pd.read_csv(r'C:\Users\Theo\PartIII\DisorderdataL15H50D20EA.csv')


#plot average spin and Q 
Plot_AverageSpin = pd.DataFrame(data, columns = ['Temperature', 'Average Spin'])
Plot_Q = pd.DataFrame(data, columns = ['Temperature', 'Average Q'])
M = Plot_AverageSpin['Average Spin']
Q = Plot_Q['Average Q']
T = Plot_AverageSpin['Temperature']
plt.scatter(T, M, c = 'k', marker = '+', s = 200, label = 'M')
plt.scatter(T, Q, c = 'k', marker = 'x', label = 'Q', s = 100)
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature", fontsize = 45) 
plt.ylabel("Average Spin", fontsize = 45) 
plt.legend()
plt.show()


#plot susceptibilities
Plot_Susceptibility = pd.DataFrame(data, columns = ['Temperature', 'Susceptibility'])
Plot_Susceptibility_Q = pd.DataFrame(data, columns = ['Temperature', 'Susceptibility Q'])
Susceptibility = Plot_Susceptibility['Susceptibility']
Susceptibility_Q = Plot_Susceptibility_Q['Susceptibility Q']
T = Plot_Susceptibility['Temperature']
plt.scatter(T, Susceptibility, c = 'k', marker = '+', s = 200, label = 'Susceptibility')
plt.scatter(T, Susceptibility_Q, c = 'k', marker = 'x', label = 'Susceptibility Q', s = 100)
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature", fontsize = 45) 
plt.ylabel("Susceptibility", fontsize = 45) 
plt.legend()
plt.show()

