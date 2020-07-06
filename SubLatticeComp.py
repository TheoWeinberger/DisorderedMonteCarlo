#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to plot sublattice magnetisations for AFM v disordered case as well as Q & M1, and Q and Susc Q      #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing thermodynamic data for AFM lattice and Disordered lattice                #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces comparison plots for the thermodynamic variables                                                 #                            
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
matplotlib.rcParams.update({'errorbar.capsize': 4}) 
plt.rcParams["font.family"] = "Times New Roman"

#read in data to be plotted
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_AFM_L20.csv')
Disorderdata = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_Disorder_L40.csv')

#Define separate data frames containing raw data and Beath data for average energy
Plot_M1 = pd.DataFrame(data, columns = ['Temperature', 'M1'])
Plot_M1Dis = pd.DataFrame(Disorderdata, columns = ['Temperature', 'Average M1'])

#plot raw data as an errorbar plot and Beath data as a line plot
plt.plot(Plot_M1['Temperature'], Plot_M1['M1'], color = 'black', linestyle = '(0,(5,10))', linewidth = 2, label = 'Pure AFM Lattice')
plt.scatter(Plot_M1Dis['Temperature'], Plot_M1Dis['Average M1'], color = 'black', marker = '+', s = 300, label = 'Disordered AFM Lattice')

#plot formatting
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel(r"Sublattice Magnetisation") 
plt.legend(fontsize = 40)
plt.show()

#Define separate data frames containing raw data and Beath data for average energy
Plot_Q = pd.DataFrame(Disorderdata, columns = ['Temperature', 'Susceptibility Q','Average Q'])

#plot raw data as an errorbar plot and Beath data as a line plot
plt.scatter(Plot_Q['Temperature'], Plot_Q['Average Q'], color = 'black', marker = 'x', label = 'Q', s = 150)
plt.scatter(Plot_M1Dis['Temperature'], Plot_M1Dis['Average M1'], color = 'black', marker = '+', s = 300, label = 'M1')

#plot formatting
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.xlabel(r"Temperature ($|J|/k_B$)") 
plt.ylabel(r"Q/M1") 
plt.legend(fontsize = 40)
plt.show()

fig, ax1 = plt.subplots()

ax1.set_xlabel(r"Temperature ($|J|/k_B$)")
ax1.set_ylabel("Q")
ax1.scatter(Plot_Q['Temperature'], Plot_Q['Average Q'], color = 'black', marker = 'x', label = 'Q', s = 150)
ax1.tick_params(axis='y')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel(r"$\chi_Q$")  # we already handled the x-label with ax1
ax2.set_ylim(-0.04,2.8)
ax2.scatter(Plot_Q['Temperature'], Plot_Q['Susceptibility Q'], color = 'black', marker = '+', s = 300, label = r'$\chi_Q$')
ax2.tick_params(axis='y')

#fig.tight_layout()  # otherwise the right y-label is slightly clipped

#plot formatting
ax1.grid(linestyle='--', linewidth='0.5', color='gray')
fig.legend(loc = 'upper right', bbox_to_anchor=(0.9, 0.88), fontsize = 40)
#fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()