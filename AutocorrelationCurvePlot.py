#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to read in autocorrelation curve data and plot the observed decay behaviour for in the disordered    #
# system to see the temperature dependence of the autocorrelation function to make qualitative observations #
# from                                                                                                      #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in 4 csv files at the desired temperatures for study: CorrelationT_.csv                         #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces a plot of the time dependence of the autcorrelation function at 4 different temperatures         #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from matplotlib import rc

#code to plot the correlation curves at 4 temperature investigated

#set plot parameters
matplotlib.rcParams.update({'font.size': 60})
plt.rcParams["font.family"] = "Times New Roman"

# activate latex text rendering
rc('text', usetex=True)

#read in T=4 data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_CorrelationT4.csv')
Plot_T4 = pd.DataFrame(data, columns = ['t_mc', 'C'])

#obtain non-zero data
t_mc4 = Plot_T4.loc[Plot_T4['C'] != 0.0, 't_mc']
CorrelationT4 = Plot_T4.loc[Plot_T4['C'] != 0.0, 'C']

#read in T=2 data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_CorrelationT2.csv')
Plot_T2 = pd.DataFrame(data, columns = ['t_mc', 'C'])

#obtain non-zero data
t_mc2 = Plot_T2.loc[Plot_T2['C'] != 0.0, 't_mc']
CorrelationT2 = Plot_T2.loc[Plot_T2['C'] != 0.0, 'C']

#read in T=1 data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_CorrelationT1.csv')
Plot_T1 = pd.DataFrame(data, columns = ['t_mc', 'C'])

#obtain non-zero data
t_mc1 = Plot_T1.loc[Plot_T1['C'] != 0.0, 't_mc']
CorrelationT1 = Plot_T1.loc[Plot_T1['C'] != 0.0, 'C']

#read in T=0.5 data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_CorrelationT05.csv')
Plot_T05 = pd.DataFrame(data, columns = ['t_mc', 'C'])

#obtain non-zero data
t_mc05 = Plot_T05.loc[Plot_T05['C'] != 0.0, 't_mc']
CorrelationT05 = Plot_T05.loc[Plot_T05['C'] != 0.0, 'C']

#plot all curves on same graph
plt.plot(t_mc4, CorrelationT4, label = r'T = 4 $|J|/k_{B}$', linewidth = '2.5')
plt.plot(t_mc2, CorrelationT2, label = r'T = 2 $|J|/k_{B}$', linewidth = '2.5')
plt.plot(t_mc1, CorrelationT1, label = r'T = 1 $|J|/k_{B}$', linewidth = '2.5')
plt.plot(t_mc05[:1500], CorrelationT05[:1500], label = r'T = 0.5 $|J|/k_{B}$', linewidth = '2.5')

plt.ylabel('Autocorrelation')
plt.xlabel(r'$\tau_{mc}$')

#show plot and legend
plt.legend(prop={'size': 40})
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.show()