#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to read in the SG and AFM structure factors along a 1D direction and fit the form of the structure   #
# factor to them in order to determine the SG and AFM correlation lengths allowing for comparison of the    #
# magnetic behaviour types within the system                                                                #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing averaged structure factor along a direction: CorrelationLengths____.csv  #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces plots of averaged structure factor along with their fitting function if desired. Calculates      #
# correlation length as a function of temperature and plots as desired                                      #                            
#############################################################################################################
#############################################################################################################



import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import math as mth
from matplotlib import rc
from scipy.optimize import curve_fit

#formatting parameters
mpl.rcParams.update({'font.size': 60})
plt.rcParams["font.family"] = "Times New Roman"
# activate latex text rendering
rc('text', usetex=True)

#reciprocal lattice system size
L = 40

#import data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_CORR_111.csv')
StructureFactor = pd.DataFrame(data, columns = ['T', 'i', 'S', 'S_SG'])

#function to be fitted to to determine p3, correlation length
def func(x, p1, p2, p3):
    return p1/(x**(1.7)+p3**(-2.3))

#produce list of temperatures studied 
T_List = StructureFactor["T"]
T_List = np.unique(T_List)

T_Plot = []
AFM_Correl = []
SG_Correl = []

for i in range(int(len(T_List)/2.2)):
    #get values from StructureFactor at a given temperature
    StructureFactorT = StructureFactor.loc[StructureFactor['T'] == T_List[i]]

    #get indices of max value in the spin glass and afm structure factor
    S_indices_max = np.where(StructureFactorT["S"] == StructureFactorT["S"].max())
    S_SG_indices_max = np.where(StructureFactorT["S_SG"] == StructureFactorT["S_SG"].max())

    S = StructureFactorT["S"].iloc[S_indices_max[0].max():len(StructureFactorT["S"])] #structure factor for afm from bragg peak to lowest value
    S_SG = StructureFactorT["S_SG"].iloc[S_SG_indices_max[0].min():int(L/2)] #structure factor data for spin glass from bragg peak to lowest value

    #k coordinates for AFM correlation length
    k_S = StructureFactorT["i"].iloc[S_indices_max[0].max():len(StructureFactorT["S"])]*(2/(L)) #coordinates in required data range
    k_S_zeroed = k_S - k_S.min() #zero coordinates

    #k coordinated for G correlation length
    k_SG = StructureFactorT["i"].iloc[S_SG_indices_max[0].min():int(L/2)]*(2/(L)) #coordinates in required data range
    k_SG_zeroed = k_SG - k_SG.min() #zero coordinates

    #k coordinates
    k = StructureFactorT["i"]*(2/(L))

    #fit curve to data for AFM data
    Const_Opt, Const_Cov = curve_fit(func, k_S_zeroed, S)
    AFM_Correl.append(Const_Opt[2])

    #plt.plot(k,StructureFactorT["S"])
    #plt.plot(k_S,func(k_S_zeroed, Const_Opt[0], Const_Opt[1], Const_Opt[2]))
    #plt.show()

    #fit curve to data for Spin Glass data
    Const_Opt, Const_Cov = curve_fit(func, k_SG_zeroed, S_SG)
    SG_Correl.append(Const_Opt[2])

    #plt.plot(k,StructureFactorT["S_SG"])
    #plt.plot(k_SG,func(k_SG_zeroed, Const_Opt[0], Const_Opt[1], Const_Opt[2]))
    #plt.show()

    T_Plot.append(T_List[i])

#plt.scatter(T_Plot,AFM_Correl, c = 'k', marker = '+')
#plt.xlabel('T')
#plt.ylabel(r'$\xi_{AFM}$ (l.u.)')
#plt.show()
#plt.scatter(T_Plot,SG_Correl, c = 'k', marker = '+')
#plt.xlabel('T')
#plt.ylabel(r'$\xi_{SG}$ (l.u.)')
#plt.show()

plt.scatter(T_Plot ,AFM_Correl, c = 'k', marker = 'x', s = 150, label = r'$\xi_{AFM}$')
plt.xlabel(r'T $|J|/k_B$')
plt.ylabel(r'$\xi$ (l.u.)')
plt.scatter(T_Plot,SG_Correl, c = 'k', marker = '+', s = 300, label = r'$\xi_{SG}$')
plt.legend()
plt.show()

#plt.plot(T_Plot,SG_Correl)
#plt.show()













