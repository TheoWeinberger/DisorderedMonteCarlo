#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to calculate entropy change in the system. Involves numerical integration of available data as well  #
# as producing high temeperature curve fitting to extend integral to infinity                               #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing thermodynamic data for lattice                                           #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces plots entropy aaginst temperature as well as showing residual entropy and error in calculation   #                            
#############################################################################################################
#############################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from scipy import integrate
import pandas as pd
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)

#plot formatting
matplotlib.rcParams.update({'font.size': 60})
matplotlib.rcParams.update({'errorbar.capsize': 6}) 
plt.rcParams["font.family"] = "Times New Roman"

#Attempt to fit curve to the heat capacity which can subsequently be integrated to find entropy

def Entropy(Specific_Heat_T):
    #find index of when temperature is above a threshold for which fitting applies
    i = np.argmax(T<2) # specify temperature range here ~11 for FM case and ~2 for disordered case

    #fit high temperature tail of curve to a 1/T^2 relationship
    #function to be fitted to
    def func(x, p1, p2):
        return p1/(x+p2)**2

    #fit parameters
    Const_Opt, Const_Cov = curve_fit(func, T[:i], Specific_Heat[:i])

    #continuous curve for comparison
    T_Linspace = np.linspace(2, 6, 100)
    Specific_Heat_Fit = func(T_Linspace, Const_Opt[0], Const_Opt[1])

    #plot data
    #Plot_SpecificHeat.plot.scatter(x = 'Temperature', y = 'Specific Heat', marker = '+', s = 100, color='black')
    #plt.grid(linestyle='--', linewidth='0.5', color='g')
    #plt.xlabel("Temperature", fontsize = 45) 
    #plt.ylabel("Specific Heat", fontsize = 45) 
    #plt.plot(T_Linspace, Specific_Heat_Fit)
    #plt.show()

    #plt.plot(T, Specific_Heat_T)
    #plt.plot(T_Linspace, Specific_Heat_Fit/T_Linspace)
    #plt.grid(linestyle='--', linewidth='0.5', color='g')
    #plt.xlabel("Temperature", fontsize = 45) 
    #plt.ylabel("Specific Heat", fontsize = 45) 
    #plt.show()

    #function for the tail to be integrated
    tail = lambda y:  (Const_Opt[0]*(1/(y+Const_Opt[1]))**2)/y 

    #calculate entropy change from largest point of temperature measured to infinity
    Delta_S_Infinity = integrate.quad(tail, T[0], np.inf)
    print(Delta_S_Infinity)

    #Calculate and plot entropy via trapezium rule on raw data. 
    Delta_S = np.log(2) - Delta_S_Infinity[0]
    Delta_S_Array = []
    Delta_S_Array.append(Delta_S)
    for i in range(len(T)-1):
        Area = (T[i]-T[i+1])*(Specific_Heat_T[i]+Specific_Heat_T[i+1])*0.5
        Delta_S_New = Delta_S_Array[i] - Area
        Delta_S_Array.append(Delta_S_New)

    return Delta_S_Array


#read in data
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_Disorder_L40.csv')
Plot_SpecificHeat = pd.DataFrame(data, columns = ['Temperature', 'Specific Heat','Specific Heat SD'])
T = Plot_SpecificHeat['Temperature']

#produce specific heat data for average, max and min cases
Specific_Heat = Plot_SpecificHeat['Specific Heat'] #+ Plot_SpecificHeat['Specific Heat SD']
Specific_Heat_T = Specific_Heat/T

Specific_Heat_Plus = Plot_SpecificHeat['Specific Heat'] + Plot_SpecificHeat['Specific Heat SD']
Specific_Heat_T_Plus = Specific_Heat_Plus/T

Specific_Heat_Minus = Plot_SpecificHeat['Specific Heat'] - Plot_SpecificHeat['Specific Heat SD']
Specific_Heat_T_Minus = Specific_Heat_Minus/T

#total change in entropy arrays
Delta_S_Array = Entropy(Specific_Heat_T)

Delta_S_Array_Plus = Entropy(Specific_Heat_T_Plus)

Delta_S_Array_Minus = Entropy(Specific_Heat_T_Minus)

Delta_S_Array = np.array(Delta_S_Array)

Delta_S_Array_Plus = np.array(Delta_S_Array_Plus)

Delta_S_Array_Minus = np.array(Delta_S_Array_Minus)

#calculate +/- errors
Delta_S_Error_Plus =  Delta_S_Array - Delta_S_Array_Plus

Delta_S_Error_Minus = Delta_S_Array_Minus - Delta_S_Array 

#concatenate error arrays
Errors = np.vstack((Delta_S_Error_Plus,Delta_S_Error_Minus))

#output data
print(Delta_S_Array)
print(Errors)

#plot data
plt.errorbar(T, Delta_S_Array, yerr = Delta_S_Error_Plus, c = 'k', markersize = 15)
plt.xlabel(r'Temperature $(|J|/k_B$)')
plt.ylabel(r'$\Delta$ S $(k_B)$')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.show()

