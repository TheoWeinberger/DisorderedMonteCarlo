#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code ro read in autocorrelation data and determine the temperature scaling parameter by fitting to        #
# the relevant exponential behaviour t = Ae^(B/T)                                                           #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in csv file containing desired autocorrlelation data: Correlationdata.csv                       #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces plot of autocorrelation data as well as data fitting. Print out the step scaling parameter       #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import rc


# activate latex text rendering
rc('text', usetex=True)

#code to plot autocorrelation function as a function of temperature an produce a fitting to it

matplotlib.rcParams.update({'font.size': 60})
plt.rcParams["font.family"] = "Times New Roman"

#read in data and assign to arrays
data = pd.read_csv(r'C:\Users\Theo\PartIII\Correlationdata.csv')
Plot = pd.DataFrame(data, columns = ['T_mc', 'T'])
T = Plot['T']
T_mc = Plot['T_mc']


#function to be fitted to
def func(x, p1, p2, p3):
    return p1*np.exp(p2/(x+p3)) 

#fit parameters
Const_Opt, Const_Cov = curve_fit(func, T, T_mc)

#print fitting parameters
print(Const_Opt)

#create fit data
T_Linspace = np.linspace(1.5, 3, 100)
T_mc_fit = func(T_Linspace, Const_Opt[0], Const_Opt[1], Const_Opt[2])

#plot data
Plot.plot.scatter(x = 'T', y = 'T_mc', marker = '+', s = 100, color='black', label = 'Autocorrelation Data')
plt.plot(T_Linspace, T_mc_fit, label = r'Fitted Curve: $\e^{}$')
plt.grid(linestyle='--', linewidth='0.5', color='g')
plt.xlabel("Temperature") 
plt.ylabel("Specific Heat") 
plt.legend()
plt.show()
