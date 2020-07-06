#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to produce histogram plots of spin flippability for whole spin system and a 3D histogram of the      #
# of the data against spin neighbour number. Also produces a scatter plot of flippability against neigbours #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Reads in CSV files containing Lattice Spin data, Lattice Spin Flip Data, Lattice Neighbour data           #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# produces 2 histograms of spin flippability and a scatter plot of flippability against neighbour number    #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import math as mth
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)

#set font rc parameters
matplotlib.rcParams.update({'font.size': 60})
#matplotlib.rcParams['axes.labelpad'] = 70 #use this for the 3D plots
plt.rcParams["font.family"] = "Times New Roman"

#choose number of bins for the histograms here
Num_Bins = 20

#import all datasets
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_LatticeSpinsT05.csv')
Lattice = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Spin'])
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_LatticeT05.csv')
LatticeSpinFlips = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Spin Flips'])
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_NearestNeighbourT05.csv')
LatticeNeighbours = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Neighbours'])
#data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_TotalNeighbourSpinT05.csv')
#SpinTotal = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Neighbour Spin'])

#list for histogram data to be stored in
HistData = []

#go through LatticeSpinFlips data frame and get data for all spins 
for i in range(len(LatticeSpinFlips)):
    if Lattice.loc[i, "Spin"] != 0:
        HistData.append(LatticeSpinFlips.loc[i, "Spin Flips"])

#assign figure & subplots
fig, ax = plt.subplots()

#Histogram parameters 
HistFitting, edges = np.histogram(HistData, bins = 20, range =  [0,1])
HistCentre = edges[:-1] + (edges[1] - edges[0])/2

#function to be fitted to
def func(x, p1):
    return x**-p1

#fit parameters
Const_Opt, Const_Cov = curve_fit(func, HistCentre, HistFitting)

#output fitting parameters
print(Const_Opt)

#continuous curve for comparison
Hist_Linspace = np.linspace(0, 1, 100)
Hist_Fit = func(Hist_Linspace, Const_Opt[0])

#plot histogram data
ax.plot(Hist_Linspace,Hist_Fit, label = r"$y=x^{-2.17}$", color = 'black', linestyle = (0,(5,10)), linewidth = 2)

#plot display add ins
ax.hist(HistData, bins = 20, color = 'lightsteelblue', edgecolor = 'black')
ax.set_yscale('log')
ax.set_ylabel('Number of Spins')
ax.set_xlabel('Flippability')
ax.legend(prop={'size': 40})
plt.show()

#lists for 3d histogram data to be stored in
Neighbours = []
Spinflips = []
NeighbourSpin = []

#go through data set and append desired data to lists for the non-zero spins
for i in range(len(LatticeSpinFlips)):
    if Lattice.loc[i, "Spin"] != 0:
        Neighbours.append(LatticeNeighbours.loc[i, "Neighbours"])
        Spinflips.append(LatticeSpinFlips.loc[i, "Spin Flips"])
        #NeighbourSpin.append(SpinTotal.loc[i, "Neighbour Spin"])

#plot a scatterplot of the data
plt.scatter(Neighbours,Spinflips, c = 'k', marker = '+', s = 300 )
plt.ylabel('Flippability')
plt.xlabel('Number of neighbouring spins')
plt.grid(linestyle='--', linewidth='0.5', color='gray')
plt.show()

#plt.scatter(NeighbourSpin,Spinflips, c = 'k', marker = '+' )
#plt.ylabel('Flipability')
#plt.xlabel('Total Spin on site')
#plt.show()

#assign fig and ax for the 3d histogram plot
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')

#assign regions of histogram to be used by 3d barplot to produce a histogram
hist, Nedges, Fedges = np.histogram2d(Neighbours, Spinflips, bins = (12, Num_Bins), range = [[0,12],[0,1]])
N, F = np.meshgrid(Nedges[:-1], Fedges[:-1], indexing="ij")
N, F = N.ravel(), F.ravel()
Z = 0
dN = 0.5 * np.ones_like(Z)
dF = 1/Num_Bins * np.ones_like(Z)
dz=hist.ravel()
dz = np.log10(dz+1)

#create array of colours so each row has a different colour
colours = np.zeros((12,Num_Bins))
for i in range(12):
    colours[i] = i

#flatten and colourmap the array
colours = plt.cm.jet(colours.flatten()/float(colours.max()))

#plot data
ax1.bar3d(N,F,Z,dN,dF,dz,zsort='average', color = colours, alpha = 0.5)
ax1.set_xlabel('Number of Neighbours')
ax1.set_ylabel('Flippability')
ax1.set_zlabel('log(Number of Spins)')
plt.show()
