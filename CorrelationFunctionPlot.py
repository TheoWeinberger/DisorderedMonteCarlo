#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to read in real space correlation functions and to plot the real space correlation of spins for      #
# the desired normalised value of correlation in order to determine the structure of the real space         #
# correlations to complement the study of the structure factor                                              #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in a csv file containing real space correlation data: CFDisorder_______.csv                     #
# Must choose the degree of correlation of which spins are plotted                                          #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Produces a 3D lattice plot of the spin that are correlated within the system                              #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math as mth
from itertools import product, combinations
from matplotlib import rc


# activate latex text rendering
rc('text', usetex=True)

#set rc parameters for plotting
matplotlib.rcParams.update({'font.size': 40})
matplotlib.rcParams['axes.labelpad'] = 40
plt.rcParams["font.family"] = "Times New Roman"
plt.set_cmap("jet")

#import data to be plotted
data = pd.read_csv(r'C:\Users\Theo\PartIII\CFDisorderT06001.csv')
Lattice = pd.DataFrame(data, columns = ['i', 'j', 'k', 'C'])

Spin_Coords = np.abs(Lattice['C']) > 0.85 #degree of correlation of the spins to be plotted. choose this parameter.

Lattice = Lattice[Spin_Coords] #reduce lattice to the subset of spin that are correlated as these are the ones to be plotted.

#assign fig and ax for full lattice plot
fig = plt.figure()
ax = plt.axes(projection="3d")

#convert lattice indexing to cartesian coordinates
x_coord = (Lattice['i']+Lattice['j'])/2
y_coord = (Lattice['k']+Lattice['i'])/2
z_coord = (Lattice['k']+Lattice['j'])/2


#assing data sets to determine colour,size etc.
Correlation = np.array(Lattice['C']) #use this for spin colour if just trying to differentiate between spins and non spins

#list to contain nearest neighbour number data to create a colour map from
Spin_Colour = []

#go through data set and assign either 0 or 1 for spins with either even or odd neighbours for the colourmap, assign the value 3 for non-spins
for i in range(len(Correlation)):
    if  Correlation[i] > 0:
        Spin_Colour.append(1)
    elif  Correlation[i] < 0:
        Spin_Colour.append(2)    
    else:
        Spin_Colour.append(3)

#plot data
ax.scatter3D(x_coord,y_coord,z_coord, c = Spin_Colour, s = 100)

# draw cube. 
#r = [9,10,11,12] # change these values depending on which spins are correlated for clarity. may be necessary to remove this section. 
#for s, e in combinations(np.array(list(product(r, r, r))), 2):
#    if np.sum(np.abs(s-e)) == r[1]-r[0]:
#        ax.plot3D(*zip(s, e), color="b", linewidth = 0.5)

#square up axes
ax.set_aspect("equal")

ax.set_xlabel('x (l.u.)')
ax.set_ylabel('y (l.u.)')
ax.set_zlabel('z (l.u.)')

plt.show()
