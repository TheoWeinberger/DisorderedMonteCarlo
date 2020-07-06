#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to produce 3D lattice plots to help for visualisation of spin behaviour in real space                #
# produces plots of spins for flippability, number of neighbours as well as sub plots of a cubic system     #
# for ease of visualisation                                                                                 #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Reads in CSV files containing Lattice Spin data, Lattice Spin Flip Data, Lattice Neighbour data           #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# produces 4 plots, 2 of flippability denoted by lattice size and spin direction on a full lattice and      #
# cubic lattice plot and 2 of lattice neighbours and flippability                                           #                            
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
matplotlib.rcParams.update({'font.size': 60})
matplotlib.rcParams['axes.labelpad'] = 70
plt.rcParams["font.family"] = "Times New Roman"
plt.set_cmap("jet")

#import data to be plotted
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_LatticeSpinsT05.csv')
Lattice = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Spin'])
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_LatticeT05.csv')
LatticeSize = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Spin Flips'])
data = pd.read_csv(r'C:\Users\Theo\PartIII\FINAL_NearestNeighbourT05.csv')
LatticeNeighbours = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Neighbours'])

#assign gif and ax for full lattice plot
fig = plt.figure()
ax = plt.axes(projection="3d")

#convert lattice indexing to cartesian coordinates
x_coord = (Lattice['i']+Lattice['j'])
y_coord = (Lattice['k']+Lattice['i'])
z_coord = (Lattice['k']+Lattice['j'])

#assing data sets to determine colour,size etc.
Neighbours = LatticeNeighbours['Neighbours']
Spin_Size = 10000*(LatticeSize['Spin Flips']+0.05)
Spin = np.abs(Lattice['Spin']) #use this for spin colour if just trying to differentiate between spins and non spins

#list to contain nearest neighbour number data to create a colour map from
Spin_Colour = []

#go through data set and assign either 0 or 1 for spins with either even or odd neighbours for the colourmap, assign the value 3 for non-spins
for i in range(len(Neighbours)):
    if  Spin[i] != 0:
        Spin_Colour.append(Neighbours[i]%2)
    else:
        Spin_Colour.append(3)

#plot data
ax.scatter3D(x_coord,y_coord,z_coord, s = Spin_Size, c = Spin_Colour)

#add text labels of number of neighbours
#for i in range(len(Neighbours)): 
#    ax.text(x_coord[i],y_coord[i],z_coord[i],  '%d' % (int(Neighbours[i])), size=20, zorder=1,  
#    color='k') 


# draw cube
r = [0,2,4,6,8]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax.plot3D(*zip(s, e), color="b", linewidth = 0.5)

#square up axes
ax.set_aspect("equal")

plt.show()

#Assing axes for the 2x2 cubic unit cell plot 
ax2 = plt.axes(projection="3d")

#lists to contain subset of coordinates and colour/size data for the 2x2 cell
x_coord_sub = []
y_coord_sub = []
z_coord_sub = []
Spin_Size_sub = []
Spin_Colour_sub = []

#append valid data values for that list
for i in range(4,9,1):
    for j in range(4,9,1):
        for k in range(4,9,1):
            for n in range(len(x_coord)):
                if x_coord[n] == i and y_coord[n] == j and z_coord[n] == k:
                    x_coord_sub.append((x_coord[n]-4)/2)
                    y_coord_sub.append((y_coord[n]-4)/2)
                    z_coord_sub.append((z_coord[n]-4)/2)
                    Spin_Size_sub.append(Spin_Size[n])
                    Spin_Colour_sub.append(Spin_Colour[n])

#plot data
ax2.scatter3D(x_coord_sub,y_coord_sub,z_coord_sub, s = 1000*Spin_Size_sub, c = Spin_Colour_sub)
ax2.set_xlabel('x (l.u.)')
ax2.set_ylabel('y (l.u.)')
ax2.set_zlabel('z (l.u.)')


ax2.set_aspect("equal")

#plot background cubes for aif of visualisation
r = [0,1,2]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax2.plot3D(*zip(s, e), color="b", linewidth = 0.5)

plt.show()        

#assing fig and axes for sublattice plot of full lattice data
fig1 = plt.figure()
ax1 = plt.axes(projection="3d")

#import data for sublattice to be used to pick required spins
data = pd.read_csv(r'C:\Users\Theo\PartIII\AFMSublattice.csv')
SubLattice = pd.DataFrame(data, columns = ['i', 'j', 'k', 'Spin'])

#booleans to be used to select data
Spin_Coords = Lattice['Spin'] != 0 
Fixed = LatticeSize['Spin Flips'] < 0.1

#use booleans to choose subsets of original data sets required for the sublattice plots
LatticeNeighbours = LatticeNeighbours[Spin_Coords & Fixed]
LatticeSize = LatticeSize[Spin_Coords & Fixed]
SubLattice = SubLattice[Spin_Coords & Fixed]
Spin_Coords = Lattice[Spin_Coords & Fixed]

#convert from FCC indexing to cartesian coords
x_coord = (Spin_Coords['i']+Spin_Coords['j'])/2
y_coord = (Spin_Coords['k']+Spin_Coords['i'])/2
z_coord = (Spin_Coords['k']+Spin_Coords['j'])/2

#colour & size arrays
Neighbours = LatticeNeighbours['Neighbours']
Spin_Size = 10000*(LatticeSize['Spin Flips']+0.05)*SubLattice['Spin']
Spin = Spin_Coords['Spin'] #use this for spin colour if just trying to differentiate between spins and non spins

#plot data
ax1.scatter3D(x_coord,y_coord,z_coord, s = Spin_Size, c = Spin)

#text labels for number of nearest neighbours
#for i in range(len(Neighbours)): 
#    ax.text(x_coord[i],y_coord[i],z_coord[i],  '%d' % (int(Neighbours[i])), size=20, zorder=1,  
#    color='k') 

#square plot
ax1.set_aspect("equal")
ax1.set_xlabel('x (l.u.)')
ax1.set_ylabel('y (l.u.)')
ax1.set_zlabel('z (l.u.)')

plt.show()

#assign axes for 2x2 cubic cell plot of sublattice
ax2 = plt.axes(projection="3d")

#lists to contain subset of coordinates and colour/size data for the 2x2 cell
x_coord_sub = []
y_coord_sub = []
z_coord_sub = []
Spin_sub = []
Spin_Colour_sub = []

#assign valid data 
for i in range(4,11,2):
    for j in range(4,11,2):
        for k in range(4,11,2):
            for n in range(len(x_coord)):
                if x_coord.iloc[n] == i and y_coord.iloc[n] == j and z_coord.iloc[n] == k:
                    x_coord_sub.append((x_coord.iloc[n]-4)/2)
                    y_coord_sub.append((y_coord.iloc[n]-4)/2)
                    z_coord_sub.append((z_coord.iloc[n]-4)/2)
                    Spin_Size_sub.append(Spin_Size.iloc[n])
                    Spin_sub.append(Spin.iloc[n])

#plot data
ax2.scatter3D(x_coord_sub,y_coord_sub,z_coord_sub, s = 1000*Spin_Size_sub, c = Spin_sub)

#square data
ax2.set_aspect("equal")
ax2.set_xlabel('x (l.u.)')
ax2.set_ylabel('y (l.u.)')
ax2.set_zlabel('z (l.u.)')

#create background cube
r = [0,1,2,3]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax2.plot3D(*zip(s, e), color="b", linewidth = 0.5)

plt.show() 
