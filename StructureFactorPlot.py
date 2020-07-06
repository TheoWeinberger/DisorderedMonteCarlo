#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to plot structure factor in a 2D plane as desired by the user                                        #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# CSV file containing the full 3D structure factor data, must specificy L the system size                   #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# plot of structure factor in desired plane as well as a 3D plot of points were the structure factor is     #
# larger than a threshold                                                                                   #                            
#############################################################################################################
#############################################################################################################

import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import math as mth
from scipy.interpolate import griddata
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

#formatting parameters
mpl.rcParams.update({'font.size': 60})
plt.rcParams["font.family"] = "Times New Roman"

# activate latex text rendering
rc('text', usetex=True)

#reciprocal lattice system size
L = 40

#import data
data = pd.read_csv(r'C:\Users\Theo\PartIII\StructureFactor8.csv')
Lattice = pd.DataFrame(data, columns = ['i', 'j', 'k', 'S'])

#array to store coordinates
HHK_Array_List = []

#reciprocal coordinates
q1_coord = Lattice['i']
q2_coord = Lattice['j']
q3_coord = Lattice['k']

#import structure factor data
Structure_Factor = pd.Series(Lattice['S'])

#specificy plane to be studied here
for i in range(len(q1_coord)):
    if (q1_coord[i] == q2_coord[i]):
        Structure_Factor_Coords = np.array([q2_coord[i], q3_coord[i],Structure_Factor[i]])
        HHK_Array_List.append(Structure_Factor_Coords)
HHK_Array = np.array(HHK_Array_List)

#meshgrid of structure factor data
X, Y = np.meshgrid((2/(L))*HHK_Array[:,0], (2/(L))*HHK_Array[:,1])
Z = griddata(((2/(L))*HHK_Array[:,0], (2/(L))*HHK_Array[:,1]), HHK_Array[:,2], (X, Y), method='cubic')
Z = Z/np.max(HHK_Array[:,2])

#plot data
fig, ax = plt.subplots()

#plot data for each brillouin zone with translations
ax.pcolormesh(X,Y,Z, cmap = 'rainbow')
ax.pcolormesh(X-np.max(X),Y,Z, cmap = 'rainbow')
ax.pcolormesh(X,Y-np.max(Y),Z,  cmap = 'rainbow')
ax.pcolormesh(X-np.max(X),Y-np.max(Y),Z,  cmap = 'rainbow')

#colormap data for colorbar
pcm = ax.pcolormesh(X, Y, Z, cmap='rainbow')

#plot colorbar
fig.colorbar(pcm,ax=ax)

#square data
ax.set_aspect("equal")

#set eaxis labels
ax.set_ylabel("(00k) r.l.u.")
ax.set_xlabel("(hh0) r.l.u.")

plt.show()

#plot data
fig1 = plt.figure()
ax1 = plt.axes(projection="3d")

#plot data for each brillouin zone with translations

Peaks = Lattice['S']/(Lattice['S'].max()) > 0.7 #define threshold for the structure factor here
Peaks = Lattice[Peaks]

#plot 3D data
ax1.scatter3D((2/L)*Peaks['i'],(2/L)*Peaks['j'],(2/L)*Peaks['k'], c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i']-2,(2/L)*Peaks['j'],(2/L)*Peaks['k'], c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i'],(2/L)*Peaks['j']-2,(2/L)*Peaks['k'], c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i'],(2/L)*Peaks['j'],(2/L)*Peaks['k']-2, c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i']-2,(2/L)*Peaks['j']-2,(2/L)*Peaks['k'], c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i']-2,(2/L)*Peaks['j'],(2/L)*Peaks['k']-2, c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i'],(2/L)*Peaks['j']-2,(2/L)*Peaks['k']-2, c = Peaks['S'])
ax1.scatter3D((2/L)*Peaks['i']-2,(2/L)*Peaks['j']-2,(2/L)*Peaks['k']-2, c = Peaks['S'])


#square plot
ax1.set_aspect("equal")
ax1.set_xlabel(r'$k_{1}$ (r.l.u.)')
ax1.set_ylabel(r'$k_{2}$ (r.l.u.)')
ax1.set_zlabel(r'$k_{3}$ (r.l.u.)')

plt.show()




