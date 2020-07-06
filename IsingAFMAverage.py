#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to average over thermodynanmic data for AFM lattice and produce error analysis                       #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in all csv files from a folder 'IsingFCCAFMData' to average over histories                      #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Averaged data of e, m, sublattice magnetisations, specific heat, susceptibility and relevant errors       #                            
#############################################################################################################
#############################################################################################################

import pandas as pd
import numpy as np
import glob

main_path = r'C:\Users\Theo\PartIII\IsingFCCAFMData'


#list to contain data for individual histories
Histories_List = []


#average data over histories
all_files = glob.glob(main_path + "/*.csv")

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    Histories_List.append(df)

Temperature = 0
for i in range(len(Histories_List)):
    Temperature += Histories_List[i]['Temperature']/len(Histories_List)

Average_Energy = 0
for i in range(len(Histories_List)):
    Average_Energy += Histories_List[i]['Average Energy']/len(Histories_List)

Average_Spin = 0
for i in range(len(Histories_List)):
    Average_Spin += Histories_List[i]['Average Spin']/len(Histories_List)

Average_Spin_Susc = 0
for i in range(len(Histories_List)):
    Average_Spin_Susc += np.absolute(Histories_List[i]['Average Spin'])/len(Histories_List)

Average_M1 = 0
for i in range(len(Histories_List)):
    Average_M1 += np.sign(Histories_List[i]['M1'].iloc[-1])*Histories_List[i]['M1']/len(Histories_List)

Average_M2 = 0
for i in range(len(Histories_List)):
    Average_M2 += np.sign(Histories_List[i]['M2'].iloc[-1])*Histories_List[i]['M2']/len(Histories_List)

Average_M3 = 0
for i in range(len(Histories_List)):
    Average_M3 += np.sign(Histories_List[i]['M3'].iloc[-1])*Histories_List[i]['M3']/len(Histories_List)

Average_M4 = 0
for i in range(len(Histories_List)):
    Average_M4 += np.sign(Histories_List[i]['M4'].iloc[-1])*Histories_List[i]['M4']/len(Histories_List)

Average_Energy_Squared = 0
for i in range(len(Histories_List)):
    Average_Energy_Squared += Histories_List[i]['Average Energy']**2/len(Histories_List)

Average_Spin_Squared = 0
for i in range(len(Histories_List)):
    Average_Spin_Squared += Histories_List[i]['Average Spin']**2/len(Histories_List)

Average_M1_Squared = 0
for i in range(len(Histories_List)):
    Average_M1_Squared += Histories_List[i]['M1']**2/len(Histories_List)

N = 0
for i in range(len(Histories_List)):
    N += Histories_List[i]['N']/len(Histories_List)

#calculate observables

Specific_Heat = N*(Average_Energy_Squared-Average_Energy**2)/Temperature**2

Susceptibility = N*(Average_Spin_Squared-Average_Spin_Susc**2)/Temperature


#calculate standard deviations

E_SD = np.sqrt(Average_Energy_Squared-Average_Energy**2)

M_SD = np.sqrt(Average_Spin_Squared-Average_Spin**2)

M1_SD = np.sqrt(Average_M1_Squared-Average_M1**2)

SH_SD = N*np.sqrt(2*E_SD**4/len(Histories_List))/Temperature**2

Susc_SD = N*np.sqrt(2*M_SD**4/len(Histories_List))/Temperature


#concatenate data

History_Averaged = pd.concat([Temperature,Average_Energy,E_SD,Average_Spin,M_SD,Average_M1,M1_SD,Average_M2,Average_M3,Average_M4,Specific_Heat,SH_SD,Susceptibility,Susc_SD], axis = 1)

History_Averaged.columns = ['Temperature','Average Energy','Energy SD','Average Spin','Spin SD','M1','M1 SD','M2','M3','M4','Specific Heat','Specific Heat SD','Susceptibility','Susceptibility SD']


#Output data

Output_Data =  History_Averaged

Output_Data.to_csv(r'C:\Users\Theo\PartIII\FINAL_AFM_L20.csv', index = False)

