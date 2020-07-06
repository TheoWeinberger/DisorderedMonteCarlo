#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to average over thermodynanmic data for FM lattice and produce error analysis                        #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in all csv files from a folder 'IsingFCCFMData' to average over histories                       #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Averaged data of e, m, specific heat, susceptibility and relevant errors                                  #                            
#############################################################################################################
#############################################################################################################

import pandas as pd
import numpy as np
import glob

main_path = r'C:\Users\Theo\PartIII\IsingFCCFMData' 

#list to contain the histories data to be averaged over in
Histories_List = []

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
    Average_Spin += np.sign(Histories_List[i]['Average Spin'].iloc[-1])*Histories_List[i]['Average Spin']/len(Histories_List)

Average_Spin_Susc = 0
for i in range(len(Histories_List)):
    Average_Spin_Susc += np.absolute(Histories_List[i]['Average Spin'])/len(Histories_List)

Average_Energy_Squared = 0
for i in range(len(Histories_List)):
    Average_Energy_Squared += Histories_List[i]['Average Energy']**2/len(Histories_List)

Average_Spin_Squared = 0
for i in range(len(Histories_List)):
    Average_Spin_Squared += Histories_List[i]['Average Spin']**2/len(Histories_List)

N = 0
for i in range(len(Histories_List)):
    N += Histories_List[i]['N']/len(Histories_List)


#calculate observables

Specific_Heat = N*(Average_Energy_Squared-Average_Energy**2)/Temperature**2

Susceptibility = N*(Average_Spin_Squared-Average_Spin_Susc**2)/Temperature

#calculate errors

E_SD = np.sqrt(Average_Energy_Squared-Average_Energy**2)

M_SD = np.sqrt(Average_Spin_Squared-Average_Spin**2)

SH_SD = N*np.sqrt(2*E_SD**4/len(Histories_List))/Temperature**2

Susc_SD = N*np.sqrt(2*M_SD**4/len(Histories_List))/Temperature


#retitle data

History_Averaged = pd.concat([Temperature,Average_Energy,E_SD,Average_Spin,M_SD,Specific_Heat,SH_SD,Susceptibility,Susc_SD], axis = 1)

History_Averaged.columns = ['Temperature','Average Energy','Energy SD','Average Spin','Spin SD','Specific Heat','Specific Heat SD','Susceptibility','Susceptibility SD']

#Output data

Output_Data =  History_Averaged

Output_Data.to_csv(r'C:\Users\Theo\PartIII\FINAL_FM_L20.csv', index = False)

