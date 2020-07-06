#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to average over thermodynanmic data for disordered lattice and produce error analysis                #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in all csv files from a folder 'DisordersEASubLattices' to average over histories and           #
# and disorders.                                                                                            #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# Averaged data of C, m, Q, Susceptibility, Q Susceptibility ,e, sublattice magnetisation as well as        #
# relevant error analysis                                                                                   #                            
#############################################################################################################
#############################################################################################################

import pandas as pd
import glob
import numpy as np

#code to average over disorders and histories

main_path = r'C:\Users\Theo\PartIII\DisordersEASubLattices' # path to data
all_files_main = glob.glob(main_path + "/*") #array of paths to data

#list to store individual history averages of disordered data in
Disorder_List = []

#for each disorder avergae over histories
for i in range(len(all_files_main)):

    #list of histories for each disroder
    Histories_List = []

    #path to data
    all_files = glob.glob(all_files_main[i] + "\\data" + "/*.csv")

    #read in data
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        Histories_List.append(df)

    #average over histories

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

    Average_Q = 0
    for i in range(len(Histories_List)):
        Average_Q += Histories_List[i]['Q']/len(Histories_List)
    
    Average_Energy_Squared = 0
    for i in range(len(Histories_List)):
        Average_Energy_Squared += Histories_List[i]['Average Energy']**2/len(Histories_List)

    Average_Spin_Squared = 0
    for i in range(len(Histories_List)):
        Average_Spin_Squared += Histories_List[i]['Average Spin']**2/len(Histories_List)

    Average_M1_Squared = 0
    for i in range(len(Histories_List)):
        Average_M1_Squared += Histories_List[i]['M1']**2/len(Histories_List)

    Average_Q_Squared = 0
    for i in range(len(Histories_List)):
        Average_Q_Squared += Histories_List[i]['Q']**2/len(Histories_List)

    N = 0
    for i in range(len(Histories_List)):
        N += Histories_List[i]['N']/len(Histories_List)

    #calculate observables

    Specific_Heat = N*(Average_Energy_Squared-Average_Energy**2)/Temperature**2

    Susceptibility = N*(Average_Spin_Squared-Average_Spin_Susc**2)/Temperature

    Susceptibility_Q = N*(Average_Q_Squared-Average_Q**2)/Temperature

    E_VAR = Average_Energy_Squared-Average_Energy**2

    M_VAR = Average_Spin_Squared-Average_Spin**2

    Q_VAR = Average_Q_Squared - Average_Q**2

    M1_VAR = Average_M1_Squared-Average_M1**2

    SH_VAR = (N**2)*(2*E_VAR**2/len(Histories_List))/Temperature**4

    Susc_VAR = (N**2)*(2*M_VAR**2/len(Histories_List))/Temperature**2

    Q_Susc_VAR = (N**2)*(2*Q_VAR**2/len(Histories_List))/Temperature**2

    #concatenate and title data into a new dataframe
    History_Averaged = pd.concat([Temperature,Average_Energy,E_VAR,Average_Spin,M_VAR,Average_M1,M1_VAR,Average_M2,Average_M3,Average_M4,Average_Q,Q_VAR,Specific_Heat,SH_VAR,Susceptibility,Susc_VAR,Susceptibility_Q,Q_Susc_VAR], axis = 1)

    History_Averaged.columns = ['Temperature','Average Energy','Energy SD','Average Spin','Spin SD','Average M1','M1 SD','Average M2','Average M3','Average M4','Average Q','Q SD','Specific Heat','Specific Heat SD','Susceptibility','Susceptibility SD','Susceptibility Q','Susceptibility Q SD']

    #append data to disorder list
    Disorder_List.append(History_Averaged)


#average over disorder
Temperature = 0
for i in range(len(Disorder_List)):
    Temperature += Disorder_List[i]['Temperature']/len(Disorder_List)

Average_Energy_Disorder = 0
for i in range(len(Disorder_List)):
    Average_Energy_Disorder += Disorder_List[i]['Average Energy']/len(Disorder_List)

Energy_SD = 0
for i in range(len(Disorder_List)):
    Energy_SD += Disorder_List[i]['Energy SD']/len(Disorder_List)

Average_Spin_Disorder = 0
for i in range(len(Disorder_List)):
    Average_Spin_Disorder += Disorder_List[i]['Average Spin']/len(Disorder_List)

Spin_SD = 0
for i in range(len(Disorder_List)):
    Spin_SD += Disorder_List[i]['Spin SD']/len(Disorder_List)

Average_M1_Disorder = 0
for i in range(len(Disorder_List)):
    Average_M1_Disorder += Disorder_List[i]['Average M1']/len(Disorder_List)

M1_SD = 0
for i in range(len(Disorder_List)):
    M1_SD += Disorder_List[i]['M1 SD']/len(Disorder_List)

Average_M2_Disorder = 0
for i in range(len(Disorder_List)):
    Average_M2_Disorder += Disorder_List[i]['Average M2']/len(Disorder_List)

Average_M3_Disorder = 0
for i in range(len(Disorder_List)):
    Average_M3_Disorder += Disorder_List[i]['Average M3']/len(Disorder_List)

Average_M4_Disorder = 0
for i in range(len(Disorder_List)):
    Average_M4_Disorder += Disorder_List[i]['Average M4']/len(Disorder_List)

Average_Q_Disorder = 0
for i in range(len(Disorder_List)):
    Average_Q_Disorder += Disorder_List[i]['Average Q']/len(Disorder_List)

Q_SD = 0
for i in range(len(Disorder_List)):
    Q_SD += Disorder_List[i]['Q SD']/len(Disorder_List)

Specific_Heat_Disorder = 0
for i in range(len(Disorder_List)):
    Specific_Heat_Disorder += Disorder_List[i]['Specific Heat']/len(Disorder_List)

SH_SD = 0
for i in range(len(Disorder_List)):
    SH_SD += Disorder_List[i]['Specific Heat SD']/len(Disorder_List)

Susceptibility_Disorder = 0
for i in range(len(Disorder_List)):
    Susceptibility_Disorder += Disorder_List[i]['Susceptibility']/len(Disorder_List)

Susceptibility_SD = 0
for i in range(len(Disorder_List)):
    Susceptibility_SD += Disorder_List[i]['Susceptibility SD']/len(Disorder_List)

Susceptibility_Q_Disorder = 0
for i in range(len(Disorder_List)):
    Susceptibility_Q_Disorder += Disorder_List[i]['Susceptibility Q']/len(Disorder_List)

Susceptibility_Q_SD = 0
for i in range(len(Disorder_List)):
    Susceptibility_Q_SD += Disorder_List[i]['Susceptibility Q SD']/len(Disorder_List)


#Convert variances to standard deviations

Energy_SD = pd.DataFrame(np.sqrt(Energy_SD))

Spin_SD = pd.DataFrame(np.sqrt(Spin_SD))

M1_SD = pd.DataFrame(np.sqrt(M1_SD))

Q_SD = pd.DataFrame(np.sqrt(Q_SD))

SH_SD = pd.DataFrame(np.sqrt(SH_SD))

Susceptibility_SD = pd.DataFrame(np.sqrt(Susceptibility_SD))

Susceptibility_Q_SD = pd.DataFrame(np.sqrt(Susceptibility_Q_SD))


#output data averaged over disorders
Output_Data = pd.concat([Temperature,Average_Energy_Disorder,Energy_SD,Average_Spin_Disorder,Spin_SD,Average_M1_Disorder,M1_SD,Average_M2_Disorder,Average_M3_Disorder,Average_M4_Disorder,Average_Q_Disorder,Q_SD,Specific_Heat_Disorder,SH_SD,Susceptibility_Disorder,Susceptibility_SD,Susceptibility_Q_Disorder,Susceptibility_Q_SD], axis = 1)

Output_Data.to_csv(r'C:\Users\Theo\PartIII\FINAL_Disorder_L40.csv', index = False)

