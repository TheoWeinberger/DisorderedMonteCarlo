#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to average structure factor over histories and disorders to be analysed                              #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in all csv files from a folder containing structure factors                                     #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# 2 csv files containing averaged structure factors. One averaged over whole set of data and one of only one#
# disorder type                                                                                             #                            
#############################################################################################################
#############################################################################################################

import pandas as pd
import glob
import numpy as np

main_path = r'C:\Users\Theo\PartIII\StructureFactors1' 
all_files_main = glob.glob(main_path + "/*")

#list to contain disorder data in it
Disorder_List = []

for i in range(len(all_files_main)):

    #list to contain histories data 
    Histories_List = []

    all_files = glob.glob(all_files_main[i] + "\\data" + "/*.csv")

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        Histories_List.append(df)


    #average over data for one history 
    i_index = 0
    for i in range(len(Histories_List)):
        i_index += Histories_List[i]['i']/len(Histories_List)

    i_index = round(i_index).astype(int)

    j_index = 0
    for i in range(len(Histories_List)):
        j_index += Histories_List[i]['j']/len(Histories_List)

    j_index = round(j_index).astype(int)

    k_index = 0
    for i in range(len(Histories_List)):
        k_index += Histories_List[i]['k']/len(Histories_List)

    k_index = round(k_index).astype(int)

    S = 0
    for i in range(len(Histories_List)):
        S += Histories_List[i]['S']/len(Histories_List)


    #concatenate and title data into a new dataframe
    History_Averaged = pd.concat([i_index,j_index,k_index,S], axis = 1)

    History_Averaged.columns = ['i','j','k','S']

    #append data to disorder list
    Disorder_List.append(History_Averaged)



#average over disorders and convert indexes back to integers
i_index = 0
for i in range(len(Disorder_List)):
    i_index += Disorder_List[i]['i']/len(Disorder_List)

i_index = round(i_index).astype(int)

j_index = 0
for i in range(len(Disorder_List)):
    j_index += Disorder_List[i]['j']/len(Disorder_List)

j_index = round(j_index).astype(int)

k_index = 0
for i in range(len(Disorder_List)):
    k_index += Disorder_List[i]['k']/len(Disorder_List)

k_index = round(k_index).astype(int)

S = 0
for i in range(len(Disorder_List)):
    S += Disorder_List[i]['S'].astype(float)/len(Disorder_List)

#output data

Output_Data =  pd.concat([i_index,j_index,k_index,S], axis = 1)

Output_Data_Single_Hist = Disorder_List[3] #can specify which disorder number here

Output_Data.to_csv(r'C:\Users\Theo\PartIII\FINAL_SF_DISORDER_1.csv', index = False)
Output_Data_Single_Hist.to_csv(r'C:\Users\Theo\PartIII\FINAL_SF_DISORDER_SINGLE_1.csv', index = False)

    
