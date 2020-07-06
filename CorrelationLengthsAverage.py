#############################################################################################################
#############################################################################################################
#                                            Author: 8230W                                                  #
#############################################################################################################
#                                           Code Description                                                #
#############################################################################################################
# Code to average over correlation lengths data (averaging over histories and disorders)                    #
#############################################################################################################
#                                               Inputs                                                      #
#############################################################################################################
# Must read in all csv files from a folder, CorrelationLengths_______ to be averaged over                   #
#############################################################################################################
#                                                Outputs                                                    #
#############################################################################################################
# produces 2 csv file outputs. One is the average over all histories and all disorders and one is the       #
# average over one disorder's histories                                                                     #                            
#############################################################################################################
#############################################################################################################


import pandas as pd
import glob
import numpy as np


main_path = r'C:\Users\Theo\PartIII\CorrelationLengths001' 
all_files_main = glob.glob(main_path + "/*")

#list to store each individual disorder in
Disorder_List = []

for i in range(len(all_files_main)):

    #list to store data for each history in 
    Histories_List = []

    #find all files
    all_files = glob.glob(all_files_main[i] + "\\data" + "/*.csv")

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        Histories_List.append(df)

    #average all data over histories

    T = 0
    for i in range(len(Histories_List)):
        T += Histories_List[i]['T']/len(Histories_List)

    i_index = 0
    for i in range(len(Histories_List)):
        i_index += Histories_List[i]['i']/len(Histories_List)

    i_index = round(i_index).astype(int)

    S = 0
    for i in range(len(Histories_List)):
        S += Histories_List[i]['S']/len(Histories_List)

    S_SG = 0
    for i in range(len(Histories_List)):
        S_SG += Histories_List[i]['S_SG']/len(Histories_List)


    #concatenate and title data into a new dataframe
    History_Averaged = pd.concat([T,i_index,S,S_SG], axis = 1)

    #retitle data columns
    History_Averaged.columns = ['T','i','S','S_SG']

    #append data to disorder list
    Disorder_List.append(History_Averaged)



# average all data over disorders

T = 0
for i in range(len(Disorder_List)):
    T += Disorder_List[i]['T']/len(Disorder_List)

i_index = 0
for i in range(len(Disorder_List)):
    i_index += Disorder_List[i]['i']/len(Disorder_List)

i_index = round(i_index).astype(int)

S = 0
for i in range(len(Disorder_List)):
    S += Disorder_List[i]['S'].astype(float)/len(Disorder_List)

S_SG = 0
for i in range(len(Disorder_List)):
    S_SG += Disorder_List[i]['S_SG'].astype(float)/len(Disorder_List)

Output_Data =  pd.concat([T,i_index,S,S_SG], axis = 1)

Output_Data_Single_Hist = Disorder_List[0] #can specify which history to be outputting here

Output_Data.to_csv(r'C:\Users\Theo\PartIII\FINAL_CORR_001.csv', index = False) #Output disorder averaged data
Output_Data_Single_Hist.to_csv(r'C:\Users\Theo\PartIII\FINAL_CORR_SINGLE_001.csv', index = False) #output history averaged data

    

    
