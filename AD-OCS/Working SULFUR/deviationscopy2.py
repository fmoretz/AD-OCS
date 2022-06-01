import numpy as np
import pandas as pd
from Influent import*
import matplotlib.pyplot as plt

t_span = np.linspace(1,200,1000) # DAYS!!
t_change_vett = T3.index.values  # Get the time when the deviations are applied
y_influent = np.zeros((len(t_span), len(y_in_0))) # Defines the y_in values which will be used in the main

index = 0
t_change_loc = t_change_vett[index]                # Get the time when the first deviation is applied
scale_loc    = T3.loc[t_change_loc].to_numpy() # Get the scale values for the first time change

i = 0
t = t_span[i]
while True:
    if t < t_change_loc:        
        y_influent[i] = y_in_0*scale_loc
        if i == (len(t_span)-1):  
            break
        else:      
            i += 1
            t = t_span[i]
    else:
        index = index+1
        t_change_loc = t_change_vett[index]
        scale_loc    = T3.loc[t_change_loc].to_numpy()
        print(f'*** CHANGE: {t_change_loc} ***')
    
        
def deviation_check(t,original):
    '''Defines the deviation to be used in the ODE integration'''
    deviated = np.zeros(len(original)) 
    i = 0
    while True: 
        if t < t_change_vett[i]:
            t_change = t_change_vett[i]
            scale = T3.loc[t_change].to_numpy()
            break
        else:
            i = i+1
            
          
    deviated = original*scale           
    return deviated