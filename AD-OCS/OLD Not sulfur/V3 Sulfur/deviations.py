import numpy as np
import pandas as pd
from Influent import*


t_span = np.linspace(15,200,1000) # DAYS!!
t_change_vett = T3.index.values  # Get the time at when the deviations are applied
y_influent = np.zeros((len(t_span), len(y_in_0))) # Defines the y_in values which will be used in the main

index = 0
t_change_loc = t_change_vett[index]                # Get the time when the first deviation is applied
scale_loc    = np.ones(len(y_in_0))   # Get the scale values for the first time change


i = 0
t = t_span[i]

while True:
    if t < t_change_loc:        
        y_influent[i] = y_in_0*scale_loc
        print(t, t_change_loc, scale_loc)
        if i == (len(t_span)-1):  
            break
        else:      
            i += 1
            t = t_span[i]
    else:        
        scale_loc    = T3.loc[t_change_loc].to_numpy()
        print(f'*** CHANGE at: {t_change_loc} ***')
        if index == (len(t_change_vett)-1):
            t_change_loc = np.inf
            print('Last Change done')
        else:
            index = index+1        
            t_change_loc = t_change_vett[index]         
                
    
        
def deviations_check(t, original, t_change_vett):
    '''Defines the deviation to be used in the ODE integration'''
    deviated = np.zeros(len(original)) 
    i = 0    
    while True: 
        if t < t_change_vett[0]: # if the time is before the first time change then the original values are used
            scale = 1
            break
        if t >= t_change_vett[-1]: # if the time is after the last time change then the last values are used
            t_change = t_change_vett[-1]
            scale = T3.loc[t_change].to_numpy()
            break
        if t > t_change_vett[i] and t < t_change_vett[i+1]:
            t_change = t_change_vett[i]
            scale = T3.loc[t_change].to_numpy()
            break
        else:
            i = i+1
            
          
    deviated = original*scale           
    return deviated

prova = deviations_check(175, y_in_0, T3.index.values)
print(prova)