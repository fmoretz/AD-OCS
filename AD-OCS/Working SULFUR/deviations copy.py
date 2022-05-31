import numpy as np
from Influent import*
import matplotlib.pyplot as plt

t_span = np.linspace(1,200,1000) # DAYS!!

y_influent = np.zeros(len(t_span), len(y_in_0))) # Defines the y_in values which will be used in the main

for i in range(len(t_span)):
    t = t_span[i]
    if t < 20 or t > 100: 
        y_influent[i] = y_in_0           
    else:
        for j in range(len(y_in_0)):
            y_influent[i,j] = y_in_0[j]
            y_influent[i,4] = 1.2*y_in_0[4]   # Assigns deviation of 1.2*XT to influent

print(y_influent.shape)
def deviation_check(t,original):
    '''Defines the deviation to be used in the ODE integration'''
    deviated = np.zeros(len(original)) 
    if t < 20 or t > 100:
        deviated[0] = original[0]
        deviated[1] = original[1]
        deviated[2] = original[2]
        deviated[3] = original[3]
        deviated[4] = original[4]                
    else:         
        deviated[0] = original[0]
        deviated[1] = original[1]
        deviated[2] = original[2]
        deviated[3] = original[3]
        deviated[4] = 1.2*original[4]   # Assigns deviation of 1.2*XT to influent
    return deviated