import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def f_RR_equilibrium(alpha, *args):
    # Unpack the arguments
    z, law, P_sat, H, P = args
    RR = np.zeros(len(z))  
    P_sp = np.zeros(len(z))
    K = np.zeros(len(z))
    # Calculate the G/L equilibrium
    for i in range(len(z)):
        if law == 'R':
            P_sp[i] = P_sat[i]
        else:
            P_sp[i] = H[i] 
        K[i] = P_sp[i]/P

        RR[i] = z[i]*(K[i]-1)/(1+alpha*(K[i]-1))

    RR_eqn = sum(RR)
    return RR_eqn