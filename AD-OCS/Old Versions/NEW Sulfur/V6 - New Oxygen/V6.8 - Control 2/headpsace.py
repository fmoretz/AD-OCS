
import numpy as np
from matplotlib import pyplot as plt
import time
from pprint import pprint as pp
from scipy.integrate import odeint

d_start = 0         # [h] - Start time
d_end   = 2      # [h] - End time
hours   = .002      # [h] - Discretization time

n_times = int((d_end-d_start)/hours)+1  # Number of time steps

print('***Intervals of {hours} hours. \n {n_times} time steps***'.format(hours=hours, n_times=n_times))

t = np.linspace(d_start,d_end,n_times) # time span
R_C = 1.8 # gS/gO2 reaction coefficient
P = 1 # atm
T = 35+273.15 # K
R = 0.0821 # L * atm/mol/K
Rgas_m3_atm_K = 0.08205746e-3 # m^3*atm/(mol*K)
V_gas = 400#m3
RC = 1.8                                                        # [gS/gO2]     - Stoichiometry of the reaction
alfa = 0.9                                                       # [-] - Power for sulfur concentration
beta = 0.2                                                      # [-] - Power for oxygen concentration
k= 0.6                                                     # [1/m3/h*(gS*gO)^proper exp] - Reaction rate constant

n_in = {
    'H2S': 4.52,
    'O2': 4.52/2,
    'SX':  0.000,
    'H2O': 178,
    'CH4': 2537.39,
    'CO2': 706.49,
} # mol/hours
PM = {
    'H2S': 34.1,
    'O2':  36.0,
    'SX':  32.1,
    'H2O': 18.0
}

Q_in = {n_in_key: n_in[n_in_key]*Rgas_m3_atm_K*T/P for n_in_key in n_in.keys()} # m3/hours
Q_sum = sum(Q_in.values())
tau_headspace =  V_gas/Q_sum
c_in = {n_in_key: n_in[n_in_key]/Q_sum for n_in_key in n_in.keys()} # mol/m3


def headspace_dyn_fanculo(x, t, Cin, *kwargs):
    global PM, 
    C_S, C_O, C_SX, C_H2O = x
    
    C_IN['H2S'] = Cin[0]
    C_IN['O2']  = Cin[1]

    w['H2S'] = c_out['H2S'] * PM['H2S'] * 1000 # g/m3
    w['O2']  = c_out['O2']  * PM['O2'] * 1000 # g/m3

    r_SOB = k * (w['H2S']**alfa) * (w['O2']**beta) # [1/h] - Reaction rate
    dC_S  =  (C_IN['H2S']-C_S)*Q/V_gas - r_SOB/PM['H2S'] # [mol/m3/h] - Sulfur concentration
    dC_O  =  (C_IN['O2']-C_O)*Q/V_gas - r_SOB/PM['O2']/RC # [mol/m3/h] - Oxygen concentration
    dC_S  = C_SX*Q/V_gas + r_SOB/PM['SX']# [mol/m3/h] - Sulfur concentration
    dC_H2O  = 0
    dCdt = [dC_S, dC_O, dC_SX, dC_H2O] 
    return dCdt

# Initial conditions
YOUT = odeint(headspace_dyn_fanculo, c_in.to_numpy, t, c_in, PM, V_gas, RC, alfa, beta, k)
plt.plot(t, YOUT[:,0], label='Sulfur')
plt.plot(t, YOUT[:,1], label='Oxygen')
plt.show()