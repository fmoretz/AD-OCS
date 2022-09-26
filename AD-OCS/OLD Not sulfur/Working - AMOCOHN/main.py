# AMOCO_HN with modified identification

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*

# SYSTEM

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]]

def f_Model_Deviations_Simple(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,y_in,k,kd,N_bac,N_S1):

    XT, X1, X2, Z, S1, S2, C = x
    print(t)
    if t < 20 or t > 100:
        S1in = y_in[0]      # [gCOD/L]
        S2in = y_in[1]      # [mmol/L]
        Cin  = y_in[2]      # [mmol/L]
        Zin  = y_in[3]      # [mmol/L]
        XTin = y_in[4]      # [gCOD/L]
    else:
        S1in = y_in[0]      # [gCOD/L]   
        S2in = y_in[1]      # [mmol/L]
        Cin  = y_in[2]      # [mmol/L]
        Zin  = y_in[3]      # [mmol/L]
        XTin = 1.2*y_in[4]  # [gCOD/L]

    mu1 = mu_max[0]*(S1/(S1+Ks[0]))                                                                  # Monod
    mu2 = mu_max[1]*(S2/(S2+Ks[1]+S2**2/KI2))                                                        # Haldane

    qM  = k[5]*mu2*X2                                                                                # [mmol/L/day] - Methane molar flow 
    CO2 = C + S2 - Z                                                                                 # [mmol/L]     - CO2 Dissolved
    phi = CO2 + KH*Pt + qM/kLa
    Pc  = (phi - (phi**2- 4*KH*Pt*CO2)**0.5)/(2*KH)                                                  # [atm] - Partial pressure CO2
    qC  = kLa*(CO2 - KH*Pc)                                                                          # [mmol/L/d] - Carbon Molar Flow
                                             

    dXT = D*(XTin - XT) - k[6]*XT                                                                    # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                                                                  # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                                                                  # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2  # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                      # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2                                                      # Evolution of VFA
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC                                                 # Evolution of inorganic carbon

    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dC]

    return dxdt

tspan = np.linspace(1,200,100)
YOUT = odeint(f_Model_Deviations_Simple,y0,tspan,args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, y_in, k, kd, N_bac, N_S1))

XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = np.array(YOUT[:,5])              # [mmol/L] - VFA dissolved
C  = YOUT[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

# Solver Output
mu1 = np.empty(len(XT))
mu2 = np.empty(len(XT))
CO2 = np.empty(len(XT))
B   = np.empty(len(XT))
phi = np.empty(len(XT))
p_C = np.empty(len(XT))
q_C = np.empty(len(XT))
q_M = np.empty(len(XT))
pH  = np.empty(len(XT))



for x in range(len(XT)):
    mu1[x] = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # [1/d]      - Specific Growth Rate for X1 (Monod)
    mu2[x] = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # [1/d]      - Specific Growth Rate for X2 (Haldane)
    CO2[x] = C[x] + S2[x] - Z[x]                                 # [mmol/L]   - Dissolved CO2
    B[x]   = Z[x] - S2[x]                                        # [mmol/L]   - Alkalinity
    phi[x] = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH) # [atm]      - CO2 Partial Pressure
    q_C[x] = kLa*(CO2[x] - KH*p_C[x])                            # [mmol/L/d] - CO2 Outlet Molar Flow
    q_M[x] = k[5]*mu2[x]*X2[x]                                   # [mmol/L/d] - CH4 Outlet Molar Flow
    pH[x]  = np.real(-np.log10(Kb*CO2[x]/B[x]))                  # [-]        - System pH

q_tot = q_C+q_M                                                  # [mmol/L/d] - Outlet global molar flow  

q_C_W = q_C*44/1000                                              # [g/L/d]    - CO2 Outlet mass flow of
q_M_W = q_M*16/1000                                              # [g/L/d]    - CH4 Outlet mass flow  
q_tot_W = q_C_W + q_M_W                                          # [g/L/d]    - Outlet global mass flow  
x_M_W   = q_M_W/q_tot_W                                          # [-]        - CH4 Weight Fraction      

y_M   = np.divide(q_M,q_tot)                                     # [-]        - CH4 Mole fraction in gas phase
y_C   = np.divide(q_C,q_tot)                                     # [-]        - CO2 Mole fraction in gas phase

print(f'mu1,max: {mu1_max}; Ks1:  {KS1}; Cd1: {C_d[0]}')
print(f'mu2,max: {mu2_max}; Ks2:  {KS2}; KI2: {KI2}; Cd2: {C_d[1]}')
print(f"Mole fractions in the gas at the end: CH4: {y_M[-1]}, CO2 {y_C[-1]}")
print(f"Mass fraction of methane in the gas at the end",float(x_M_W[-1]))
plt.close("all")
