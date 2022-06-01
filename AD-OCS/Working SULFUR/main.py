# AMOCO_HN with modified identification

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from deviationscopy2 import*
from SS_Algebraic import*

# SYSTEM

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]]

def gompertz(x,a,b,c):
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))
def growth_SRB(x,a,b,c):
    return b*np.exp((np.exp(1)*b/a*(c-x))-np.exp(np.exp(1)*b/a*(c-x)+1)+2)

def f_Model_Deviations_Simple(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,k,kd,N_bac,N_S1,X2_0,t_0,y_in):

    XT, X1, X2, Z, S1, S2, C = x
    
    X2_0 = 0.0
    XT_in_0 =y_in_0[4]
    y_in = deviation_check(t,y_in)
    
    S1in = y_in[0]      # [gCOD/L]
    S2in = y_in[1]      # [mmol/L]
    Cin  = y_in[2]      # [mmol/L]
    Zin  = y_in[3]      # [mmol/L]
    XTin = y_in[4]      # [gCOD/L]
   
    mu1 = mu_max[0]*(S1/(S1+Ks[0]))                                                                  # Monod
    mu2 = mu_max[1]*(S2/(S2+Ks[1]+S2**2/KI2))                                                        # Haldane

    qM  = k[5]*mu2*X2                                                                                # [mmol/L/day] - Methane molar flow 
    CO2 = C + S2 - Z                                                                                 # [mmol/L]     - CO2 Dissolved
    phi = CO2 + KH*Pt + qM/kLa
    Pc  = (phi - (phi**2- 4*KH*Pt*CO2)**0.5)/(2*KH)                                                  # [atm] - Partial pressure CO2
    qC  = kLa*(CO2 - KH*Pc)                                                                          # [mmol/L/d] - Carbon Molar Flow

    Ss_max = 0.02*XT_in*1000/64*S2/XT_in_0
    Xs_max = Y_srb/(1-Y_srb)*Ss_max              # maximum sulfur concentration g/L
    mu_max_srb = np.nan_to_num((- X2 + X2)/(t_span[i] - t_span[0]), nan=0, neginf=0)                                         
    rho_srb = growth_SRB(t, Xs_max, mu_max_srb, 0)

    dXT = D*(XTin - XT) - k[6]*XT                                                                    # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                                                                  # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                                                                  # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2  # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                      # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2 - rho_srb/Y_srb                                                 # Evolution of VFA
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC                                             # Evolution of inorganic carbon

    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dC]

    return dxdt

YOUT = odeint(f_Model_Deviations_Simple,y0,t_span,args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0))

XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = YOUT[:,5]              # [mmol/L] - VFA dissolved
C  = YOUT[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

# Solver Output
mu1 = np.zeros(len(XT))
mu2 = np.zeros(len(XT))
CO2 = np.zeros(len(XT))
B   = np.zeros(len(XT))
phi = np.zeros(len(XT))
p_C = np.zeros(len(XT))
q_C = np.zeros(len(XT))
q_M = np.zeros(len(XT))
pH  = np.zeros(len(XT))

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

# Calculations for the gas flows and mass
q_tot = q_C+q_M                                                  # [mmol/L/d] - Outlet global molar flow  

q_C_W = q_C*MW_co/1000                                              # [g/L/d]    - CO2 Outlet mass flow of
q_M_W = q_M*MW_met/1000                                              # [g/L/d]    - CH4 Outlet mass flow  
q_tot_W = q_C_W + q_M_W                                          # [g/L/d]    - Outlet global mass flow  
x_M_W   = q_M_W/q_tot_W                                          # [-]        - CH4 Weight Fraction      

y_M   = np.divide(q_M,q_tot)                                     # [-]        - CH4 Mole fraction in gas phase
y_C   = np.divide(q_C,q_tot)                                     # [-]        - CO2 Mole fraction in gas phase

# Sulfur 

Xs = np.zeros(len(XT))
lam = np.zeros(len(XT))
Ss  = np.zeros(len(XT))
y_S  = np.zeros(len(XT))
Ss_max= np.zeros(len(XT))
Xs_max = np.zeros(len(XT))
mu_max_srb = np.zeros(len(XT))

for i in range(len(XT)):
    # Species differences
    
    Ss_max[i] = frac_sulfur*y_influent[i,4]*1000/64*S2[i]/y_influent[0,4]
    Xs_max[i] = Y_srb/(1-Y_srb)*Ss_max[i]   # maximum sulfur concentration g/L
    mu_max_srb[i] = np.nan_to_num((- X2[0] + X2[i])/(t_span[i] - t_span[0]), nan=0, neginf=0)
    lam[i]    = 0
    
    # Gompertz function
    Xs[i]    = gompertz(t_span[i], Xs_max[i], mu_max_srb[i], lam[i])
    Ss[i]    = (1-Y_srb)/(Y_srb)*Xs[i]
    y_S[i]   = (H_S*Ss[i])/1
    
print(f'mu1,max: {mu1_max}; Ks1:  {KS1}; Cd1: {C_d[0]}')
print(f'mu2,max: {mu2_max}; Ks2:  {KS2}; KI2: {KI2}; Cd2: {C_d[1]}')
print(f"Mole fractions in the gas at the end: CH4: {y_M[-1]}, CO2 {y_C[-1]}")
print(f"Mass fraction of methane in the gas at the end",float(x_M_W[-1]))

# Recalculation of and flows and molar fractions
q_S = np.zeros(len(XT))
y = np.transpose(np.vstack((y_M,y_C,y_S)))
q = np.transpose(np.vstack((q_M,q_C,q_S)))
q_W = np.transpose(np.vstack((q_M_W,q_C_W,q_S)))
y_new   = np.zeros([len(y),3])
q_new   = np.zeros([len(y),3])
q_W_new  = np.zeros([len(y),3])
q_W_tot_new   = np.zeros(len(y))
x_W_new = np.zeros([len(y),3])


for i in range(len(XT)):   
    y_sum = sum(y[i])
    y_new[i] = y[i]/y_sum
    q_new[i] = y_new[i]*q_tot[i]
    
    q_W_new[i] = q_new[i]*MW_gas/1000

    q_W_tot_new[i] = sum(q_W_new[i])   
    x_W_new[i] = q_W_new[i]/q_W_tot_new[i]

plt.figure(100)
plt.plot(t_span, y, label=("CH4","CO2","H2S"))
plt.plot(t_span, y_new, '--', label=("new CH4","new CO2","new H2S"))
plt.ylabel('Molar frac')
plt.legend()
plt.grid(True)

plt.figure(101)
plt.plot(t_span, q, label=("CH4","CO2","H2S"))
plt.plot(t_span, q_new, '--', label=("new CH4","new CO2","new H2S"))
plt.ylabel('Molar flows [mmol/L/d]')
plt.legend()
plt.grid(True)

plt.figure(102)
plt.plot(t_span, x_W_new, label=("CH4","CO2","H2S"))
plt.ylabel('Mass Fractions') 
plt.legend()
plt.grid(True)

plt.figure(103)
plt.plot(t_span, q_W_new, label=("CH4","CO2","H2S"))
plt.ylabel('Mass flows [g/L/d]')
plt.legend()
plt.grid(True)

print(f'flows: {q_W_new[-1]}')
print(f'fractions: {x_W_new[-1]}')
# plt.show()