# AMOCO_HN with modified identification - CONTAINS S2 AFFECTED BY SULFUR AND NOT

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*
from Influent import*

from functions_V3 import gompertz, growth_SRB, AD_OCS_Model, AMOCO_HN, f_deviations, deviations_check


# System definition
t_span = np.linspace(0,200,453) # time span

y_influent = f_deviations(t_span, T3.index.values, y_in_0) # Get the deviated influent values at each timestamp

# --------------------------------------------------------------------------------------------
# Ode Integration of AMOCO_HN: use to get X2 
# --------------------------------------------------------------------------------------------

y0= [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions established from SS

YOUT_pre = odeint(AMOCO_HN, y0, t_span, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values))
X2_pre = YOUT_pre[:,2]
print('************** AMOCOHN OK *******************')

# --------------------------------------------------------------------------------------------
# Ode Integration of the AD_OCS_Model: uses the previous X2 to get rho
# --------------------------------------------------------------------------------------------

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[5], SSTATE[6]] # initial conditions established from SS

YOUT = odeint(AD_OCS_Model, y0, t_span, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values, X2_pre, t_span, rates_df))

# Get results
XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = YOUT[:,5]              # [mmol/L] - VFA dissolved
S2_new = YOUT[:,6]          # [mmol/L] - VFA dissolved
C  = YOUT[:,7]              # [mmol/L] - Inorganic Carbon Dissolved

print('************** AD_OCS_Model OK *******************')

# Solver Output: from all the variables from the ones of the ODE
mu1 = np.zeros(len(t_span))
mu2 = np.zeros(len(t_span))
CO2 = np.zeros(len(t_span))
B   = np.zeros(len(t_span))
phi = np.zeros(len(t_span))
p_C = np.zeros(len(t_span))
q_C = np.zeros(len(t_span))
q_M = np.zeros(len(t_span))
pH  = np.zeros(len(t_span))

for x in range(len(t_span)):
    mu1[x] = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # [1/d]      - Specific Growth Rate for X1 (Monod)
    mu2[x] = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # [1/d]      - Specific Growth Rate for X2 (Haldane)
    CO2[x] = C[x] + S2[x] - Z[x]                                 # [mmol/L]   - Dissolved CO2
    B[x]   = Z[x] - S2[x]                                        # [mmol/L]   - Alkalinity
    phi[x] = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH) # [atm]      - CO2 Partial Pressure
    q_C[x] = kLa*(CO2[x] - KH*p_C[x])                            # [mmol/L/d] - CO2 Outlet Molar Flow
    q_M[x] = k[5]*mu2[x]*X2[x]                                   # [mmol/L/d] - CH4 Outlet Molar Flow
    pH[x]  = np.real(-np.log10(Kb*CO2[x]/B[x]))                  # [-]        - System pH

q_tot = q_C + q_M                                                   # [mmol/L/d] - Outlet global molar flow  

# Compute Mass Flows

q_C_W   = q_C*44/1000                                            # [g/L/d]    - CO2 Outlet mass flow of
q_M_W   = q_M*16/1000                                            # [g/L/d]    - CH4 Outlet mass flow  
q_tot_W = q_C_W + q_M_W                                          # [g/L/d]    - Outlet global mass flow  
x_M_W   = q_M_W/q_tot_W                                          # [-]        - CH4 Weight Fraction      

y_M   = np.divide(q_M,q_tot)                                     # [-]        - CH4 Mole fraction in gas phase
y_C   = np.divide(q_C,q_tot)                                     # [-]        - CO2 Mole fraction in gas phase

# Sulfur Influence Evaluation

Xs     = np.zeros(len(t_span))
lam    = np.zeros(len(t_span))
Ss     = np.zeros(len(t_span))
y_S    = np.zeros(len(t_span))
Ss_max = np.zeros(len(t_span))
Xs_max = np.zeros(len(t_span))
mu_srb = np.zeros(len(t_span))
growth_rate = np.zeros(len(t_span))

# mu_srb_loc = np.zeros([len(t_span),len(t_span)])
dXsdt = np.zeros([len(t_span),len(t_span)])
lam = 0 

for j in range(len(t_span)):
    # Iterate over each time step

    Ss_max[j] = frac_sulfur*y_influent[j,4]*1000/64*S2[j]/y_influent[0,4]
    Xs_max[j] = Y_srb/(1-Y_srb)*Ss_max[j] 
    
    
    # Gompertz function for microbial population and dissolved sulfur
    Xs[j]  = gompertz(t_span[j], Xs_max[j], mu_srb[j], 0)
    Ss[j]  = Xs[j]*(1-Y_srb)/(Y_srb)                                  # maximum sulfur concentration g/L
        
    mu_srb[j] = (- X2[0] + X2[j])/(t_span[j] - t_span[0]) 
    for snapshot in range(len(t_span)):
        
        mu_srb_loc = np.nan_to_num((- X2[0] + X2[snapshot])/(t_span[snapshot] - t_span[0]), nan=0, neginf=0)   
        # Gompertz function derivative for growth rate   
        dXsdt[j, snapshot] = growth_SRB(t_span[snapshot], Xs_max[j], mu_srb_loc, 0)
    
    # Retrieve the true growth rate
    growth_rate[j] = np.nanmax(dXsdt[j])

    # Sulfur G/L equilibrium
    y_S[j]   = (H_S*Ss[j])/1

# print(f'mu_srb: {mu_srb}')
# print(f'Xs_max: {Xs_max}')
# print(f'mu1,max: {mu1_max}; Ks1:  {KS1}; Cd1: {C_d[0]}')
# print(f'mu2,max: {mu2_max}; Ks2:  {KS2}; KI2: {KI2}; Cd2: {C_d[1]}')
# print(f"Mole fractions in the gas at the end: CH4: {y_M[-1]}, CO2 {y_C[-1]}")
# print(f"Mass fraction of methane in the gas at the end",float(x_M_W[-1]))

# Recalculation of y_i and flows
y_sum = y_M + y_C + y_S

y_M_new = y_M/y_sum
y_C_new = y_C/y_sum
y_S_new = y_S/y_sum

q_M_new = y_M_new*q_tot
q_C_new = y_C_new*q_tot
q_S_new = y_S_new*q_tot

q_C_W = q_C_new*MW_co/1000                                               # [g/L/d]    - CO2 Outlet mass flow of
q_M_W = q_M_new*MW_met/1000                                              # [g/L/d]    - CH4 Outlet mass flow  
q_S_W = q_S_new*MW_H2S/1000                                              # [g/L/d]    - CH4 Outlet mass flow  
q_W_tot = q_C_W+q_M_W+q_S_W

x_M_W = q_M_W/q_W_tot                                                    # [-]        - CH4 Weight Fraction
x_C_w = q_C_W/q_W_tot                                                    # [-]        - CO2 Weight Fraction
x_S_W = q_S_W/q_W_tot                                                    # [-]        - H2S Weight Fraction

plt.figure(100)
plt.subplot(4,1,1)
plt.plot(t_span, Xs, label="SRB")
plt.plot(t_span, Xs_max,'--', label="Max")
plt.legend()
plt.ylabel('Xs [g/L]')
plt.grid(True)

plt.subplot(4,1,2)
plt.grid(True)
plt.plot(t_span, Ss, label="Sulfur")
plt.plot(t_span, S2, label="S2")
plt.plot(t_span, Ss_max,'--', label="Max")
plt.xlabel('Time [d]')
plt.ylabel('S [mmol/L]')
plt.legend()

plt.subplot(4,1,3)
plt.plot(t_span, growth_rate)
plt.xlabel('Time [d]')
plt.ylabel('dX/dt [g/L/d]')

plt.grid(True)

plt.subplot(4,1,4)
plt.grid(True)
plt.plot(t_span, mu_srb)
plt.xlabel('Time [d]')
plt.ylabel('mu_max [g/L/d]')

# Create a figure to compare S2 and S2 new
plt.figure(103)
plt.subplot(3,1,1)
plt.plot(t_span, S2, label="S2")
plt.plot(t_span, S2_new,'--', label="new")
plt.ylabel('S2 [mmol/L]')
plt.legend()

plt.subplot(3,1,2)
plt.plot(t_span, (S2-S2_new))
plt.ylabel('delta S2')

delta = np.zeros(len(t_span))
for ii in range (len(t_span)):
    if ii > 0:
        delta[ii] = Xs[ii] - Xs[ii-1]
    else:
        delta[ii] = 0
    
plt.subplot(3,1,3)
plt.ylabel(' Xs growth [g/L/d]')
# plt.plot(t_span, mu_srb, label="mu_srb")
plt.plot(t_span, growth_rate, label="rho_srb")
plt.xlabel('Time [d]')
plt.legend()

plt.figure(104)
plt.subplot(2,1,1)
plt.plot(t_span, X2_pre, label="X2_pre")
plt.plot(t_span, X2, '--', label="X2")
plt.ylabel('X2 [g/L]')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t_span, X2_pre-X2)
plt.xlabel('Time [d]')
plt.ylabel('delta X2 [g/L]')

plt.figure(105)
# plot the growth rate
plt.subplot(2,1,1)
plt.plot(rates_df['time'], rates_df['rho'], label="growth rate")
plt.plot(t_span, growth_rate, '--', label="growth rate new")
plt.legend()
print(rates_df['rho'].head())
plt.show()