''' AD_OCS Model with sulfur addition'''


import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*
from Influent import*

from functions_V3 import gompertz, growth_SRB, AD_OCS_Model, AMOCO_HN, f_deviations, deviations_check


# System definition
t_span = np.linspace(0,200,300) # time span [d]

y_influent = f_deviations(t_span, T3.index.values, y_in_0[:5]) # Get the deviated influent values at each timestamp

# --------------------------------------------------------------------------------------------
# Ode Integration of AMOCO_HN: used to get X2 
# --------------------------------------------------------------------------------------------

y0= [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions established from SS

YOUT_pre = odeint(AMOCO_HN, y0, t_span, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values))
X2_pre = YOUT_pre[:,2]
print('************** AMOCOHN OK *******************')

# --------------------------------------------------------------------------------------------
# Ode Integration of the AD_OCS_Model: uses the previous X2 to get rho
# --------------------------------------------------------------------------------------------

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[5], SSTATE[6]] # initial conditions established from SS

YOUT = odeint(AD_OCS_Model, y0, t_span, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values, X2_pre, t_span))

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

# Sulfur and Oxygen Influence Evaluation

q_O2   = np.zeros(len(t_span))                                   # [mmol/L/d] - O2 Outlet Molar Flow
r_O2   = np.zeros(len(t_span))                                   # [mmol/L/d] - Reaction rate O2-H2S

Xs     = np.zeros(len(t_span))                                   # [g/L]      - Sulfate Reducing Bacteria
Ss     = np.zeros(len(t_span))                                   # [mmol/L]      - Sulfur dissolved
y_S    = np.zeros(len(t_span))                                   # [-]        - Sulfur Mole fraction in gas phase
Ss_max = np.zeros(len(t_span))                                   # [mmol/L]      - Maximum Sulfur dissolved concentration
Xs_max = np.zeros(len(t_span))                                   # [g/L]      - Maximum Sulfate Reducing Bacteria concentration (Gompertz Asymptote)
mu_srb = np.zeros(len(t_span))                                   # [g/L/d]    - Gompertz parameter for SRB growth
growth_rate = np.zeros(len(t_span))                              # [g/L/d]    - Gompertz derivative: states the growth rate of SRB

dXsdt = np.zeros([len(t_span),len(t_span)])                     # [g/L/d]    - SRB growth rate matrix, preallocation
lam = 0                                                         # [-]        - Lag Time as defined by Gompertz

for j in range(len(t_span)):
    # Iterate over each time step

    Ss_max[j] = frac_sulfur*y_influent[j,4]*1000/64*S2[j]/y_influent[0,4]
    Xs_max[j] = Y_srb/(1-Y_srb)*Ss_max[j]     
    
    # Gompertz function for microbial population and dissolved sulfur
    mu_srb[j] = (- X2[0] + X2[j])/max(1e-14,t_span[j]-t_span[0])            # [g/L/d]    - Gompertz parameter for SRB growth
    Xs[j]  = gompertz(t_span[j], Xs_max[j], mu_srb[j], lam)           # [g/L]      - Sulfate Reducing Bacteria - Gompertz
    Ss[j]  = Xs[j]*(1-Y_srb)/(Y_srb)                                  # [g/L]      - Sulfur dissolved concentration
        
    for snapshot in range(len(t_span)): # this ruotine is used just to show the growth rate of SRB which intervenes in the ODE
        
        mu_srb_loc = np.nan_to_num((- X2[0] + X2[snapshot])/max(1e-14,(t_span[snapshot] - t_span[0])), nan=0, neginf=0)    # [g/L/d]    - Gompertz parameter for SRB growth, local  
        dXsdt[j, snapshot] = growth_SRB(t_span[snapshot], Xs_max[j], mu_srb_loc, lam)                           # [g/L/d]    - SRB growth rate matrix, local (Gompertz Derivative)
    
        growth_rate[j] = np.nanmax(dXsdt[j])                                  # [g/L/d]    - Get the growth rate of SRB at each time step as the maximum of the possible rates

    y_S[j]   = (KH_S*Ss[j])/P_dig                                    # [-]        - Sulfur Mole fraction in gas phase (G/L equilibrium)

y_sum = y_M + y_C + y_S

y_M_new = y_M/y_sum                                               # [-]        - Untreated CH4 Mole fraction in gas phase
y_C_new = y_C/y_sum                                               # [-]        - Untreated CO2 Mole fraction in gas phase
y_S_new = y_S/y_sum                                               # [-]        - Untreated Sulfur Mole fraction in gas phase

q_M_new = y_M_new*q_tot
q_C_new = y_C_new*q_tot
q_S_new = y_S_new*q_tot

# plt.figure()
# plt.subplot(2,1,1)
# plt.plot(t_span, q_M_new, label='CH4')
# plt.plot(t_span, q_C_new, label='CO2')
# plt.plot(t_span, q_S_new, label='Sulfur')
# plt.legend()
# plt.xlabel('Time [d]')
# plt.ylabel('Molar Flow [mmol/L/d]')
# plt.subplot(2,1,2)
# plt.plot(t_span, y_M_new, label='CH4')
# plt.plot(t_span, y_C_new, label='CO2')
# plt.plot(t_span, y_S_new, label='Sulfur')
# plt.legend()
# plt.xlabel('Time [d]')
# plt.ylabel('Mole Fraction [-]')

# plt.figure()
# plt.suptitle('Sulfur Kinetics')
# plt.subplot(4,1,1)
# plt.plot(t_span, Xs, label='SRB')
# plt.ylabel(ylabel='Concentration [g/L]')
# plt.subplot(4,1,2)
# plt.plot(t_span, Ss, label='Sulfur')
# plt.plot(t_span, Ss_max, '--', label='Sulfur Max')
# plt.legend()
# plt.ylabel('Concentration [mmol/L]')
# plt.subplot(4,1,3)
# plt.plot(t_span, growth_rate, label='Growth Rate')
# plt.ylabel('Growth Rate [g/L/d]')
# plt.legend()
# plt.subplot(4,1,4)
# plt.plot(t_span, mu_srb, label='Gompertz Parameter')
# plt.ylabel('Gompertz Parameter [g/L/d]')
# plt.legend()


# plt.figure()
# plt.suptitle('Sulfur Kinetics effect on S2')
# plt.subplot(2,1,1)
# plt.plot(t_span, S2, label='S2')
# plt.plot(t_span, S2_new, label='S2_new')
# plt.legend()
# plt.ylabel('Concentration [mmol/L]')
# plt.subplot(2,1,2)
# plt.plot(t_span, S2-S2_new, label='S2-S2_new')
# plt.ylabel('Difference [mmol/L]')
# # plt.show()
