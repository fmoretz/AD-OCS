
''' AD_OCS Model with sulfur addition and oxygen implementation
    Version 3.2: Headspace divided in two regions: 
        1. intephase, where H2S is in equilibrium with the liquid;
        2. gas phase, where Oxygen is present and reacts with H2S.'''

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*
from Influent import*

from functions_V3 import gompertz, growth_SRB, AD_OCS_Model, AMOCO_HN, f_deviations, deviations_check

# System definition
t_span = np.linspace(10,30,40) # time span

y_influent = f_deviations(t_span, T3.index.values, y_in_0) # Get the deviated influent values at each timestamp

# --------------------------------------------------------------------------------------------
# Ode Integration of AMOCO_HN: use to get X2 
# --------------------------------------------------------------------------------------------

y0= [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

YOUT_pre = odeint(AMOCO_HN, y0, t_span, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values))
X2_pre = YOUT_pre[:,2]
print('************** AMOCOHN OK *******************')

# --------------------------------------------------------------------------------------------
# Ode Integration of the AD_OCS_Model: uses the previous X2 to get rho
# --------------------------------------------------------------------------------------------

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

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

q_tot = q_C + q_M                                                # [mmol/L/d] - Outlet global molar flow  

# Compute Mass Flows

q_C_W   = q_C*44/1000                                            # [g/L/d]    - CO2 Outlet mass flow of
q_M_W   = q_M*16/1000                                            # [g/L/d]    - CH4 Outlet mass flow  
q_tot_W = q_C_W + q_M_W                                          # [g/L/d]    - Outlet global mass flow  
x_M_W   = q_M_W/q_tot_W                                          # [-]        - CH4 Weight Fraction      

y_M   = np.divide(q_M,q_tot)                                     # [-]        - CH4 Mole fraction in gas phase
y_C   = np.divide(q_C,q_tot)                                     # [-]        - CO2 Mole fraction in gas phase

# Sulfur and Oxygen Influence Evaluation
# Oxygen test


Xs     = np.zeros(len(t_span))                                   # [g/L]      - Sulfate Reducing Bacteria
Ss     = np.zeros(len(t_span))                                   # [mmol/L]      - Sulfur dissolved
y_S_int = np.zeros(len(t_span))                                   # [-]        - Sulfur Mole fraction in gas phase
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
    Xs[j]  = gompertz(t_span[j], Xs_max[j], mu_srb[j], lam)           # [g/L]      - Sulfate Reducing Bacteria - Gompertz
    Ss[j]  = Xs[j]*(1-Y_srb)/(Y_srb)                                  # [g/L]      - Sulfur dissolved concentration
        
    mu_srb[j] = (- X2[0] + X2[j])/(t_span[j] - t_span[0])             # [g/L/d]    - Gompertz parameter for SRB growth
    for snapshot in range(len(t_span)):
        
        mu_srb_loc = np.nan_to_num((- X2[0] + X2[snapshot])/(t_span[snapshot] - t_span[0]), nan=0, neginf=0)    # [g/L/d]    - Gompertz parameter for SRB growth, local  
        dXsdt[j, snapshot] = growth_SRB(t_span[snapshot], Xs_max[j], mu_srb_loc, lam)                     # [g/L/d]    - SRB growth rate matrix, local (Gompertz Derivative)
    
        growth_rate[j] = np.nanmax(dXsdt[j])                                  # [g/L/d]    - Get the growth rate of SRB at each time step as the maximum of the possible rates

    y_S_int[j]   = (H_S*Ss[j])/1                                    # [-]        - Sulfur Mole fraction in gas phase (G/L equilibrium)


# Recalculation of y_i and flows by normalization - Version 1
# y_sum = y_M + y_C + y_S_int

# y_M_new = y_M/y_sum                                               # [-]        - Untreated CH4 Mole fraction in gas phase
# y_C_new = y_C/y_sum                                               # [-]        - Untreated CO2 Mole fraction in gas phase
# y_S_int = y_S_int/y_sum                                               # [-]        - Untreated Sulfur Mole fraction in gas phase

# q_M_new = q_tot*y_M_new                                           # [mmol/L/d] - Untreated CH4 Outlet Molar Flow
# q_C_new = q_tot*y_C_new                                           # [mmol/L/d] - Untreated CO2 Outlet Molar Flow
# q_S_new = q_tot*y_S_int                                           # [mmol/L/d] - Untreated Sulfur Outlet Molar Flow

# n_C = V_liq*q_C_new                                               # [mmol/d]   - CO2 Outlet Molar Flow
# n_M = V_liq*q_M_new                                               # [mmol/d]   - CH4 Outlet Molar Flow
# n_S = V_liq*q_S_new                                               # [mmol/d]   - Sulfur Outlet Molar Flow

# n_out = n_C + n_M + n_S                                           # [mmol/d]   - Total Outlet Molar Flow

# Recalculation of y_i and flows by normalization - Version 2 (Should we switch to mol instead of mmol?)
n_C = V_liq*q_C                                              # [mmol/d]   - CO2 Outlet Molar Flow
n_M = V_liq*q_M                                              # [mmol/d]   - CH4 Outlet Molar Flow

n_out = (n_C + n_M)/(1-y_S_int)                                   # [mmol/d]   - Total Outlet Molar Flow

n_S_int = n_out*y_S_int                                           # [mmol/d]   - Sulfur Outlet Molar Flow

coeff_S = -1
coeff_O2 = -0.5
alfa = -0.8
beta = 1
k_SOB = 0.1

# Vectors allocation
y_S = np.zeros(len(t_span))                                   # [-]        - Sulfur Mole fraction in gas phase
y_O2 = np.zeros(len(t_span))                                   # [-]        - Oxygen Mole fraction in gas phase
n_S_out = np.zeros(len(t_span))                                   # [mmol/d]   - Sulfur Outlet Molar Flow
n_O2_out = np.zeros(len(t_span))                                   # [mmol/d]   - Oxygen Outlet Molar Flow
n_M_out = np.zeros(len(t_span))                                   # [mmol/d]   - Methane Outlet Molar Flow
n_C_out = np.zeros(len(t_span))                                   # [mmol/d]   - Carbon Dioxide Outlet Molar Flow

for i in range(len(t_span)):
    n_O2_in = 800 # [mmol/d] - Influent O2 flowrate PROVA!!!!
    
    def Headspace_reactions(y):
        global coeff_S, coeff_O2, V_gas, k_SOB, alfa, beta
        y_S_loc = y[0]
        y_O_loc = y[1]

        n_S = n_out[i]*y_S_loc
        n_O2 = n_out[i]*y_O_loc

        r_SOB = k_SOB*y_O_loc**alfa*y_S_loc**beta                          # [mmol/L/d] - Reaction rate SOB
        BM_S  = n_S_int[i] - n_S     + coeff_S*V_gas*r_SOB                               # [mmol/d]   - Sulfur Balance
        BM_O2 = n_O2_in    - n_O2    + coeff_O2*V_gas*r_SOB                                  # [mmol/d]   - Oxygen Balance
        return [BM_S, BM_O2]

    y0 = [0.02, 0.015]
    [y_S[i], y_O2[i]] = fsolve(Headspace_reactions, y0)
    print(y_S[i], y_O2[i])

y_M = n_M/n_out
y_C = n_C/n_out

sumy = y_S + y_O2 + y_M + y_C
y_S_out = y_S/sumy
y_O2_out = y_O2/sumy
y_M_out = y_M/sumy
y_C_out = y_C/sumy

n_S_out = n_out*y_S_out
n_O2_out = n_out*y_O2_out
n_M_out = n_out*y_M_out
n_C_out = n_out*y_C_out

color_S = 'goldenrod'
color_O2 = 'skyblue'
color_M = 'violet'
color_C = 'orange'


plt.figure(1)
plt.title('Treatment Effectiveness')
plt.plot(t_span, y_S, color= color_S, label='H2S')
plt.plot(t_span, y_S_int, '--', color= color_S, label='H2S Int')
plt.plot(t_span, y_S_out, ':', color= color_S, label='H2S Norm.')
plt.plot(t_span, y_O2, color = color_O2,label='O2')
plt.plot(t_span, y_O2_out, ':', color = color_O2,label='O2 Norm.')
plt.legend()

plt.figure(2)
plt.suptitle('Out flows and moalr fractions with normalization effect (dotted lines)')

plt.subplot(2,1,1)
plt.title('Outlet Gas Molar Flow')
plt.plot(t_span, n_S_out, color= color_S, label='H2S')
plt.plot(t_span, n_O2_out, color = color_O2,label='O2')
plt.plot(t_span, n_M_out, color = color_M,label='CH4')
plt.plot(t_span, n_C_out, color = color_C,label='CO2')
plt.plot(t_span, n_C, ':', color = color_C)
plt.plot(t_span, n_M, ':', color = color_M)
plt.ylabel('Molar Flow [mmol/d]')
plt.legend()

plt.subplot(2,1,2)
plt.title('Gas Molar Fractions')
plt.plot(t_span, y_S_out, color= color_S, label='H2S')
plt.plot(t_span, y_O2_out, color = color_O2,label='O2')
plt.plot(t_span, y_M_out, color = color_M,label='CH4')
plt.plot(t_span, y_C_out, color = color_C,label='CO2')
plt.plot(t_span, y_S, ':', color= color_S)
plt.plot(t_span, y_O2, ':', color = color_O2)
plt.plot(t_span, y_M, ':', color = color_M)
plt.plot(t_span, y_C, ':', color = color_C)
plt.ylabel('Molar Faction [-]')
plt.legend()

plt.show()

