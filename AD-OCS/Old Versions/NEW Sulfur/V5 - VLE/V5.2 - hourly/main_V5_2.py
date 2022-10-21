''' V5.0 AD_OCS Model with sulfur addition and VLE equilibrium
    V5.1 - Level: Level control of the digester
    V5.2 - hourly: hourly control of the digester'''

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*
from Influent import*
from GLConstants import*

from functions_V5 import gompertz, growth_SRB, AD_OCS_Model, AMOCO_HN, f_deviations, deviations_check, level_t
from GL_Equilibrium import f_RR_equilibrium


# System definition
d_start = 0    # [h] - Start time
d_end   = 60  # [h] - End time
hours   = 1    # [h] - Discretization time
n_times = int((d_end-d_start)*hours)+1 # Number of time steps

print('***Intervals of {hours} hours. \n {n_times} time steps***'.format(hours=hours, n_times=n_times))

t_span = np.linspace(d_start,d_end,n_times) # time span
t_span_d = t_span/24 # time span in days

y_influent_changes = f_deviations(t_span_d, T3.index.values, y_in_0) # Get the deviated influent values at each timestamp
y_influent = y_influent_changes[:,0:5] # Remove the Q_in values

Q_in = y_influent_changes[:,5] # Remove the y_in values

# --------------------------------------------------------------------------------------------
# Ode Integration of AMOCO_HN: use to get X2 
# --------------------------------------------------------------------------------------------

y0= [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

YOUT_pre = odeint(AMOCO_HN, y0, t_span_d, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values))
X2_pre = YOUT_pre[:,2]
print('************** AMOCOHN OK *******************')

# --------------------------------------------------------------------------------------------
# Ode Integration of the AD_OCS_Model: uses the previous X2 to get rho
# --------------------------------------------------------------------------------------------

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

YOUT = odeint(AD_OCS_Model, y0, t_span_d, args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values, X2_pre, t_span))

# Get results
XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = YOUT[:,5]              # [mmol/L] - VFA dissolved
C  = YOUT[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

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

# Sulfur and Oxygen Influence Evaluation

Xs     = np.zeros(len(t_span))                                   # [g/L]      - Sulfate Reducing Bacteria
Ss     = np.zeros(len(t_span))                                   # [mmol/L]      - Sulfur dissolved
y_S = np.zeros(len(t_span))                                      # [-]        - Sulfur Mole fraction in gas phase
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
    mu_srb[j] = (- X2[0] + X2[j])/max(1e-14,t_span[j]-t_span[0])      # [g/L/d]    - Gompertz parameter for SRB growth
    Xs[j]  = gompertz(t_span[j], Xs_max[j], mu_srb[j], lam)           # [g/L]      - Sulfate Reducing Bacteria - Gompertz
    Ss[j]  = Xs[j]*(1-Y_srb)/(Y_srb)                                  # [g/L]      - Sulfur dissolved concentration
        
    for snapshot in range(len(t_span)):        
        mu_srb_loc = np.nan_to_num((- X2[0] + X2[snapshot])/(max(1e-14,t_span[snapshot] - t_span[0])), nan=0, neginf=0)    # [g/L/d]    - Gompertz parameter for SRB growth, local  
        dXsdt[j, snapshot] = growth_SRB(t_span[snapshot], Xs_max[j], mu_srb_loc, lam)                     # [g/L/d]    - SRB growth rate matrix, local (Gompertz Derivative)
    
        growth_rate[j] = np.nanmax(dXsdt[j])                  # [g/L/d]    - Get the growth rate of SRB at each time step as the maximum of the possible rates

    y_S[j]   = (KH_S*Ss[j])/P_dig                             # [-]        - Sulfur Mole fraction in gas phase (G/L equilibrium)
q_S = y_S*q_tot 
plt.figure()                                              # [mmol/L/d] - Sulfur Outlet Specific Molar Flow
plt.plot(t_span, q_S, label='Sulfur')
plt.plot(t_span, q_C, label='Carbon')
plt.plot(t_span, q_M, label='Methane')
plt.xlabel('Time [d]')
plt.show()
### Assess the liquid volume dynamics
h = np.zeros(len(t_span))
print(T3.index.values)

h0 = h_SS
t_change = t_span[0]
index = 0
V_liq = np.zeros(len(t_span))                                 # [m3] - Liquid Volume

for i in range(len(t_span)):
    t = t_span[i]    
    if t < T3.index.values[index]*24:
        pass

    elif t >= T3.index.values[-1]*24:
        t_change = T3.index.values[index]*24 
        h0 = h[len(t_span[t_span < t_change])-1]
        
    else:    
        t_change = T3.index.values[index]*24       
        h0 = h[i-1]
        index = min(index + 1, len(T3.index.values)-1)  
                
    h[i] = level_t(t, D, Q_in[i], SR, h0, t_change)
    V_liq[i] = np.pi*Dr**2*h[i]/4                             # [m3] - Liquid Volum
    if h[i] > hmax:
        print('Kittemmuort lo stai riempiendo troppo')
        input('Press Enter to continue')
    elif h[i] < hmin:
        print('Kittemmuort lo stai svuotando troppo')
        input('Press Enter to continue')
        
#### Evaluation of the G/L equilibrium effects on the system ####

n_species = 4                                                 # [-]       - Number of species in the system

F_i = np.zeros([len(t_span),n_species])                       # [-] - Outlet Molar Fraction - In of the flash
F_M = np.zeros(len(t_span))                                   # [-] - Methane Flow
F_C = np.zeros(len(t_span))                                   # [-] - Carbon Flow
F_S = np.zeros(len(t_span))                                   # [-] - Sulfur Flow
F_W = np.zeros(len(t_span))                                   # [-] - Water Flow

for i in range(len(t_span)):  
    F_W[i] = water_percentage*Q_in[i]*rho_water/18*1000/24         # [mol/h]    - Water Flow 
    F_M[i] = q_M[i]*V_liq[i]/24                                    # [mol/h]   - Methane Flow 
    F_C[i] = q_C[i]*V_liq[i]/24                                    # [mol/h]   - CO2 Flow
    F_S[i] = q_S[i]*V_liq[i]/24                                    # [mol/h]   - Sulfur Flow

    F_i[i] = np.array([F_M[i], F_C[i], F_S[i], F_W[i]])         # [mol/h] - Outlet Molar Flow - In of the flash

z_i = np.zeros([len(t_span),n_species])                         # [-] - Outlet Molar Fraction - In of the flash
x_i = np.zeros([len(t_span),n_species])                         # [-] - Outlet Molar Fraction - Liquid phase of the flash
y_i = np.zeros([len(t_span),n_species])                         # [-] - Outlet Molar Fraction - Vapor phase of the flash

K = np.zeros([len(t_span),n_species])                         # [-] - Equilibrium Constant
num = np.zeros([len(t_span),n_species])                       # [-] - Numerator of the equilibrium equation
den = np.zeros([len(t_span),n_species])                       # [-] - Denominator of the equilibrium equation

N_L = np.zeros([len(t_span),n_species])                       # [-] - Liquid Phase Mole Fraction
N_V = np.zeros([len(t_span),n_species])                       # [-] - Vapor Phase Mole Fraction

alpha = np.zeros(len(t_span))                                 # [-]     - 
law = ('H', 'H', 'H', 'R')                                    # [-]     - Define the law to be used for the equilibrium calculation (R for Raoult's Law, H for Henry's Law)
for t in range(len(t_span)):  
    F_in = sum(F_i[t])                                        # [mol/h] - Total Inlet Molar Flow
    for i in range(n_species):        
        if law[i] == 'R':
            P_sp = P_sat[i]
        else:
            P_sp = H_atm[i]             
        K[t,i] = P_sp/P_dig
        z_i[t,i] = F_i[t,i]/F_in
            
    guess = sum(z_i[t,0:2])
    alpha[t] = fsolve(f_RR_equilibrium, guess , args=(z_i[t,:], law, P_sat, H_atm, P_dig)) # [-] - Solve the G/L equilibrium for the alpha factor
    N_V_tot = F_in*alpha[t]                                       # [mol/h] - Vapor Molar Flow
    N_L_tot = F_in - N_V_tot                                      # [mol/h] - Liquid Molar Flow
    for i in range(n_species):
        x_i[t,i] = z_i[t,i]/(1+alpha[t]*(K[t,i]-1))
        y_i[t,i] = K[t,i]*x_i[t,i]
        N_L[t,i] = N_L_tot*x_i[t,i]
        N_V[t,i] = N_V_tot*y_i[t,i]

entrance = sum(F_i[0,0:2])
out = sum(N_V[0,0:2])
delta = entrance - out
print('Entrance: ', entrance)
print('Out: ', out)
print('Delta: ', delta)
print('Relative Error: ', delta/entrance*100,'%')

print(N_V[0])
print(N_L[0])

print('Alpha',alpha[0])
print('Liquid:', x_i[0,:])
print('Vapour:',y_i[0,:])
plt.figure()
plt.subplot(2,1,1)
plt.plot(t_span, x_i[:,0], label='Methane')
plt.plot(t_span, x_i[:,1], label='CO2')
plt.plot(t_span, x_i[:,2], label='Sulfur')
plt.plot(t_span, x_i[:,3], label='Water')
plt.legend()
plt.xlabel('Time [h]')
plt.ylabel('Molar Fraction [-]')
plt.title('Liquid Phase')
plt.subplot(2,1,2)
plt.plot(t_span, y_i[:,0], label='Methane')
plt.plot(t_span, y_i[:,1], label='CO2')
plt.plot(t_span, y_i[:,2], label='Sulfur')
plt.plot(t_span, y_i[:,3], label='Water')
plt.legend()
plt.xlabel('Time [h]')
plt.ylabel('Molar Fraction [-]')
plt.title('Vapor Phase')

plt.figure()
plt.subplot(2,1,1)
plt.plot(t_span, h)
plt.title('Level of liquid in the digester')
plt.grid()
plt.xlabel('Time [h]')
plt.ylabel('Liquid Level [m]')
plt.subplot(2,1,2)
plt.plot(t_span, Q_in)
plt.title('Influent flow rate')
plt.grid()
plt.xlabel('Time [h]')
plt.ylabel('Flow rate [m3/h]')

plt.figure()
plt.plot(t_span, N_V[:,0], label='Methane')
plt.plot(t_span, N_V[:,1], label='CO2')

plt.show()

