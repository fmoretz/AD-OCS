''' Current version of AD_OCS model. Define the input according to the defined unit of measure from the spreadsheet.
    Define the timespan in hours.'''
import time 
import math
from pprint import pprint as pp
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*
from dataimport import*
from PhysConstants  import*

from functions import gompertz, growth_SRB, f_deviations, deviations_check, level_t, headspace_dynamics_discr, logistic_deviations
from models import AD_OCS, AMOCO_HN
from mix_real_gas import f_VL_realgas


# System definition
d_start = 1        # [h] - Start time
d_end   = 24*31     # [h] - End time
hours   = 0.1  # [h] - Discretization time

n_times = int((d_end-d_start)/hours)+1  # Number of time steps

print('***Intervals of {hours} hours. \n {n_times} time steps***'.format(hours=hours, n_times=n_times))

t_span = np.linspace(d_start,d_end,n_times) # time span
t_span_d = t_span/24 # time span in days

y_influent_changes = f_deviations(t_span_d, T3.index.values, y_in_0) # Get the deviated influent values at each timestamp
y_influent = y_influent_changes[:,0:5] # Remove the Q_in values

Q_in = y_influent_changes[:,5] # Remove the y_in values
gammaS_in = y_influent_changes[:,6] # Remove the y_in values
xW_in = y_influent_changes[:,7] # Remove the y_in values

# --------------------------------------------------------------------------------------------
# Ode Integration of AMOCO_HN: use to get X2 
# --------------------------------------------------------------------------------------------

y0= [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

YOUT_pre = odeint(AMOCO_HN, y0, t_span_d, hmax = (t_span_d[1]-t_span_d[0]), args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values))
# X2_pre = YOUT_pre[:,2]
S2_pre = YOUT_pre[:,5]
print('************** AMOCOHN OK *******************')

# --------------------------------------------------------------------------------------------
# Ode Integration of the AD_OCS_Model: uses the previous X2
# --------------------------------------------------------------------------------------------
y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions established from SS

YOUT1 = odeint(AD_OCS, y0, np.linspace(0,3,500), hmax = (t_span_d[1]-t_span_d[0]), args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values, t_span_d))

y0 = YOUT1[-1,:] # use the initial condition from the results of the previous integration

YOUT= odeint(AD_OCS, y0, t_span_d, hmax = (t_span_d[1]-t_span_d[0]), args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values, t_span_d))

# Get results
XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = YOUT[:,5]              # [mmol/L] - VFA dissolved
C  = YOUT[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

print('************** AD_OCS OK *******************')

# Solver Output: from all the variables from the ones of the ODE
mu1 = np.empty(len(XT))
mu2 = np.empty(len(XT))
CO2 = np.empty(len(XT))
B   = np.empty(len(XT))
phi = np.empty(len(XT))
p_C = np.empty(len(XT))
q_C = np.empty(len(XT))
q_M_u = np.empty(len(XT))
pH  = np.empty(len(XT))

for x in range(len(t_span)):
    mu1[x] = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # [1/d]      - Specific Growth Rate for X1 (Monod)
    mu2[x] = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # [1/d]      - Specific Growth Rate for X2 (Haldane)
    CO2[x] = C[x] + S2[x] - Z[x]                                 # [mmol/L]   - Dissolved CO2
    B[x]   = Z[x] - S2[x]                                        # [mmol/L]   - Alkalinity
    phi[x] = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH) # [atm]      - CO2 Partial Pressure
    q_C[x] = kLa*(CO2[x] - KH*p_C[x])                            # [mmol/L/d] - CO2 Outlet Molar Flow
    q_M_u[x] = k[5]*mu2[x]*X2[x]                                 # [mmol/L/d] - CH4 Outlet Molar Flow
    pH[x]  = np.real(-np.log10(Kb*CO2[x]/B[x]))                  # [-]        - System pH


# Sulfur and Oxygen Influence Evaluation

Xs     = np.zeros(len(t_span))                                   # [g/L]      - Sulfate Reducing Bacteria
Ss     = np.zeros(len(t_span))                                   # [mmol/L]      - Sulfur dissolved
y_S    = np.zeros(len(t_span))                                   # [-]        - Sulfur Mole fraction in gas phase
Ss_max = np.zeros(len(t_span))                                   # [mmol/L]      - Maximum Sulfur dissolved concentration
Xs_max = np.zeros(len(t_span))                                   # [g/L]      - Maximum Sulfate Reducing Bacteria concentration (Gompertz Asymptote)
mu_srb = np.zeros(len(t_span))                                   # [g/L/d]    - Gompertz parameter for SRB growth
growth_rate = np.zeros(len(t_span))                              # [g/L/d]    - Gompertz derivative: states the growth rate of SRB

dXsdt = np.zeros([len(t_span),len(t_span)])                     # [g/L/d]    - SRB growth rate matrix, preallocation
lam = 0                                                         # [-]        - Lag Time as defined by Gompertz

index = 0
gammaS_in_real = np.ones(len(gammaS_in))*gammaS_in[0]
xW_in_real = np.ones(len(gammaS_in))*xW_in[0]
loc_index = 0
k_logi = 0.1
for i in range(len(gammaS_in)):
    t=t_span[i]    
    if t < T3.index.values[0]*24:
        gammaS_in_real[i] = gammaS_in[i]
        xW_in_real[i]     = xW_in[i]
    elif t >= T3.index.values[-1]*24:
        gammaS_in_0_loc = gammaS_in_real[t_span<T3.index.values[-1]*24][-1]
        xW_in_0_loc = xW_in_real[t_span<T3.index.values[-1]*24][-1]
        gammaS_in_real[i] = logistic_deviations(t, T3.gammaSin.iloc[-1]*gammaS_in[0], k_logi, T3.index.values[-1]*24, gammaS_in_0_loc)
        xW_in_real[i] = logistic_deviations(t, T3.xWin.iloc[-1]*xW_in[0], k_logi, T3.index.values[-1]*24, xW_in_0_loc)
    else:
        gammaS_in_0_loc = gammaS_in_real[t_span<T3.index.values[loc_index]*24][-1]
        xW_in_0_loc = xW_in_real[t_span<T3.index.values[loc_index]*24][-1]
        gammaS_in_real[i] = logistic_deviations(t, T3.gammaSin.iloc[loc_index]*gammaS_in[0], k_logi, T3.index.values[loc_index]*24, gammaS_in_0_loc)
        xW_in_real[i] = logistic_deviations(t, T3.xWin.iloc[loc_index]*xW_in[0], k_logi, T3.index.values[loc_index]*24, xW_in_0_loc)
        if t >= T3.index.values[loc_index+1]*24:
            loc_index += 1
xW_in = np.ones(len(t_span))*xW_in_0
for j in range(len(t_span)):
    # Iterate over each time step        
    Ss_max[j] = gammaS_in_real[j]*1000/64*S2[j]                   #*y_influent[j,4]/y_influent[0,4]
    Xs_max[j] = Y_srb/(1-Y_srb)*Ss_max[j]     
    
    # Gompertz function for microbial population and dissolved sulfur
    mu_srb[j] = abs((- X2[0] + X2[j])/max(1e-16,t_span[j]-t_span[0]))      # [g/L/d]    - Gompertz parameter for SRB growth
    Xs[j]  = gompertz(t_span[j], Xs_max[j], mu_srb[j], lam)           # [g/L]      - Sulfate Reducing Bacteria - Gompertz
    Ss[j]  = Xs[j]*(1-Y_srb)/(Y_srb)                                  # [mmol/L]      - Sulfur dissolved concentration
    
    growth_rate[j] = mu2[j]*Xs[j]/X2[j]                      # [g/L/d]    - Get the growth rate of SRB at each time step as the maximum of the possible rates

    y_S[j]   = (KH_S*Ss[j])/P_dig                             # [-]        - Sulfur Mole fraction in gas phase (G/L equilibrium)

I_Ss = 1-KI_SS*Ss                                             # [-]   - Sulfate Ion concentration  
q_M  = q_M_u*I_Ss
q_tot = q_C + q_M                                             # [mmol/L/d] - Outlet global molar flow  
q_S = y_S*q_tot                                               # [mmol/L/d] - Sulfur Outlet Specific Molar Flow

### Assess the liquid volume dynamics ###
# Smooth the step change of the liquid influent
Q_in_real = np.ones(len(Q_in))*Q_in[0]
loc_index = 0
k_logi = 0.1
for i in range(len(Q_in)):
    t=t_span[i]    
    if t < T3.index.values[0]*24:
        Q_in_real[i] = Q_in[i]
    elif t >= T3.index.values[-1]*24:
        Q_in_0_loc = Q_in_real[t_span<T3.index.values[-1]*24][-1]
        Q_in_real[i] = logistic_deviations(t, T3.Qin.iloc[-1]*Q_in[0], k_logi, T3.index.values[-1]*24, Q_in_0_loc)
    else:
        Q_in_0_loc = Q_in_real[t_span<T3.index.values[loc_index]*24][-1]
        Q_in_real[i] = logistic_deviations(t, T3.Qin.iloc[loc_index]*Q_in[0], k_logi, T3.index.values[loc_index]*24, Q_in_0_loc)
        if t >= T3.index.values[loc_index+1]*24:
            loc_index += 1
Q_df = pd.DataFrame(index=t_span)
Q_df = pd.DataFrame()
Q_df["Qin"] = Q_in
Q_df["Qin_real"] = Q_in_real
Q_df["Qin_smooth"] = Q_df["Qin_real"].rolling(window=10).mean()
Q_df["Qin_smooth"].fillna(method='bfill', inplace=True)
Q_in_real = Q_df["Qin_smooth"].values

h = np.zeros(len(t_span))
h0 = h_SS
t_change = t_span[0]
index = 0
V_liq = np.zeros(len(t_span))                                 # [m3] - Liquid Volume
contr = 0
for i in range(len(t_span)):
    t = t_span[i]  
    if t < T3.index.values[0]*24:
        h[i] = h0
    elif t >= T3.index.values[-1]*24:
        if contr == 0:
            h0_loc = h[t_span<T3.index.values[index]*24][-1]
            h0_loc = level_t(t, D, Q_in_real[i], SR, h0_loc, t_change)  
            contr = 1
        t_change = T3.index.values[-1]*24
        h[i] = level_t(t, D, Q_in_real[i], SR, h0_loc, t_change)     
    else:    
        h0_loc = h[t_span<T3.index.values[index]*24][-1]
        t_change = T3.index.values[index]*24
        h[i] = level_t(t, D, Q_in_real[i], SR, h0_loc, t_change)     
        if t >= T3.index.values[index+1]*24:
            index += 1             
   
    
    if h[i] > hmax:
        print('!!! Level is too high !!!')
        input('Press Enter to continue')
    elif h[i] < hmin:
        print('!!! Level is too low !!!')
        input('Press Enter to continue')

h_df = pd.DataFrame(index=t_span_d)
h_df = pd.DataFrame()
h_df["h"] = h
h_df["h_smooth"] = h_df["h"].rolling(window=5).mean()
h_df["h_smooth"].fillna(method='bfill', inplace=True)
h = h_df["h_smooth"].values
V_liq = np.pi*Dr**2*h/4                             # [m3] - Liquid Volum
#### Evaluation of the G/L equilibrium effects on the system ####

n_species = 4                                                 # [-]       - Number of species in the system

F_i = np.zeros([len(t_span),n_species])                       # [-] - Outlet Molar Fraction - In of the flash
F_M = np.zeros(len(t_span))                                   # [-] - Methane Flow
F_C = np.zeros(len(t_span))                                   # [-] - Carbon Flow
F_S = np.zeros(len(t_span))                                   # [-] - Sulfur Flow
F_W = np.zeros(len(t_span))                                   # [-] - Water Flow

for i in range(len(t_span)):  
    F_W[i] = xW_in[i]*Q_in_real[i]*rho_W/18*1000/24             # [mol/h]   - Water Flow 
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

N_L_tot = np.zeros(len(t_span))                               # [mol/h] - Total Liquid Molar Flow
N_V_tot = np.zeros(len(t_span))                               # [mol/h] - Total Vapor Molar Flow

N_L = np.zeros([len(t_span),n_species])                       # [mol/h] - Liquid Phase Molar Flow
N_V = np.zeros([len(t_span),n_species])                       # [mol/h] - Vapor Phase Molar Flow
Q_L = np.zeros([len(t_span),n_species])                       # [m3/h] - Liquid Phase volumetric flow
Q_V = np.zeros([len(t_span),n_species])                       # [m3/h] - Vapor Phase volumetric flow

alpha_flash = np.zeros(len(t_span))                           # [-]     - 
MW_vett = [MW_M, MW_C, MW_S, MW_W]                            # [g/mol] - Molecular Weight Vector
rho_vett = [rho_M, rho_C, rho_S, rho_W]                       # [g/L]   - Density Vector
law = ('H', 'H', 'H', 'R')                                    # [-]     - Define the law to be used for the equilibrium calculation (R for Raoult's Law, H for Henry's Law)
for t in range(len(t_span)):  
    F_in = sum(F_i[t])                                        # [mol/h] - Total Inlet Molar Flow
    for i in range(n_species):        
        z_i[t,i] = F_i[t,i]/F_in                              # [-]     - Outlet Molar Fraction - In of the flash        
      
    VLsolution = f_VL_realgas(z_i[t], T, P_dig)               # [-]     - Solve the G/L equilibrium
    alpha_flash[t] = VLsolution[0]                            # [-] - Vapor fraction      
    K[t]= VLsolution[1]         
    
    N_V_tot[t] = F_in*alpha_flash[t]                                # [mol/h] - Vapor Molar Flow
    N_L_tot[t] = F_in - N_V_tot[t]                                  # [mol/h] - Liquid Molar Flow
    for i in range(n_species):
        x_i[t,i] = z_i[t,i]/(1+alpha_flash[t]*(K[t,i]-1))          # [-]     - Outlet Molar Fraction - Liquid phase of the flash
        y_i[t,i] = K[t,i]*x_i[t,i]                                 # [-]     - Outlet Molar Fraction - Vapor phase of the flash
        N_L[t,i] = N_L_tot[t]*x_i[t,i]                             # [mol/h] - Liquid Molar Flow of species i
        N_V[t,i] = N_V_tot[t]*y_i[t,i]                             # [mol/h] - Vapor Molar Flow of species i
        Q_L[t,i] = N_L[t,i]*MW_vett[i]/rho_vett[i]/1000            # [m3/h]  - Liquid Volumetric Flow of species i
        Q_V[t,i] = N_V[t,i]*(Rgas_m3_atm_K*(T+273.15)/P_norm)      # [m3/h]  - Vapor Volumetric Flow of species i


#### Evaluation of the Oxygen effect equilibrium effects on the system: CSTR behaviour ####
# IN: N_V[M, C, S, W] + N_in_O2
# OUT: N_out[M, C, S, W, O]


RC = 1.8                                                      # [gS/gO2]     - Stoichiometry of the reaction
alfa = 1                                                        # [-] - Power for sulfur concentration
beta = 0.25                                                      # [-] - Power for oxygen concentration
k_SOB = 1                                                     # [1/m3/h*(gS*gO)^proper exp] - Reaction rate constant
r_vett = np.zeros(len(t_span))                                  # [mg/L/min] - Reaction rate
S_vett = np.zeros((len(t_span),4))                              # [mg/L] - Sulfur concentration

V_gas = ((V_reactor - V_liq) + V_headspace)

tau_headspace = V_gas/Q_V.sum(1)                            # [h] - Time to fill the headspace

Q_in_air = 10 # [m3/h] - Air flow rate
Q_in_O2 = Q_in_air*0.21 # [m3/h] - Oxygen flow rate
Q_in_N2 = Q_in_air*0.79 # [m3/h] - Nitrogen flow rate
N_in_O2 = Q_in_O2/Rgas_m3_atm_K/(T_norm+273.15)*P_norm*np.ones(len(N_V)) # [mol/h] - Oxygen molar flow rate
N_in_N2 = Q_in_N2/Rgas_m3_atm_K/(T_norm+273.15)*P_norm*np.ones(len(N_V)) # [mol/h] - Oxygen molar flow rate

#N_in_O2 = 0.5*N_V[:,2] # [mol/h] - Oxygen Flow with respect to sulfur flow
headspace_dict = {}
digester_out = pd.DataFrame()

for ind in range(len(t_span)):
    t_cstr = np.arange(0,tau_headspace[ind],hours)
    while ind < (len(t_span)-len(t_cstr)):                
        N, r = headspace_dynamics_discr(N_V[ind:(ind+len(t_cstr)),:], N_in_O2[ind:(ind+len(t_cstr))], N_in_N2[ind:(ind+len(t_cstr))], P_dig, T, V_gas, t_cstr, k_SOB, RC, alfa, beta)
        headspace_dict[ind]= {'t_in': t_span[ind],'t_cstr': t_cstr, 'H2S': N['H2S'], 'H2O': N['H2O'], 'O2': N['O2'], 'SX': N['SX'], 'r_sob': r} # stores the results of the CSTR dynamics at each iteration, r [mg/L/h], N [mol/h]
        digester_out = digester_out.append(pd.DataFrame({
        't_in': t_span[ind], 't': (t_span[ind]+tau_headspace[ind]),
        'F_CH4': N_V[ind,0], 'F_CO2': N_V[ind,1], 'F_H2S': N['H2S'][-1], 'F_H2O': N['H2O'][-1], 'F_O2': N['O2'][-1], 'F_N2': N_in_N2[ind],
        'r_avg': (r.mean()*24)
        }, index=[ind]))
        break

digester_out['efficiency'] = (1-digester_out['F_H2S'] / N_V[0:len(digester_out),2])*100

# EFFECTS OF OXYGEN INHIBITION ON THE DIGESTER

x_V_O2 = digester_out['F_O2']/(digester_out['F_CH4']+digester_out['F_CO2']+digester_out['F_H2S']+digester_out['F_O2']+digester_out['F_H2O']+digester_out['F_N2'])
x_L_O2 = x_V_O2/H_O2_atm                                          # [-] - Mole fraction of oxygen in the liq
w_fr_O2 =  x_L_O2* P_dig/Rgas_L_atm_K/(T+273.15)*0.032         # [kg/m3] - Mass concentration of oxygen in the liq

I_O2 = 1/(1+w_fr_O2/KI_O2) # [-] - Inhibition factor of oxygen
Ch4_inh = digester_out['F_CH4']*I_O2

# Elaboration of the results: molar fractions, volumetric flows, etc.
digester_out['F_tot'] = digester_out['F_CH4']+digester_out['F_CO2']+digester_out['F_H2S']+digester_out['F_O2']+digester_out['F_H2O']+digester_out['F_N2']

# Create columns in digester_out accounting for molar fractions
digester_out['x_CH4'] = digester_out['F_CH4']/digester_out['F_tot']
digester_out['x_CO2'] = digester_out['F_CO2']/digester_out['F_tot']
digester_out['x_H2S'] = digester_out['F_H2S']/digester_out['F_tot']
digester_out['x_H2O'] = digester_out['F_H2O']/digester_out['F_tot']
digester_out['x_O2']  = digester_out['F_O2']/digester_out['F_tot']
digester_out['x_N2']  = digester_out['F_N2']/digester_out['F_tot']

# Create columns in digester_out for volumetric flowrates at given temperature and pressure

digester_out['Q_CH4'] = digester_out['F_CH4']*Rgas_m3_atm_K*(T_norm+273.15)/P_norm          # [m3/h] - Volumetric flowrate of CH4
digester_out['Q_CO2'] = digester_out['F_CO2']*Rgas_m3_atm_K*(T_norm+273.15)/P_norm          # [m3/h] - Volumetric flowrate of CO2
digester_out['Q_H2S'] = digester_out['F_H2S']*Rgas_m3_atm_K*(T_norm+273.15)/P_norm          # [m3/h] - Volumetric flowrate of H2S
digester_out['Q_H2O'] = digester_out['F_H2O']*Rgas_m3_atm_K*(T_norm+273.15)/P_norm          # [m3/h] - Volumetric flowrate of H2O
digester_out['Q_O2']  = digester_out['F_O2']*Rgas_m3_atm_K*(T_norm+273.15)/P_norm           # [m3/h] - Volumetric flowrate of O2
digester_out['Q_N2']  = digester_out['F_N2']*Rgas_m3_atm_K*(T_norm+273.15)/P_norm           # [m3/h] - Volumetric flowrate of N2
digester_out['Q_total'] = digester_out['Q_CH4']+digester_out['Q_CO2']+digester_out['Q_H2S']+digester_out['Q_H2O']+digester_out['Q_O2'] +digester_out['Q_N2']# [m3/h] - Total volumetric flowrate

digester_out['x_v_CH4'] = digester_out['Q_CH4']/digester_out['Q_total']          # [-] - Vol. Fraction of CH4 in the total flowrate
digester_out['x_v_CO2'] = digester_out['Q_CO2']/digester_out['Q_total']          # [-] - Vol. Fraction of CO2 in the total flowrate
digester_out['x_v_H2S'] = digester_out['Q_H2S']/digester_out['Q_total']          # [-] - Vol. Fraction of H2S in the total flowrate
digester_out['x_v_H2O'] = digester_out['Q_H2O']/digester_out['Q_total']          # [-] - Vol. Fraction of H2O in the total flowrate
digester_out['x_v_O2']  = digester_out['Q_O2']/digester_out['Q_total']           # [-] - Vol. Fraction of O2 in the total flowrate
digester_out['x_v_N2']  = digester_out['Q_N2']/digester_out['Q_total']           # [-] - Vol. Fraction of N2 in the total flowrate


# plt.figure(figsize=(20,10))
# plt.subplot(2,1,1)
# plt.plot(digester_out['t'], digester_out['x_v_CH4'], label='$CH_4$')
# plt.plot(digester_out['t'], digester_out['x_v_CO2'], label='$CO_2$')
# plt.plot(digester_out['t'], digester_out['x_v_H2O'], label='$H_2O$')
# plt.ylabel(ylabel='Volumetric fraction')
# plt.ylim([0,0.88])
# plt.legend()
# plt.xlim([digester_out['t'].min(), digester_out['t'].max()])
# plt.xticks(rotation=45)
# plt.subplot(2,1,2)
# plt.plot(digester_out['t'], digester_out['x_v_H2S']*1e+6, label='$H_2S$')
# plt.plot(digester_out['t'], digester_out['x_v_O2']*1e+6, label='$O_2$')
# plt.legend()
# plt.xlim([digester_out['t'].min(), digester_out['t'].max()])
# plt.xticks(rotation=45)
# plt.ylim([0,max(digester_out['x_v_H2S']*1e+6*1.1)])
# plt.ylabel('Volumetric fraction [ppm]')
# plt.suptitle('Outlet Fractions')
# plt.subplots_adjust(top=0.93)
# plt.show()

plt.figure(figsize=(20,10))
plt.subplot(2,1,1)
plt.plot(digester_out['t'], digester_out['x_v_CH4'], label='$CH_4$')
#plt.plot(digester_out['t'], digester_out['x_v_CH4'].rolling(10).mean(), label='$CH_4$ - Rolling mean')
# plt.plot(digester_out['t'], digester_out['x_v_CO2'], label='$CO_2$')
# plt.plot(digester_out['t'], digester_out['x_v_H2O'], label='$H_2O$')
plt.ylabel(ylabel='Volumetric fraction')

plt.legend()
plt.xlim([digester_out['t'].min(), digester_out['t'].max()])
plt.xticks(rotation=45)
plt.subplot(2,1,2)
plt.plot(digester_out['t'], digester_out['x_v_H2S']*1e+6, label='$H_2S$')
# plt.plot(digester_out['t'], digester_out['x_v_O2']*1e+6, label='$O_2$')
plt.legend()
plt.xlim([digester_out['t'].min(), digester_out['t'].max()])
plt.xticks(rotation=45)
plt.ylim([0,max(digester_out['x_v_H2S']*1e+6*1.1)])
plt.ylabel('Volumetric fraction [ppm]')
plt.suptitle('Outlet Fractions')
plt.subplots_adjust(top=0.93)
#plt.show()


# # plt.figure()
# # samples = [0, 250, 300, 350, 400, 450]
# # for i in samples:
# #     plt.plot(headspace_dict[i]['t_cstr'], headspace_dict[i]['H2S'], label = 't = %.4s h' % t_span[i])
# # plt.legend()


# plt.figure()
# plt.plot(t_span,growth_rate)
# plt.figure()
# plt.plot(t_span,S2-S2_pre)
# plt.show()


# plt.figure()
# plt.plot(digester_out['t'], digester_out['CH4']/(digester_out['CH4']+digester_out['CO2']), label='CH4')
# plt.plot(digester_out['t'], digester_out['CO2']/(digester_out['CH4']+digester_out['CO2']), label='CO2')
